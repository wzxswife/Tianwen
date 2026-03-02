"""
Tianwen-1 磁场数据分析模块

功能:
    - 火星磁场边界模型 (弓激波、磁层顶)
    - 32Hz磁场数据加载与预处理
    - 小波功率谱分析
    - 多面板可视化

使用示例:
    include("Tianwen_plot.jl")
    using .TianwenPlot

    # 一键加载和分析
    df = TianwenPlot.load_data("data.dat")
    result = TianwenPlot.analyze(df)
    TianwenPlot.plot_all(result, "output.png")
"""

module TianwenPlot

using Dates
using DelimitedFiles
using DataFrames
using Interpolations
using LinearAlgebra
using DSP
using Wavelets
using CairoMakie
using GeometryBasics
using LaTeXStrings

export bowshock, magnetopause, calculate_model_curves
export load_data, analyze, plot_overview, plot_trajectory
export run_analysis

# ============================================================
# 常量定义
# ============================================================

const Rm = 3390.0  # 火星半径 (km)

# ============================================================
# 模型函数
# ============================================================

"""
    bowshock(x)

火星弓激波模型 (Edberg et al., 2008)

参数:
    x: 火星固连坐标系X坐标 (Rm)

返回值:
    弓激波半径 (Rm), Inf64 表示超出模型范围
"""
function bowshock(xshock)
    xF = 0.55
    ϵ = 1.05
    L = 2.10
    temp = (ϵ^2-1.0)*(xshock-xF)^2 - 2ϵ*L*(xshock-xF) + L^2
    return temp >= 0 ? sqrt(temp) : Inf64
end

"""
    magnetopause(x)

火星磁层顶模型 (Edberg et al., 2008)

参数:
    x: 火星固连坐标系X坐标 (Rm)

返回值:
    磁层顶半径 (Rm), Inf64 表示超出模型范围
"""
function magnetopause(xmp)
    rSD = 1.33
    xF = 0.86
    ϵ = 0.92
    L = 0.90
    temp = (ϵ^2-1.0)*(xmp-xF)^2 - 2ϵ*L*(xmp-xF) + L^2
    return temp >= 0 ? sqrt(temp) : Inf64
end

"""
    calculate_model_curves(; np=5001, xmin=-10, xmax=3)

计算火星磁场边界模型曲线

关键字参数:
    np: 曲线分辨率
    xmin: X轴范围下限
    xmax: X轴范围上限

返回值:
    NamedTuple: (xshock, ryzshock, xmp, ryzmp)
"""
function calculate_model_curves(; np::Integer=5001, xmin::Real=-10, xmax::Real=3)
    # 弓激波
    xsh = range(xmin, xmax, length=np)
    ryzsh = bowshock.(xsh)
    ind = findall(isfinite, ryzsh)
    xshock = [xsh[ind]; reverse(xsh[ind])]
    ryzshock = [ryzsh[ind]; -reverse(ryzsh[ind])]

    # 磁层顶
    xmp_range = range(xmin, xmax, length=np)
    ryzmp_val = magnetopause.(xmp_range)
    ind = findall(isfinite, ryzmp_val)
    xmp = [xmp_range[ind]; reverse(xmp_range[ind])]
    ryzmp = [ryzmp_val[ind]; -reverse(ryzmp_val[ind])]

    return (xshock=xshock, ryzshock=ryzshock, xmp=xmp, ryzmp=ryzmp)
end

# ============================================================
# 数据加载与分析函数
# ============================================================

"""
    load_data(filepath; Rm_val=Rm)

一键加载和预处理32Hz磁场数据文件

参数:
    filepath: 数据文件路径
    Rm_val: 火星半径常数 (用于位置归一化)

返回值:
    DataFrame: 预处理后的数据框
"""
function load_data(filepath::String; Rm_val::Real=Rm)
    # 读取原始数据
    df = DataFrame(readdlm(filepath, skipstart=19), :auto)
    rename!(df, ["Time", "Sampling_Rate", "X_MSO", "Y_MSO", "Z_MSO",
                "Probe_Position_X_MSO", "Probe_Position_Y_MSO", "Probe_Position_Z_MSO",
                "Roll", "Pitch", "Yaw", "Quality_Flags"])

    # 时间转换
    df[!, :Time] = map(x -> DateTime(x[begin:end-4], DateFormat("y-m-dTH:M:S.s")), df[!, :Time])
    df[!, :JulUT] = datetime2julian.(df[!, :Time])

    # 位置归一化
    df[!, [:Probe_Position_X_MSO, :Probe_Position_Y_MSO, :Probe_Position_Z_MSO]] ./= Rm_val

    # 去重排序
    unique!(df)
    sort!(df)

    # 插值缺失值
    if any(isnan, df[!, :X_MSO])
        col_data = df[:, [:X_MSO, :Y_MSO, :Z_MSO]]
        time_data = df[:, :JulUT]
        ind0 = findall(!isnan, col_data[:, 1])
        ind1 = findall(isnan, col_data[:, 1])

        if length(ind1) > 0 && length(ind0) > 1
            for col in [:X_MSO, :Y_MSO, :Z_MSO]
                interp = linear_interpolation(time_data[ind0], col_data[ind0, col]; extrapolation_bc=Line())
                col_data[ind1, col] = interp(time_data[ind1])
            end
        end
    end

    return df
end

"""
    analyze(df; dt=1.0/32.0, mother="MORLET")

对磁场数据进行小波分析

参数:
    df: 预处理后的数据框
    dt: 采样间隔 (秒)
    mother: 母小波类型

返回值:
    NamedTuple: 包含功率谱、时间等分析结果
"""
function analyze(df::DataFrame; dt::Float64=1.0/32.0, mother::String="MORLET")
    # 包含小波模块
    include(joinpath(@__DIR__, "..", "src", "scripts", "TCWavelet.jl"))

    # 提取数据
    BMSO = Matrix{Float64}(df[:, [:X_MSO, :Y_MSO, :Z_MSO]])
    magut = df[!, :JulUT] .- df[1, :JulUT]
    ns = size(BMSO, 1)

    # 计算三分量小波功率
    wave, period, scale, coi = wavelet(reshape(BMSO[:, 1], ns), dt; pad=1, mother=mother)
    xpower = abs.(wave).^2

    wave, period, scale, coi = wavelet(reshape(BMSO[:, 2], ns), dt; pad=1, mother=mother)
    ypower = abs.(wave).^2

    wave, period, scale, coi = wavelet(reshape(BMSO[:, 3], ns), dt; pad=1, mother=mother)
    zpower = abs.(wave).^2

    total_power = xpower + ypower + zpower

    # 下采样 (32Hz -> 1Hz)
    factor = 32
    len = round(Int, size(magut, 1) / factor)
    Bpowertemp = zeros(size(total_power, 1), len)
    maguttemp = zeros(len)

    for i in 1:len
        Bpowertemp[:, i] = total_power[:, (i-1)*factor+1]
        maguttemp[i] = magut[(i-1)*factor+1]
    end

    magtemp = df[1:factor:(len-1)*factor+1, :]
    BMSOtemp = Matrix{Float64}(df[1:factor:(len-1)*factor+1, [:X_MSO, :Y_MSO, :Z_MSO]])

    return (
        df=df, BMSO=BMSO, magut=magut,
        period=period, scale=scale, coi=coi,
        total_power=total_power,
        Bpowertemp=Bpowertemp, maguttemp=maguttemp,
        magtemp=magtemp, BMSOtemp=BMSOtemp
    )
end

# ============================================================
# 绑图函数
# ============================================================

"""
    setup_theme()

设置绑图主题 (深色+LaTeX字体)
"""
function setup_theme()
    CairoMakie.activate!()
    dark_latexfonts = merge(theme_dark(), theme_latexfonts())
    set_theme!(dark_latexfonts)
end

"""
    plot_trajectory!(ax, df, model_curves; kwargs...)

在轴上绑制航天器轨迹

参数:
    ax: Makie轴对象
    df: 数据框
    model_curves: 模型曲线

关键字参数:
    xlims: X轴范围
    ylims: Y轴范围
"""
function plot_trajectory!(ax, df::DataFrame, model_curves;
                         xlims=(-10, 10), ylims=(-5, 5))
    xshock = model_curves.xshock
    ryzshock = model_curves.ryzshock
    xmp = model_curves.xmp
    ryzmp = model_curves.ryzmp

    # 轨迹 (XZ投影)
    lines!(ax, df[!, :Probe_Position_X_MSO], df[!, :Probe_Position_Z_MSO],
           color=:white, linewidth=1.5)

    # 边界模型
    lines!(ax, xshock, ryzshock, color=:cyan, linewidth=2, label="Bow Shock")
    lines!(ax, xmp, ryzmp, color=:orange, linewidth=2, label="Magnetopause")

    # 火星表面
    theta = range(0, 2π, length=100)
    lines!(ax, cos.(theta), sin.(theta), color=:white, linewidth=1, linestyle=:dash)

    # 坐标范围
    xlims!(ax, xlims...)
    ylims!(ax, ylims...)

    ax.xlabel = L"$x$ ($R_{\mathrm{M}}$)"
    ax.ylabel = L"$z$ ($R_{\mathrm{M}}$)"
    axislegend(ax)

    return ax
end

"""
    plot_magnetic_field!(ax, df; kwargs...)

绑制磁场时间序列

参数:
    ax: Makie轴对象
    df: 数据框

关键字参数:
    linewidth: 线宽
"""
function plot_magnetic_field!(ax, df::DataFrame; linewidth=1.0)
    time = df[!, :JulUT] .- df[1, :JulUT]

    lines!(ax, time, df[!, :X_MSO], color=:red, linewidth=linewidth, label=L"$B_x$")
    lines!(ax, time, df[!, :Y_MSO], color=:green, linewidth=linewidth, label=L"$B_y$")
    lines!(ax, time, df[!, :Z_MSO], color=:blue, linewidth=linewidth, label=L"$B_z$")

    ax.xlabel = "Time (s)"
    ax.ylabel = L"$\mathbf{B}$ (nT)"
    axislegend(ax)

    return ax
end

"""
    plot_wavelet_power!(ax, result; kwargs...)

绑制小波功率谱

参数:
    ax: Makie轴对象
    result: analyze() 返回的结果
"""
function plot_wavelet_power!(ax, result; colorrange=(2e-2, 3e2))
    ax.yscale = log2
    ax.yreversed = true
    ylims!(ax, 1.0/16.0, 128)

    heatmap!(ax, result.maguttemp, result.period, result.Bpowertemp',
             colorscale=log10, colormap=:gist_earth, colorrange=colorrange)

    ax.xlabel = "Time (s)"
    ax.ylabel = L"Period $T$ (s)"

    return ax
end

"""
    plot_overview(result, model_curves, datestr, output_path; kwargs...)

创建并保存完整概览图

参数:
    result: analyze() 返回的结果
    model_curves: 模型曲线
    datestr: 日期字符串
    output_path: 输出目录

关键字参数:
    figure_size: 图片尺寸
    fontsize: 字体大小
"""
function plot_overview(result, model_curves, datestr, output_path;
                      figure_size=(30, 17), fontsize=25)
    setup_theme()
    size_pt = 72 .* figure_size

    df = result.magtemp
    dte = df[1, :Time]

    fig = Figure(size=figure_size, fontsize=fontsize)

    # 轨迹面板
    ax1 = Axis(fig[1, 1], aspect=DataAspect(),
               xlabel=L"$x$ ($R_{\mathrm{M}}$)", ylabel=L"$z$ ($R_{\mathrm{M}}$)",
               title=datestr)
    xlims!(ax1, -10, 10)
    ylims!(ax1, -5, 5)

    theta = range(0, 2π, length=100)
    lines!(ax1, cos.(theta), sin.(theta), color=:white, linewidth=1, linestyle=:dash)
    lines!(ax1, model_curves.xshock, model_curves.ryzshock, color=:cyan, linewidth=2)
    lines!(ax1, model_curves.xmp, model_curves.ryzmp, color=:orange, linewidth=2)
    lines!(ax1, df[!, :Probe_Position_X_MSO], df[!, :Probe_Position_Z_MSO],
           color=:white, linewidth=1)

    # 磁场时间序列
    time = df[!, :JulUT] .- df[1, :JulUT]
    ax2 = Axis(fig[2, 1], xlabel="Time (s)", ylabel=L"$\mathbf{B}$ (nT)")
    lines!(ax2, time, df[!, :X_MSO], color=:red, linewidth=0.8, label=L"$B_x$")
    lines!(ax2, time, df[!, :Y_MSO], color=:green, linewidth=0.8, label=L"$B_y$")
    lines!(ax2, time, df[!, :Z_MSO], color=:blue, linewidth=0.8, label=L"$B_z$")
    axislegend(ax2, position=:rt)

    # 小波功率谱
    ax3 = Axis(fig[3, 1], xlabel="Time (s)", ylabel=L"Period $T$ (s)", yscale=log2)
    ax3.yreversed = true
    ylims!(ax3, 1.0/16.0, 128)
    heatmap!(ax3, result.maguttemp, result.period, result.Bpowertemp',
             colorscale=log10, colormap=:gist_earth, colorrange=(2e-2, 3e2))
    Colorbar(fig[3, 2], label=L"Wavelet Power $P_{\mathrm{B}}$")

    # 保存
    if !isdir(output_path); mkpath(output_path); end
    filename = joinpath(output_path, "overview_$(datestr).png")
    save(filename, fig, pt_per_unit=1)

    return fig, filename
end

# ============================================================
# 主函数
# ============================================================

"""
    run_analysis(datestr, datapath, outputpath; kwargs...)

运行完整的分析流程

参数:
    datestr: 日期字符串 (yyyymmdd格式)
    datapath: 数据目录路径
    outputpath: 输出目录路径

关键字参数:
    date_offset: 日期偏移
    num_figures: 生成图片数量
    hours_per_figure: 每张图覆盖的小时数
    figure_size: 图片尺寸
"""
function run_analysis(datestr::String, datapath::String, outputpath::String;
                     date_offset::Integer=0,
                     num_figures::Integer=6,
                     hours_per_figure::Integer=4,
                     figure_size::Tuple{<:Real,<:Real}=(30, 17))
    # 模型曲线
    model_curves = calculate_model_curves()

    # 数据文件
    filename = "TW1_MOMAG_MSO_32Hz_" * datestr * "_2C_v03.dat"
    filepath = joinpath(datapath, filename)
    println("读取: $filepath")

    # 加载数据
    df = load_data(filepath)

    # 分析
    result = analyze(df)

    # 绑图
    setup_theme()

    date = DateTime(2022, 10, 24) .+ Dates.Day(1) .* range(0, 10)
    dte = date[date_offset + 1]

    for hr in 1:num_figures
        println("绑图 $hr / $num_figures")
        start_time = dte + Dates.Hour((hr-1)*hours_per_figure)

        fig = Figure(size=figure_size, fontsize=25)

        # 轨迹
        ax1 = Axis(fig[1, 1], aspect=DataAspect(),
                   xlabel=L"$x$ ($R_{\mathrm{M}}$)", ylabel=L"$z$ ($R_{\mathrm{M}}$)")
        xlims!(ax1, -10, 10)
        ylims!(ax1, -5, 5)
        theta = range(0, 2π, length=100)
        lines!(ax1, cos.(theta), sin.(theta), color=:white, linewidth=1, linestyle=:dash)
        lines!(ax1, model_curves.xshock, model_curves.ryzshock, color=:cyan, linewidth=2)
        lines!(ax1, model_curves.xmp, model_curves.ryzmp, color=:orange, linewidth=2)
        lines!(ax1, result.magtemp[!, :Probe_Position_X_MSO],
               result.magtemp[!, :Probe_Position_Z_MSO], color=:white, linewidth=1)

        # 磁场
        ax2 = Axis(fig[2, 1], xlabel="Time (s)", ylabel=L"$\mathbf{B}$ (nT)")
        time = result.magtemp[!, :JulUT] .- result.magtemp[1, :JulUT]
        lines!(ax2, time, result.magtemp[!, :X_MSO], color=:red, linewidth=0.8, label=L"$B_x$")
        lines!(ax2, time, result.magtemp[!, :Y_MSO], color=:green, linewidth=0.8, label=L"$B_y$")
        lines!(ax2, time, result.magtemp[!, :Z_MSO], color=:blue, linewidth=0.8, label=L"$B_z$")
        axislegend(ax2, position=:rt)

        # 小波
        ax3 = Axis(fig[3, 1], xlabel="Time (s)", ylabel=L"Period $T$ (s)", yscale=log2)
        ax3.yreversed = true
        ylims!(ax3, 1.0/16.0, 128)
        heatmap!(ax3, result.maguttemp, result.period, result.Bpowertemp',
                 colorscale=log10, colormap=:gist_earth, colorrange=(2e-2, 3e2))
        Colorbar(fig[3, 2], label=L"Wavelet Power $P_{\mathrm{B}}$")

        resize_to_layout!(fig)

        # 保存
        if !isdir(outputpath); mkpath(outputpath); end
        save(joinpath(outputpath, "$(datestr)$hr.png"), fig, pt_per_unit=1)
    end

    return result
end

end  # module
