"""
Tianwen-1 弓激波IMF法向夹角分析模块

功能:
    - 计算弓激波表面法向量
    - 计算IMF与弓激波法向的夹角
    - 识别弓激波穿越事件
    - 可视化磁场方向与法向关系

使用示例:
    include("bowshock_angle.jl")
    using .BowshockAngle

    # 加载数据
    df = BowshockAngle.load_data(filepath)

    # 分析弓激波穿越
    crossings = BowshockAngle.analyze_crossings(df)

    # 绘制3D图
    BowshockAngle.plot_3d_analysis(df, crossings, "output.png")
"""

module BowshockAngle

using Dates
using DelimitedFiles
using DataFrames
using LinearAlgebra
using Statistics
using CairoMakie
using GeometryBasics
using LaTeXStrings

export bowshock, calculate_bowshock_normal, distance_to_bowshock
export calculate_imf_normal_angle, analyze_crossings
export plot_3d_analysis, plot_2d_slices

# ============================================================
# 常量定义
# ============================================================

const Rm = 3390.0  # 火星半径 (km)

# 弓激波模型参数 (Edberg et al., 2008)
const BS_XF = 0.55
const BS_EPSILON = 1.05
const BS_L = 2.10

# ============================================================
# 弓激波模型函数
# ============================================================

"""
    bowshock_r(x)

计算弓激波在给定x坐标处的半径

参数:
    x: MSO坐标系X坐标 (Rm)

返回值:
    弓激波半径 (Rm), Inf64 表示超出模型范围
"""
function bowshock_r(x)
    temp = (BS_EPSILON^2 - 1.0) * (x - BS_XF)^2 - 
           2 * BS_EPSILON * BS_L * (x - BS_XF) + BS_L^2
    return temp >= 0 ? sqrt(temp) : Inf64
end

"""
    bowshock_surface(x, y, z)

弓激波隐函数: F(x,y,z) = r^2 - (epsilon^2-1)(x-xF)^2 + 2*epsilon*L*(x-xF) - L^2 = 0

参数:
    x, y, z: MSO坐标 (Rm)

返回值:
    隐函数值 (0表示在表面上)
"""
function bowshock_surface(x, y, z)
    r = sqrt(y^2 + z^2)
    r_bs = bowshock_r(x)
    return isfinite(r_bs) ? r - r_bs : Inf64
end

"""
    calculate_bowshock_normal(x, y, z; dx=0.001)

计算弓激波表面法向量（通过数值微分）

参数:
    x, y, z: 位置坐标 (Rm)
    dx: 数值微分步长

返回值:
    单位法向量 (指向太阳方向为外法向)
"""
function calculate_bowshock_normal(x, y, z; dx::Float64=0.001)
    # 数值计算梯度
    # F(x,y,z) = sqrt(y^2+z^2) - r_bs(x)

    # 计算r_bs(x)的导数
    r_bs = bowshock_r(x)
    r_bs_plus = bowshock_r(x + dx)
    r_bs_minus = bowshock_r(x - dx)

    dr_bs_dx = (r_bs_plus - r_bs_minus) / (2 * dx)

    # 计算 sqrt(y^2+z^2) 的偏导数
    r_perp = sqrt(y^2 + z^2)
    if r_perp < 1e-10
        # 在x轴上，法向量沿x方向
        return [1.0, 0.0, 0.0]
    end

    # 梯度: ∇F = [-dr_bs_dx, y/r, z/r]
    grad = [-dr_bs_dx, y / r_perp, z / r_perp]

    # 归一化
    norm_grad = norm(grad)
    if norm_grad < 1e-10
        return [1.0, 0.0, 0.0]
    end

    n = grad / norm_grad

    # 确保法向量指向太阳方向（外法向）
    # 在弓激波处，外法向应大致指向-X方向
    if n[1] > 0
        n = -n
    end

    return n
end

"""
    distance_to_bowshock(x, y, z)

计算位置到弓激波表面的距离

参数:
    x, y, z: 位置坐标 (Rm)

返回值:
    距离 (Rm), 正值表示在弓激波外（太阳风侧），负值表示在内部
"""
function distance_to_bowshock(x, y, z)
    r = sqrt(y^2 + z^2)
    r_bs = bowshock_r(x)

    if !isfinite(r_bs)
        return Inf64
    end

    # 距离 = r - r_bs(x)
    # 正值: r > r_bs, 在弓激波外
    # 负值: r < r_bs, 在弓激波内
    return r - r_bs
end

"""
    calculate_imf_normal_angle(n, B)

计算IMF与弓激波法向的夹角

参数:
    n: 单位法向量
    B: IMF矢量 (nT)

返回值:
    夹角 (度, 0-180)
"""
function calculate_imf_normal_angle(n::Vector{Float64}, B::Vector{Float64})
    B_mag = norm(B)
    if B_mag < 1e-10
        return NaN
    end

    B_unit = B / B_mag

    # 点积: n · B = |n||B|cos(theta)
    cos_theta = dot(n, B_unit)

    # 限制在[-1, 1]范围内避免数值误差
    cos_theta = clamp(cos_theta, -1.0, 1.0)

    # 转换为角度
    theta = acosd(cos_theta)

    return theta
end

# ============================================================
# 数据加载与分析函数
# ============================================================

"""
    load_data(filepath; Rm_val=Rm)

加载32Hz磁场数据文件

参数:
    filepath: 数据文件路径
    Rm_val: 火星半径常数

返回值:
    DataFrame: 预处理后的数据
"""
function load_data(filepath::String; Rm_val::Real=Rm)
    df = DataFrame(readdlm(filepath, skipstart=19), :auto)
    rename!(df, ["Time", "Sampling_Rate", "X_MSO", "Y_MSO", "Z_MSO",
                "Probe_Position_X_MSO", "Probe_Position_Y_MSO", "Probe_Position_Z_MSO",
                "Roll", "Pitch", "Yaw", "Quality_Flags"])

    df[!, :Time] = map(x -> DateTime(x[begin:end-4], DateFormat("y-m-dTH:M:S.s")), df[!, :Time])
    df[!, [:Probe_Position_X_MSO, :Probe_Position_Y_MSO, :Probe_Position_Z_MSO]] ./= Rm_val

    # 计算IMF矢量大小
    df[!, :B_mag] = sqrt.(df[!, :X_MSO].^2 .+ df[!, :Y_MSO].^2 .+ df[!, :Z_MSO].^2)

    unique!(df)
    sort!(df)

    return df
end

"""
    analyze_crossings(df; threshold=0.5, window=300)

分析弓激波穿越事件

参数:
    df: 数据框
    threshold: 弓激波附近距离阈值 (Rm)
    window: 事件窗口大小（数据点数）

返回值:
    DataFrame: 穿越事件信息
"""
function analyze_crossings(df::DataFrame; threshold::Float64=0.5, window::Int=300)
    # 计算每个点到弓激波的距离
    distances = Float64[]
    angles = Float64[]
    normals = Vector{Float64}[]

    for i in 1:nrow(df)
        x = df[i, :Probe_Position_X_MSO]
        y = df[i, :Probe_Position_Y_MSO]
        z = df[i, :Probe_Position_Z_MSO]

        # 距离
        dist = distance_to_bowshock(x, y, z)
        push!(distances, dist)

        # 法向量
        n = calculate_bowshock_normal(x, y, z)
        push!(normals, n)

        # IMF与法向夹角
        B = [df[i, :X_MSO], df[i, :Y_MSO], df[i, :Z_MSO]]
        angle = calculate_imf_normal_angle(n, B)
        push!(angles, angle)
    end

    df[!, :dist_to_bs] = distances
    df[!, :imf_normal_angle] = angles
    df[!, :bs_normal_x] = [n[1] for n in normals]
    df[!, :bs_normal_y] = [n[2] for n in normals]
    df[!, :bs_normal_z] = [n[3] for n in normals]

    # 识别弓激波穿越事件（距离在阈值内）
    near_bs = abs.(distances) .< threshold

    # 找到连续的事件段
    events = DataFrame(
        start_idx = Int[],
        end_idx = Int[],
        start_time = DateTime[],
        end_time = DateTime[],
        min_dist = Float64[],
        mean_angle = Float64[],
        x_pos = Float64[],
        y_pos = Float64[],
        z_pos = Float64[]
    )

    in_event = false
    start_idx = 0

    for i in 1:length(near_bs)
        if near_bs[i] && !in_event
            # 事件开始
            in_event = true
            start_idx = i
        elseif !near_bs[i] && in_event
            # 事件结束
            in_event = false
            if i - start_idx > 10  # 至少10个数据点
                event_df = df[start_idx:i-1, :]
                min_dist_idx = argmin(abs.(event_df[!, :dist_to_bs]))

                push!(events, (
                    start_idx = start_idx,
                    end_idx = i-1,
                    start_time = event_df[1, :Time],
                    end_time = event_df[end, :Time],
                    min_dist = event_df[min_dist_idx, :dist_to_bs],
                    mean_angle = mean(filter(!isnan, event_df[!, :imf_normal_angle])),
                    x_pos = event_df[min_dist_idx, :Probe_Position_X_MSO],
                    y_pos = event_df[min_dist_idx, :Probe_Position_Y_MSO],
                    z_pos = event_df[min_dist_idx, :Probe_Position_Z_MSO]
                ))
            end
        end
    end

    return df, events
end

# ============================================================
# 可视化函数
# ============================================================

"""
    setup_theme()

设置绑图主题
"""
function setup_theme()
    CairoMakie.activate!()
    dark_latexfonts = merge(theme_dark(), theme_latexfonts())
    set_theme!(dark_latexfonts)
end

"""
    plot_3d_analysis(df, events, output_path; filename="bowshock_angle_3d.png")

绘制3D分析图

参数:
    df: 数据框
    events: 穿越事件信息
    output_path: 输出目录
    filename: 文件名
"""
function plot_3d_analysis(df::DataFrame, events::DataFrame, output_path::String;
                         filename::String="bowshock_angle_3d.png")
    setup_theme()

    fig = Figure(size=(16, 12))

    # 3D轴
    ax = Axis3(fig[1, 1],
               xlabel=L"$x$ ($R_{\mathrm{M}}$)",
               ylabel=L"$y$ ($R_{\mathrm{M}}$)",
               zlabel=L"$z$ ($R_{\mathrm{M}}$)",
               title="弓激波IMF法向夹角分析")

    # 生成弓激波表面网格
    x_range = range(-10, 3, length=50)
    phi_range = range(0, 2π, length=50)

    x_surf = Float64[]
    y_surf = Float64[]
    z_surf = Float64[]

    for x in x_range
        r_bs = bowshock_r(x)
        if isfinite(r_bs)
            for phi in phi_range
                push!(x_surf, x)
                push!(y_surf, r_bs * cos(phi))
                push!(z_surf, r_bs * sin(phi))
            end
        end
    end

    # 绑制弓激波表面（半透明）
    scatter!(ax, x_surf, y_surf, z_surf,
             color=:cyan, markersize=2, alpha=0.3, label="Bow Shock")

    # 绑制航天器轨迹
    lines!(ax, df[!, :Probe_Position_X_MSO],
           df[!, :Probe_Position_Y_MSO],
           df[!, :Probe_Position_Z_MSO],
           color=:white, linewidth=1, label="Trajectory")

    # 标记穿越事件
    if !isempty(events)
        for i in 1:min(nrow(events), 5)  # 最多显示5个事件
            event = events[i, :]

            # 事件位置
            pos = [event.x_pos, event.y_pos, event.z_pos]

            # 标记位置
            scatter!(ax, [pos[1]], [pos[2]], [pos[3]],
                    color=:red, markersize=15, label=i == 1 ? "Crossing" : nothing)

            # 计算该位置的法向量和IMF
            n = calculate_bowshock_normal(pos[1], pos[2], pos[3])

            # 找到最近的数据点的IMF
            idx = event.start_idx + div(event.end_idx - event.start_idx, 2)
            B = [df[idx, :X_MSO], df[idx, :Y_MSO], df[idx, :Z_MSO]]

            # 缩放因子
            scale = 1.0

            # 绑制法向量（绿色）
            arrows!(ax, [pos[1]], [pos[2]], [pos[3]],
                   [n[1]*scale], [n[2]*scale], [n[3]*scale],
                   color=:green, linewidth=3, arrowsize=0.3,
                   label=i == 1 ? "Normal" : nothing)

            # 绑制IMF（黄色）
            B_norm = B / (norm(B) + 1e-10) * scale
            arrows!(ax, [pos[1]], [pos[2]], [pos[3]],
                   [B_norm[1]], [B_norm[2]], [B_norm[3]],
                   color=:yellow, linewidth=3, arrowsize=0.3,
                   label=i == 1 ? "IMF" : nothing)
        end
    end

    # 火星
    theta = range(0, 2π, length=50)
    phi = range(0, π, length=25)
    for p in phi
        x_mars = cos(p)
        y_mars = sin(p) .* cos.(theta)
        z_mars = sin(p) .* sin.(theta)
        lines!(ax, x_mars, y_mars, z_mars, color=:white, linewidth=0.5, alpha=0.5)
    end

    # 设置坐标范围
    xlims!(ax, -8, 2)
    ylims!(ax, -4, 4)
    zlims!(ax, -4, 4)

    axislegend(ax, position=:rt)

    # 保存
    if !isdir(output_path)
        mkpath(output_path)
    end
    save(joinpath(output_path, filename), fig, pt_per_unit=1)

    return fig
end

"""
    plot_2d_slices(df, events, output_path; filename="bowshock_angle_2d.png")

绘制2D切面分析图

参数:
    df: 数据框
    events: 穿越事件信息
    output_path: 输出目录
    filename: 文件名
"""
function plot_2d_slices(df::DataFrame, events::DataFrame, output_path::String;
                       filename::String="bowshock_angle_2d.png")
    setup_theme()

    fig = Figure(size=(18, 6))

    # ========== X-R平面 ==========
    ax1 = Axis(fig[1, 1],
               xlabel=L"$x$ ($R_{\mathrm{M}}$)",
               ylabel=L"$r = \sqrt{y^2+z^2}$ ($R_{\mathrm{M}}$)",
               title="X-R平面",
               aspect=DataAspect())

    # 弓激波曲线
    x_range = range(-10, 3, length=200)
    r_bs_curve = [bowshock_r(x) for x in x_range]
    ind = findall(isfinite, r_bs_curve)
    lines!(ax1, x_range[ind], r_bs_curve[ind], color=:cyan, linewidth=2, label="Bow Shock")
    lines!(ax1, x_range[ind], -r_bs_curve[ind], color=:cyan, linewidth=2)

    # 火星
    theta = range(0, 2π, length=100)
    lines!(ax1, cos.(theta), sin.(theta), color=:white, linewidth=1, linestyle=:dash)

    # 轨迹
    r_traj = sqrt.(df[!, :Probe_Position_Y_MSO].^2 .+ df[!, :Probe_Position_Z_MSO].^2)
    lines!(ax1, df[!, :Probe_Position_X_MSO], r_traj,
           color=:white, linewidth=1, alpha=0.5)
    lines!(ax1, df[!, :Probe_Position_X_MSO], -r_traj,
           color=:white, linewidth=1, alpha=0.5)

    # 标记穿越事件
    if !isempty(events)
        for i in 1:min(nrow(events), 5)
            event = events[i, :]
            r_pos = sqrt(event.y_pos^2 + event.z_pos^2)

            # 位置
            scatter!(ax1, [event.x_pos], [r_pos], color=:red, markersize=12)
            scatter!(ax1, [event.x_pos], [-r_pos], color=:red, markersize=12)

            # 法向量和IMF（简化显示）
            n = calculate_bowshock_normal(event.x_pos, event.y_pos, event.z_pos)
            idx = event.start_idx + div(event.end_idx - event.start_idx, 2)
            B = [df[idx, :X_MSO], df[idx, :Y_MSO], df[idx, :Z_MSO]]

            # 在X-R平面的投影
            scale = 0.5
            arrows!(ax1, [event.x_pos], [r_pos],
                   [n[1]*scale], [n[2]*scale],
                   color=:green, linewidth=2, arrowsize=0.2)

            B_scale = B / (norm(B) + 1e-10) * scale
            r_B = sqrt(B[2]^2 + B[3]^2)
            arrows!(ax1, [event.x_pos], [r_pos],
                   [B_scale[1]], [r_B/(norm(B)+1e-10)*scale],
                   color=:yellow, linewidth=2, arrowsize=0.2)
        end
    end

    xlims!(ax1, -8, 2)
    ylims!(ax1, -5, 5)

    # ========== XY平面 ==========
    ax2 = Axis(fig[1, 2],
               xlabel=L"$x$ ($R_{\mathrm{M}}$)",
               ylabel=L"$y$ ($R_{\mathrm{M}}$)",
               title="XY平面 (z=0)",
               aspect=DataAspect())

    # 弓激波曲线（z=0截面）
    y_bs = [bowshock_r(x) for x in x_range]
    ind = findall(isfinite, y_bs)
    lines!(ax2, x_range[ind], y_bs[ind], color=:cyan, linewidth=2)
    lines!(ax2, x_range[ind], -y_bs[ind], color=:cyan, linewidth=2)

    # 火星
    lines!(ax2, cos.(theta), sin.(theta), color=:white, linewidth=1, linestyle=:dash)

    # 轨迹
    lines!(ax2, df[!, :Probe_Position_X_MSO], df[!, :Probe_Position_Y_MSO],
           color=:white, linewidth=1, alpha=0.5)

    # 标记事件
    if !isempty(events)
        for i in 1:min(nrow(events), 5)
            event = events[i, :]
            scatter!(ax2, [event.x_pos], [event.y_pos], color=:red, markersize=12)

            n = calculate_bowshock_normal(event.x_pos, event.y_pos, event.z_pos)
            idx = event.start_idx + div(event.end_idx - event.start_idx, 2)
            B = [df[idx, :X_MSO], df[idx, :Y_MSO], df[idx, :Z_MSO]]

            scale = 0.5
            arrows!(ax2, [event.x_pos], [event.y_pos],
                   [n[1]*scale], [n[2]*scale],
                   color=:green, linewidth=2, arrowsize=0.2)

            B_scale = B / (norm(B) + 1e-10) * scale
            arrows!(ax2, [event.x_pos], [event.y_pos],
                   [B_scale[1]], [B_scale[2]],
                   color=:yellow, linewidth=2, arrowsize=0.2)
        end
    end

    xlims!(ax2, -8, 2)
    ylims!(ax2, -5, 5)

    # ========== XZ平面 ==========
    ax3 = Axis(fig[1, 3],
               xlabel=L"$x$ ($R_{\mathrm{M}}$)",
               ylabel=L"$z$ ($R_{\mathrm{M}}$)",
               title="XZ平面 (y=0)",
               aspect=DataAspect())

    # 弓激波曲线（y=0截面）
    z_bs = [bowshock_r(x) for x in x_range]
    ind = findall(isfinite, z_bs)
    lines!(ax3, x_range[ind], z_bs[ind], color=:cyan, linewidth=2)
    lines!(ax3, x_range[ind], -z_bs[ind], color=:cyan, linewidth=2)

    # 火星
    lines!(ax3, cos.(theta), sin.(theta), color=:white, linewidth=1, linestyle=:dash)

    # 轨迹
    lines!(ax3, df[!, :Probe_Position_X_MSO], df[!, :Probe_Position_Z_MSO],
           color=:white, linewidth=1, alpha=0.5)

    # 标记事件
    if !isempty(events)
        for i in 1:min(nrow(events), 5)
            event = events[i, :]
            scatter!(ax3, [event.x_pos], [event.z_pos], color=:red, markersize=12)

            n = calculate_bowshock_normal(event.x_pos, event.y_pos, event.z_pos)
            idx = event.start_idx + div(event.end_idx - event.start_idx, 2)
            B = [df[idx, :X_MSO], df[idx, :Y_MSO], df[idx, :Z_MSO]]

            scale = 0.5
            arrows!(ax3, [event.x_pos], [event.z_pos],
                   [n[1]*scale], [n[3]*scale],
                   color=:green, linewidth=2, arrowsize=0.2)

            B_scale = B / (norm(B) + 1e-10) * scale
            arrows!(ax3, [event.x_pos], [event.z_pos],
                   [B_scale[1]], [B_scale[3]],
                   color=:yellow, linewidth=2, arrowsize=0.2)
        end
    end

    xlims!(ax3, -8, 2)
    ylims!(ax3, -5, 5)

    # 添加总标题
    Label(fig[0, :], "弓激波穿越事件分析 - IMF与法向夹角",
          fontsize=20, font=:bold)

    # 添加图例说明
    legend_elements = [
        (LineElement(color=:cyan, linewidth=2), "Bow Shock"),
        (LineElement(color=:white, linewidth=1), "Trajectory"),
        (MarkerElement(color=:red, marker=:circle, markersize=12), "Crossing"),
        (LineElement(color=:green, linewidth=2), "Normal"),
        (LineElement(color=:yellow, linewidth=2), "IMF")
    ]
    Legend(fig[2, 2], legend_elements, framevisible=true)

    colgap!(fig.layout, 20)

    # 保存
    if !isdir(output_path)
        mkpath(output_path)
    end
    save(joinpath(output_path, filename), fig, pt_per_unit=1)

    return fig
end

"""
    plot_angle_timeseries(df, output_path; filename="angle_timeseries.png")

绘制夹角时间序列图

参数:
    df: 数据框（包含夹角信息）
    output_path: 输出目录
    filename: 文件名
"""
function plot_angle_timeseries(df::DataFrame, output_path::String;
                              filename::String="angle_timeseries.png")
    setup_theme()

    fig = Figure(size=(16, 10))

    # 距离时间序列
    ax1 = Axis(fig[1, 1], ylabel=L"Distance ($R_{\mathrm{M}}$)")
    lines!(ax1, 1:nrow(df), df[!, :dist_to_bs], color=:cyan, linewidth=1)
    hlines!(ax1, [0], color=:white, linestyle=:dash)
    ax1.title = "到弓激波的距离"

    # 夹角时间序列
    ax2 = Axis(fig[2, 1], ylabel=L"Angle (degree)", xlabel="Time index")
    valid_idx = findall(!isnan, df[!, :imf_normal_angle])
    lines!(ax2, valid_idx, df[valid_idx, :imf_normal_angle], color=:orange, linewidth=1)
    ax2.title = "IMF与弓激波法向夹角"
    ylims!(ax2, 0, 180)

    # 磁场大小
    ax3 = Axis(fig[3, 1], ylabel=L"$|B|$ (nT)", xlabel="Time index")
    lines!(ax3, 1:nrow(df), df[!, :B_mag], color=:white, linewidth=1)
    ax3.title = "磁场强度"

    linkxaxes!(ax1, ax2, ax3)

    # 保存
    if !isdir(output_path)
        mkpath(output_path)
    end
    save(joinpath(output_path, filename), fig, pt_per_unit=1)

    return fig
end

# ============================================================
# 主函数
# ============================================================

"""
    run_analysis(datestr, datapath, outputpath; threshold=0.5)

运行完整的弓激波夹角分析

参数:
    datestr: 日期字符串
    datapath: 数据目录
    outputpath: 输出目录
    threshold: 弓激波附近距离阈值 (Rm)
"""
function run_analysis(datestr::String, datapath::String, outputpath::String;
                     threshold::Float64=0.5)
    # 加载数据
    filename = "TW1_MOMAG_MSO_32Hz_" * datestr * "_2C_v03.dat"
    filepath = joinpath(datapath, filename)

    println("="^60)
    println("弓激波IMF法向夹角分析")
    println("="^60)
    println("\n数据文件: $filepath")

    if !isfile(filepath)
        error("文件不存在: $filepath")
    end

    df = load_data(filepath)
    println("✓ 数据加载完成: $(nrow(df)) 个数据点")
    println("  时间范围: $(first(df.Time)) 至 $(last(df.Time))")

    # 分析穿越事件
    println("\n分析弓激波穿越事件...")
    df, events = analyze_crossings(df; threshold=threshold)
    println("✓ 发现 $(nrow(events)) 个穿越事件")

    # 显示事件统计
    if !isempty(events)
        println("\n穿越事件统计:")
        for i in 1:nrow(events)
            event = events[i, :]
            println("  事件 $i:")
            println("    时间: $(event.start_time) - $(event.end_time)")
            println("    最近距离: $(round(event.min_dist, digits=3)) Rm")
            println("    平均夹角: $(round(event.mean_angle, digits=1))°")
            println("    位置: ($(round(event.x_pos, digits=2)), " *
                    "$(round(event.y_pos, digits=2)), " *
                    "$(round(event.z_pos, digits=2))) Rm")
        end
    end

    # 绘制3D图
    println("\n绑制3D分析图...")
    plot_3d_analysis(df, events, outputpath)
    println("✓ 3D图已保存")

    # 绘制2D切面图
    println("\n绑制2D切面图...")
    plot_2d_slices(df, events, outputpath)
    println("✓ 2D切面图已保存")

    # 绘制时间序列
    println("\n绑制时间序列图...")
    plot_angle_timeseries(df, outputpath)
    println("✓ 时间序列图已保存")

    println("\n" * "="^60)
    println("分析完成!")
    println("="^60)
    println("输出目录: $outputpath")

    return df, events
end

end  # module
