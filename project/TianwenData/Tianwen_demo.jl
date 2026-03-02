"""
Tianwen-1 磁场数据分析脚本

本脚本演示如何使用 Tianwen_plot.jl 模块进行火星磁场数据分析。

使用方法:
    julia Tianwen_demo.jl

功能:
    1. 数据加载与预处理
    2. 小波分析
    3. 可视化绑图
    4. 完整分析流程

输出:
    - 绑图文件保存至指定输出目录
"""

# ============================================================
# 导入模块
# ============================================================

using Dates
using DelimitedFiles
using DataFrames
using CairoMakie
using LaTeXStrings

include("Tianwen_plot.jl")
using .TianwenPlot

# ============================================================
# 配置参数
# ============================================================

const DATA_PATH = joinpath(@__DIR__, "..", "..", "Data", "32Hz")
const OUTPUT_PATH = joinpath(@__DIR__, "Picture", "32Hz")
const DATE_STR = "20221025"
const NUM_FIGURES = 2

# ============================================================
# 1. 常量与模型演示
# ============================================================

println("="^60)
println("1. 常量与模型演示")
println("="^60)

println("火星半径 Rm = ", TianwenPlot.Rm, " km")

x_test = -5.0
r_shock = TianwenPlot.bowshock(x_test)
r_mp = TianwenPlot.magnetopause(x_test)
println("在 x = $x_test Rm 处:")
println("  弓激波半径 = $r_shock Rm")
println("  磁层顶半径 = $r_mp Rm")

model_curves = TianwenPlot.calculate_model_curves()
println("\n模型曲线计算完成")
println("  弓激波点数: ", length(model_curves.xshock))
println("  磁层顶点数: ", length(model_curves.xmp))

# ============================================================
# 2. 数据加载演示
# ============================================================

println("\n" * "="^60)
println("2. 数据加载演示")
println("="^60)

filename = "TW1_MOMAG_MSO_32Hz_" * DATE_STR * "_2C_v03.dat"
filepath = joinpath(DATA_PATH, filename)
println("读取数据文件: $filepath")
println("文件存在: ", isfile(filepath))

df = TianwenPlot.load_and_preprocess(filepath)
println("\n数据加载完成")
println("  数据形状: ", size(df))
println("  时间范围: ", first(df.Time), " 至 ", last(df.Time))

println("\n磁场统计 (nT):")
for col in [:X_MSO, :Y_MSO, :Z_MSO]
    println("  $col: mean=$(mean(df[!, col])), std=$(std(df[!, col]))")
end

# ============================================================
# 3. 小波分析演示
# ============================================================

println("\n" * "="^60)
println("3. 小波分析演示")
println("="^60)

BMSO32hz = Matrix{Float64}(df[:, [:X_MSO, :Y_MSO, :Z_MSO]])
magut32hz = df[!, :JulUT] .- df[1, :JulUT]
dt = 1.0/32.0

println("开始小波分析...")
wavelet_result = TianwenPlot.calculate_wavelet_power(BMSO32hz, dt)
println("小波分析完成!")
println("  功率矩阵: ", size(wavelet_result.total_power))
println("  周期范围: $(minimum(wavelet_result.period)) - $(maximum(wavelet_result.period)) s")

ds_result = TianwenPlot.downsample_data(wavelet_result.total_power, magut32hz, df; factor=32)
println("\n下采样完成")
println("  压缩后矩阵: ", size(ds_result.Bpowertemp))

# ============================================================
# 4. 绑图演示
# ============================================================

println("\n" * "="^60)
println("4. 绑图演示")
println("="^60)

TianwenPlot.setup_plot_theme()

println("创建轨迹面板...")
fig_traj = Figure(size=(18, 6))
traj_result = TianwenPlot.plot_trajectory_panels!(fig_traj, df, model_curves)
traj_result.ax1.title = DATE_STR * " - 轨迹"
TianwenPlot.setup_figure_layout!(fig_traj)
println("轨迹面板创建完成")

# 保存轨迹图
if !isdir(OUTPUT_PATH); mkpath(OUTPUT_PATH); end
TianwenPlot.save_figure(fig_traj, OUTPUT_PATH, "trajectory_$DATE_STR.png")
println("轨迹图已保存")

# ============================================================
# 5. 完整分析流程
# ============================================================

println("\n" * "="^60)
println("5. 完整分析流程")
println("="^60)

println("\n开始运行完整分析...")
println("  日期: $DATE_STR")
println("  输出路径: $OUTPUT_PATH")
println("  生成图片数: $NUM_FIGURES")

start_time_all = time()

results = TianwenPlot.run_analysis(
    DATE_STR,
    DATA_PATH,
    OUTPUT_PATH;
    date_offset=0,
    num_figures=NUM_FIGURES,
    hours_per_figure=4,
    figure_size=(30, 17)
)

elapsed_time = time() - start_time_all
println("\n分析完成!")
println("  总耗时: $(round(elapsed_time, digits=2)) 秒")
println("  生成图片数: $NUM_FIGURES")

# 列出生成的文件
output_files = filter(f -> endswith(f, ".png"), readdir(OUTPUT_PATH))
if !isempty(output_files)
    println("\n生成的文件:")
    for f in output_files
        fsize = filesize(joinpath(OUTPUT_PATH, f)) / 1024
        println("  - $f ($(round(fsize, digits=1)) KB)")
    end
end

# ============================================================
# 6. 自定义绑图示例
# ============================================================

println("\n" * "="^60)
println("6. 自定义绑图示例")
println("="^60)

dte = df[1, :Time]

println("创建自定义概览图...")
fig_custom = TianwenPlot.create_overview_figure(
    df,
    Matrix{Float64}(ds_result.BMSOtemp),
    ds_result.Bpowertemp,
    ds_result.maguttemp,
    wavelet_result.period,
    model_curves,
    DATE_STR,
    dte;
    figure_size=(25, 14),
    fontsize=20
)
TianwenPlot.save_figure(fig_custom, OUTPUT_PATH, "custom_overview_$DATE_STR.png")
println("自定义概览图已保存")

println("\n" * "="^60)
println("所有分析完成!")
println("="^60)
