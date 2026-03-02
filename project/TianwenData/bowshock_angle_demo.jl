"""
弓激波IMF法向夹角分析演示脚本

使用方法:
    julia bowshock_angle_demo.jl

功能:
    演示如何使用 bowshock_angle.jl 模块计算和可视化
    行星际磁场与弓激波法向的夹角
"""

# ============================================================
# 导入模块
# ============================================================

using Dates
using DataFrames
using Statistics
using LinearAlgebra

include("bowshock_angle.jl")
using .BowshockAngle

# ============================================================
# 配置参数
# ============================================================

const DATA_PATH = joinpath(@__DIR__, "..", "..", "Data", "32Hz")
const OUTPUT_PATH = joinpath(@__DIR__, "Picture", "BowshockAngle")
const DATE_STR = "20221025"
const THRESHOLD = 0.5  # 弓激波附近距离阈值 (Rm)

# ============================================================
# 1. 基础函数演示
# ============================================================

println("="^60)
println("1. 弓激波模型函数演示")
println("="^60)

# 测试弓激波半径计算
println("\n弓激波半径计算:")
test_x = [-8.0, -5.0, -2.0, 0.0]
for x in test_x
    r = BowshockAngle.bowshock_r(x)
    println("  x = $x Rm → r_bs = $r Rm")
end

# 测试法向量计算
println("\n弓激波法向量计算 (在 x=-5, y=1, z=0):")
n = BowshockAngle.calculate_bowshock_normal(-5.0, 1.0, 0.0)
println("  法向量 n = ($n)")
println("  |n| = $(norm(n))")

# 测试距离计算
println("\n到弓激波的距离计算:")
test_positions = [(-5.0, 1.0, 0.0), (-5.0, 3.0, 0.0), (-5.0, 2.0, 0.0)]
for pos in test_positions
    x, y, z = pos
    dist = BowshockAngle.distance_to_bowshock(x, y, z)
    println("  位置 $pos → 距离 = $(round(dist, digits=3)) Rm")
end

# 测试夹角计算
println("\nIMF与法向夹角计算:")
n_test = [1.0, 0.0, 0.0]  # 沿x轴
B_test1 = [1.0, 0.0, 0.0]  # 平行
B_test2 = [0.0, 1.0, 0.0]  # 垂直
B_test3 = [-1.0, 0.0, 0.0] # 反向

angle1 = BowshockAngle.calculate_imf_normal_angle(n_test, B_test1)
angle2 = BowshockAngle.calculate_imf_normal_angle(n_test, B_test2)
angle3 = BowshockAngle.calculate_imf_normal_angle(n_test, B_test3)

println("  n = [1,0,0], B = [1,0,0] → 夹角 = $(round(angle1, digits=1))°")
println("  n = [1,0,0], B = [0,1,0] → 夹角 = $(round(angle2, digits=1))°")
println("  n = [1,0,0], B = [-1,0,0] → 夹角 = $(round(angle3, digits=1))°")

# ============================================================
# 2. 数据加载与分析
# ============================================================

println("\n" * "="^60)
println("2. 数据加载与分析")
println("="^60)

# 加载数据
filename = "TW1_MOMAG_MSO_32Hz_" * DATE_STR * "_2C_v03.dat"
filepath = joinpath(DATA_PATH, filename)

println("\n加载数据文件: $filepath")
df = BowshockAngle.load_data(filepath)
println("✓ 加载完成: $(nrow(df)) 个数据点")
println("  时间范围: $(first(df.Time)) 至 $(last(df.Time))")

# 分析穿越事件
println("\n分析弓激波穿越事件 (阈值: $THRESHOLD Rm)...")
df, events = BowshockAngle.analyze_crossings(df; threshold=THRESHOLD)
println("✓ 发现 $(nrow(events)) 个穿越事件")

# 显示统计信息
if !isempty(events)
    println("\n穿越事件详情:")
    for i in 1:nrow(events)
        event = events[i, :]
        println("\n  事件 $i:")
        println("    时间段: $(event.start_time) - $(event.end_time)")
        println("    数据点: $(event.end_idx - event.start_idx + 1)")
        println("    最近距离: $(round(event.min_dist, digits=3)) Rm")
        println("    平均夹角: $(round(event.mean_angle, digits=1))°")
        println("    位置: ($(round(event.x_pos, digits=2)), " *
                "$(round(event.y_pos, digits=2)), " *
                "$(round(event.z_pos, digits=2))) Rm")

        # 计算该事件的夹角统计
        event_data = df[event.start_idx:event.end_idx, :]
        angles = filter(!isnan, event_data[!, :imf_normal_angle])
        if !isempty(angles)
            println("    夹角范围: [$(round(minimum(angles), digits=1))°, " *
                    "$(round(maximum(angles), digits=1))°]")
            println("    夹角标准差: $(round(std(angles), digits=1))°")
        end
    end
else
    println("\n未检测到穿越事件。可能原因:")
    println("  - 航天器距离弓激波较远")
    println("  - 阈值设置过小")
    println("  - 该时间段内无穿越事件")
end

# ============================================================
# 3. 可视化
# ============================================================

println("\n" * "="^60)
println("3. 可视化")
println("="^60)

# 绑制3D图
println("\n绑制3D分析图...")
fig_3d = BowshockAngle.plot_3d_analysis(df, events, OUTPUT_PATH)
println("✓ 3D图已保存到: $OUTPUT_PATH/bowshock_angle_3d.png")

# 绑制2D切面图
println("\n绑制2D切面图...")
fig_2d = BowshockAngle.plot_2d_slices(df, events, OUTPUT_PATH)
println("✓ 2D切面图已保存到: $OUTPUT_PATH/bowshock_angle_2d.png")

# 绑制时间序列图
println("\n绑制时间序列图...")
fig_ts = BowshockAngle.plot_angle_timeseries(df, OUTPUT_PATH)
println("✓ 时间序列图已保存到: $OUTPUT_PATH/angle_timeseries.png")

# ============================================================
# 4. 统计分析
# ============================================================

println("\n" * "="^60)
println("4. 统计分析")
println("="^60)

# 整体统计
valid_angles = filter(!isnan, df[!, :imf_normal_angle])
if !isempty(valid_angles)
    println("\n整体夹角统计:")
    println("  数据点数: $(length(valid_angles))")
    println("  平均夹角: $(round(mean(valid_angles), digits=1))°")
    println("  中位数: $(round(median(valid_angles), digits=1))°")
    println("  标准差: $(round(std(valid_angles), digits=1))°")
    println("  范围: [$(round(minimum(valid_angles), digits=1))°, " *
            "$(round(maximum(valid_angles), digits=1))°]")

    # 夹角分布
    println("\n  夹角分布:")
    bins = [0, 30, 60, 90, 120, 150, 180]
    for i in 1:length(bins)-1
        count = sum(bins[i] .<= valid_angles .< bins[i+1])
        percentage = 100 * count / length(valid_angles)
        println("    $(bins[i])°-$(bins[i+1])°: $count ($(round(percentage, digits=1))%)")
    end
end

# 距离统计
println("\n距离统计:")
valid_dist = filter(!isinf, df[!, :dist_to_bs])
if !isempty(valid_dist)
    println("  平均距离: $(round(mean(valid_dist), digits=3)) Rm")
    println("  最小距离: $(round(minimum(valid_dist), digits=3)) Rm")
    println("  最大距离: $(round(maximum(valid_dist), digits=3)) Rm")
end

# ============================================================
# 5. 运行完整分析流程
# ============================================================

println("\n" * "="^60)
println("5. 运行完整分析流程")
println("="^60)

# 运行另一个日期的分析作为演示
println("\n运行完整分析流程 (日期: 20221026)...")
df2, events2 = BowshockAngle.run_analysis("20221026", DATA_PATH, OUTPUT_PATH;
                                           threshold=THRESHOLD)

# ============================================================
# 6. 总结
# ============================================================

println("\n" * "="^60)
println("分析完成!")
println("="^60)

println("\n输出文件:")
output_files = filter(f -> endswith(f, ".png"), readdir(OUTPUT_PATH))
for f in output_files
    filepath = joinpath(OUTPUT_PATH, f)
    fsize = filesize(filepath) / 1024
    println("  - $f ($(round(fsize, digits=1)) KB)")
end

println("\n" * "="^60)
println("使用说明:")
println("="^60)
println("1. 修改 DATE_STR 变量可分析不同日期")
println("2. 调整 THRESHOLD 变量可改变事件检测灵敏度")
println("3. 输出图片位于: $OUTPUT_PATH")
println("4. 查看生成的PNG文件了解分析结果")
println("="^60)
