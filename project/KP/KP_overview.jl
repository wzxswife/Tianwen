using Glob
using CommonDataFormat
using Dates
using Statistics
using CairoMakie
include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_plot.jl")
using .MAVEN_load
using .MAVEN_plot

# ========== 配置：日期范围 ==========
start_date = Date(2022, 08, 18)
end_date   = Date(2022, 08, 28)
dates = start_date:Dates.Day(1):end_date
# ===================================

# ========== 太阳风筛选配置 ==========
# false = 画全部数据；true = 只画太阳风区间
solar_wind_only = true

# 弓激波位置筛选（几何上排除磁鞘/磁层数据）
sw_use_position = true
# 卫星需在弓激波外多少 Rm 才视为太阳风（0.0 = 刚出弓激波就算）
sw_bowshock_margin = 0.2

# 物理参数阈值筛选（在位置筛选基础上，进一步排除非太阳风数据）
# 设 0.0 或 Inf 表示不启用该条件
sw_min_velocity = 300.0      # 最小 Vp (km/s)，太阳风通常 >300 km/s
sw_max_density  = 100.0       # 最大 Np (cm⁻³)，排除磁鞘高密度区
sw_min_tem      = 0.0        # 最小 Tp (eV)，0=不限制
sw_max_tem      = 50.0     # 最大 Tp (eV)，排除异常低温点
sw_max_mag      = 20.0       # 最大 |B| (nT)，排除磁鞘高磁场区
# ===================================

# ========== 波动时间标记（黑色虚线）==========
wave_times = DateTime[
    DateTime(2022,08,22,07,30,18)
]
# ============================================

data_dir = "E:\\Data"
output_dir = "E:\\Output"
MAVEN_path = joinpath(data_dir, "MAVEN")
MAVEN_output_path = joinpath(output_dir, "MAVEN")
KP_output_path = joinpath(MAVEN_output_path, "KP")

# 输出目录用起始日期
Overview_path = joinpath(KP_output_path, "Overview", 
    Dates.format(start_date, "yyyy"), Dates.format(start_date, "mm"))
mkpath(Overview_path)

# 遍历日期，加载并拼接 KP 数据
all_epoch = DateTime[]
all_HDensity = Float32[]
all_HVelox = Float32[]
all_HVeloy = Float32[]
all_HVeloz = Float32[]
all_HTem = Float32[]
all_magMSOX = Float32[]
all_magMSOY = Float32[]
all_magMSOZ = Float32[]
all_MSO_x = Float32[]
all_MSO_y = Float32[]
all_MSO_z = Float32[]

for date in dates
    datestr = Dates.format(date, "yyyymmdd")
    year_num = Dates.format(date, "yyyy")
    month_num = Dates.format(date, "mm")
    KP_path = joinpath(MAVEN_path, "KP", "l3", year_num, month_num)
    KP_file = glob("mvn_kp*$datestr*.jld2", KP_path)
    if isempty(KP_file)
        @warn "No KP file for $datestr, skipping"
        continue
    end
    println("Loading: ", KP_file[1])
    data = MAVEN_load.load_kp_l3(KP_file[1])
    nt = length(data[:time])
    append!(all_epoch, data[:time])
    append!(all_HDensity, data[:vars][:var_41])
    append!(all_HVelox, data[:vars][:var_43])
    append!(all_HVeloy, data[:vars][:var_45])
    append!(all_HVeloz, data[:vars][:var_47])
    append!(all_HTem, data[:vars][:var_49])
    append!(all_magMSOX, data[:vars][:var_128])
    append!(all_magMSOY, data[:vars][:var_130])
    append!(all_magMSOZ, data[:vars][:var_132])
    # MSO 位置（单位：km，需除以 Rm 归一化）
    if haskey(data[:vars], :var_190) && haskey(data[:vars], :var_191) && haskey(data[:vars], :var_192)
        append!(all_MSO_x, data[:vars][:var_190])
        append!(all_MSO_y, data[:vars][:var_191])
        append!(all_MSO_z, data[:vars][:var_192])
    else
        @warn "No MSO position vars in $datestr, solar wind filter won't work"
    end
end

if isempty(all_epoch)
    error("No data loaded for any date in the range!")
end

epoch = all_epoch
timeJUL = datetime2julian.(epoch)
HDensity = all_HDensity
velocity = sqrt.(all_HVelox.^2 .+ all_HVeloy.^2 .+ all_HVeloz.^2)
HTem = all_HTem
magMSOX = all_magMSOX
magMSOY = all_magMSOY
magMSOZ = all_magMSOZ
mag = sqrt.(magMSOX.^2 .+ magMSOY.^2 .+ magMSOZ.^2)

# ---------- 弓激波模型 ----------
function bowshock(xshock)
    xF = 0.55
    ϵ  = 1.05
    L  = 2.10
    temp = (ϵ^2-1.0)*(xshock-xF)^2 - 2ϵ*L*(xshock-xF) + L^2
    return temp >= 0 ? sqrt(temp) : Inf64
end

# ---------- 太阳风筛选 ----------
if solar_wind_only
    n_total = length(epoch)
    sw_mask = trues(n_total)

    # 第 1 层：弓激波位置筛选
    if sw_use_position
        if isempty(all_MSO_x)
            @warn "No MSO position data available, skipping position filter"
        else
            Rm = 3393.5  # 火星半径 km
            mso_x_Rm = all_MSO_x ./ Rm
            mso_y_Rm = all_MSO_y ./ Rm
            mso_z_Rm = all_MSO_z ./ Rm
            rho_sc = sqrt.(mso_y_Rm.^2 .+ mso_z_Rm.^2)     # 航天器距 X 轴距离
            rho_bs_raw = bowshock.(mso_x_Rm)                 # 弓激波半径（无余量）
            is_upstream = isinf.(rho_bs_raw)                  # Inf = 日出方向上游，无条件为太阳风
            rho_bs = rho_bs_raw .+ sw_bowshock_margin        # 含余量的弓激波半径
            # 非上游区域：X-R、X-Y、X-Z 三个方向都需在弓激波外
            mask_outside = (rho_sc .> rho_bs) .& (abs.(mso_y_Rm) .> rho_bs) .& (abs.(mso_z_Rm) .> rho_bs)
            sw_mask .&= is_upstream .| mask_outside
            println("  Position (3D bow shock + $(sw_bowshock_margin)Rm): $(count(sw_mask)) / $n_total remaining")
        end
    end

    # 第 2 层：物理参数阈值筛选
    if sw_min_velocity > 0
        sw_mask .&= velocity .>= sw_min_velocity
        # println("  Vp >= $(sw_min_velocity) km/s: $(count(sw_mask)) / $n_total remaining")
    end
    if sw_max_density < Inf
        sw_mask .&= HDensity .<= sw_max_density
        # println("  Np <= $(sw_max_density) cm⁻³: $(count(sw_mask)) / $n_total remaining")
    end
    if sw_min_tem > 0
        sw_mask .&= HTem .>= sw_min_tem
        # println("  Tp >= $(sw_min_tem) eV: $(count(sw_mask)) / $n_total remaining")
    end
    if sw_max_tem < Inf
        sw_mask .&= HTem .<= sw_max_tem
        # println("  Tp <= $(sw_max_tem) eV: $(count(sw_mask)) / $n_total remaining")
    end
    if sw_max_mag < Inf
        sw_mask .&= mag .<= sw_max_mag
        # println("  |B| <= $(sw_max_mag) nT: $(count(sw_mask)) / $n_total remaining")
    end

    epoch = epoch[sw_mask]
    timeJUL = timeJUL[sw_mask]
    HDensity = HDensity[sw_mask]
    velocity = velocity[sw_mask]
    HTem = HTem[sw_mask]
    magMSOX = magMSOX[sw_mask]
    magMSOY = magMSOY[sw_mask]
    magMSOZ = magMSOZ[sw_mask]
    mag = mag[sw_mask]
end

# 添加文件名后缀以区分
sw_suffix = solar_wind_only ? "_SW" : ""

# 检查筛选后是否还有数据
if isempty(timeJUL)
    error("No data remaining after filtering! Check your filter thresholds.")
end

# 构造多天的 x 轴刻度（每天正午 12:00）
xd = DateTime[]
for date in dates
    push!(xd, DateTime(date, Time(0, 0, 0)))
end
# 过滤超出数据范围的刻度
xmin, xmax = minimum(timeJUL), maximum(timeJUL)
xd = filter(x -> xmin <= datetime2julian(x) <= xmax, xd)
xtk = datetime2julian.(xd)
xtk_labels = Dates.format.(xd, "dd")

# 动态 xlabel：当月名称不跨月时显示 "Day of Nov 2021"
if Dates.month(start_date) == Dates.month(end_date)
    xlabel_str = "Day of $(Dates.format(start_date, "U")) $(Dates.format(start_date, "yyyy"))"
else
    xlabel_str = "Day of month ($(Dates.format(start_date, "yyyy")))"
end

# 波动时间转儒略日
wave_JUL = datetime2julian.(wave_times)

name = "KP_Overview$(sw_suffix)_" * Dates.format(start_date, "yyyymmdd") * 
       "_to_" * Dates.format(end_date, "yyyymmdd") * ".pdf"

CairoMakie.activate!()
size_inches = (20, 30)
size_pt = 72 .* size_inches
fig = Figure(size = size_pt, fontsize = 24, figure_padding = (5,20, 9, 10), rowgap = 20)
set_theme!(;
    Axis=(backgroundcolor = :white,
        bordercolor=:black, bordersize=15, spinewidth = 2, 
        xticksize=8, xminorticksize=4, yticksize=8, yminorticksize=4,
        xtickwidth = 2, ytickwidth = 2, xminortickwidth = 1, yminortickwidth = 1,
        xlabelsize=26, ylabelsize=26, xticklabelsize=22, yticklabelsize=22, 
        xminorticksvisible = true, yminorticksvisible = true,
        xminorticks = IntervalsBetween(10), yminorticks = IntervalsBetween(10))
)

ax1 = Axis(fig[1,1], xticks = (xtk, xtk_labels), ylabel=L"$B (nT)$")
ax2 = Axis(fig[2,1], xticks = (xtk, xtk_labels), 
    ylabel=L"$V_p$ ($km$ $s^{-1}$)")
    # ylabel=rich(rich("V"),subscript("p"),rich("(km s"),superscript("-1"),rich(")")))
ax3 = Axis(fig[3,1], xticks = (xtk, xtk_labels), 
    ylabel=L"$N_p$ ($cm^{-3}$)")
    # ylabel=rich(rich("N"),subscript("p"),rich("(cm"),superscript("-3"),rich(")")))
ax4 = Axis(fig[4,1], xticks = (xtk, xtk_labels), 
    xlabel=xlabel_str, ylabel=L"$T_p$ ($eV$)")
    # xlabel=xlabel_str, ylabel=rich(rich("T"),subscript("p"),rich("(eV)")))
ax1.xticklabelsvisible=false
ax2.xticklabelsvisible=false
ax3.xticklabelsvisible=false

xlims!(ax1, datetime2julian(epoch[1]), datetime2julian(epoch[end]))
ylims!(ax1, -10, 10)
ylims!(ax2, 320, 600)
ylims!(ax3, 0, 30)
ylims!(ax4, 0, 40)
linkxaxes!(ax1, ax2, ax3, ax4)

# 子图序号标注
text!(ax1, 0.02, 0.98, text="(a)", space=:relative, fontsize=30, align=(:left, :top))
text!(ax2, 0.02, 0.98, text="(b)", space=:relative, fontsize=30, align=(:left, :top))
text!(ax3, 0.02, 0.98, text="(c)", space=:relative, fontsize=30, align=(:left, :top))
text!(ax4, 0.02, 0.98, text="(d)", space=:relative, fontsize=30, align=(:left, :top))

scatter!(ax1, timeJUL, mag; markersize=2, color="black")
scatter!(ax1, timeJUL, magMSOX; markersize=2)
scatter!(ax1, timeJUL, magMSOY; markersize=2)
scatter!(ax1, timeJUL, magMSOZ; markersize=2)
# 图例：用文字颜色代替色块
lcn = Makie.wong_colors()
axislegend(ax1,
    [LineElement(color = (:white, 0.0)) for _ in 1:4],
    [rich(rich("B"), subscript("t"), rich("MSO"), color = "black"),
     rich(rich("B"), subscript("x"), rich("MSO"), color = lcn[1]),
     rich(rich("B"), subscript("y"), rich("MSO"), color = lcn[2]),
     rich(rich("B"), subscript("z"), rich("MSO"), color = lcn[3])];
    framevisible = false, position = :rt, labelsize=22, 
    orientation = :horizontal)
scatter!(ax2, timeJUL, velocity; markersize=2, color = "black")
scatter!(ax3, timeJUL, HDensity; markersize=2, color = "black")
scatter!(ax4, timeJUL, HTem; markersize=2, color = "black")

# 波动时间黑色虚线
if !isempty(wave_JUL)
    vlines!(ax1, wave_JUL; color=:black, linewidth=2, linestyle=:dash)
    vlines!(ax2, wave_JUL; color=:black, linewidth=2, linestyle=:dash)
    vlines!(ax3, wave_JUL; color=:black, linewidth=2, linestyle=:dash)
    vlines!(ax4, wave_JUL; color=:black, linewidth=2, linestyle=:dash)
end

save(joinpath(Overview_path, name), fig, pt_per_unit = 1)
println("Saved figure: ", joinpath(Overview_path, name))