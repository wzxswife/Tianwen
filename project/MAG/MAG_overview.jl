using Glob
using CommonDataFormat
using Dates
using Statistics
using CairoMakie
include("../../src/scripts/MAG_plot.jl")
include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_plot.jl")
using .MAG_plot
using .MAVEN_load
using .MAVEN_plot

date = Date(2022, 8, 24)
datestr = Dates.format.(date, "yyyymmdd")
year_num = Dates.format(date, "yyyy")
month_num = Dates.format(date, "mm")

data_dir = "E:\\Data"
output_dir = "E:\\Output"
MAVEN_path = joinpath(data_dir, "MAVEN")
MAG_path = joinpath(MAVEN_path, "MAG", "l2", year_num, month_num)
KP_path = joinpath(MAVEN_path, "KP", "l3", year_num, month_num)
SWIA_path = joinpath(MAVEN_path, "SWIA", "l2", year_num, month_num)
STATIC_path = joinpath(MAVEN_path, "STATIC", "l2", year_num, month_num)
MAVEN_output_path = joinpath(output_dir, "MAVEN")
STATIC_output_path = joinpath(MAVEN_output_path, "STATIC")
KP_output_path = joinpath(MAVEN_output_path, "KP")
MAG_output_path = joinpath(MAVEN_output_path, "MAG")
Overview_path = joinpath(MAG_output_path, "Overview", year_num, month_num)
mkpath(Overview_path)
MAG_file = glob("mvn_mag_l2*$datestr*.sts", MAG_path)
KP_file = glob("mvn_kp*$datestr*.jld2", KP_path)
SWIA_file = glob("mvn_swi_l2_coarsesvy3d*$datestr*.cdf", SWIA_path)
STA_file = glob("mvn_sta_l2*32e4d16a8m*$datestr*.cdf", STATIC_path)
println("MAG file: ", MAG_file)
println("KP file: ", KP_file)
println("SWIA file: ", SWIA_file)
println("STATIC file: ", STA_file)
println("Output path: ", Overview_path)

mag_data = load_mag_l2(MAG_file[1])
println("Data Loaded.")
const Rm = 3389.5
mag_data[:position] = mag_data[:position]./Rm

time_range = DateTime(date, Time(7, 24, 00)).+Dates.Minute(2) .* range(0, 6)

CairoMakie.activate!()
size_inches = (20, 18)
size_pt = 72 .* size_inches
fig = Figure(size = size_pt, fontsize = 24, figure_padding = (30, 10, 20, 20), rowgap = 20)
set_theme!(;
    Axis=(backgroundcolor = :white,
        bordercolor=:black, bordersize=15, spinewidth = 2, 
        xlabelsize=26, ylabelsize=26, xticklabelsize=22, yticklabelsize=22, 
        xminorticksvisible = true, yminorticksvisible = true,
        xminorticks = IntervalsBetween(10), yminorticks = IntervalsBetween(10))
)

xlim1 = (6,-12)
ax1 = Axis(fig[1,1:3], aspect=DataAspect(), 
    xlabel=L"$x$ ($R_{\mathrm{M}}$)", 
    ylabel=L"$\sqrt{y^2+z^2}$ ($R_{\mathrm{M}}$)")
MAG_plot.bowshock_plot_xr(ax1;xlim=xlim1,ylim=(0,5))
MAG_plot.plot_spacecraft_orbit_xr(ax1,mag_data,[time_range[1]];xlim=xlim1,ylim=(0,5))

ax2 = Axis(fig[2, 1:3], ylabel = L"Magnetic $\mathbf{B}$ (nT)")
MAG_plot.plot_B_total(ax2, mag_data, time_range)
MAG_plot.plot_B(ax2, mag_data, time_range)
ax2.xticklabelsvisible=false
Legend(fig[2, 4], ax2)

ax3 = Axis(fig[3, 1:3], ylabel = L"Period T(s)", yscale = log2)
hm = MAG_plot.plot_wavelet(ax3, mag_data, time_range, 1/32.0)
Colorbar(fig[3, 4], hm, label = L"Wavelet Power $P_{\mathrm{B}}$")

# ---------- 对齐 ax2 与 ax3 时间轴 ----------
full_jd = datetime2julian.(time_range)
jd_lo, jd_hi = extrema(full_jd)
# 计算 plot_B 内部用的 t0 (第一个有效数据的 JulUTtime)
idx_r = findall(minimum(time_range) .<= mag_data[:epoch] .<= maximum(time_range))
if haskey(mag_data, :position)
    nan_r = any(isnan.(mag_data[:position][idx_r, :]), dims=2)[:]
    idx_r = idx_r[.!nan_r]
end
if haskey(mag_data, :B)
    nan_r = any(isnan.(mag_data[:B][idx_r, :]), dims=2)[:]
    idx_r = idx_r[.!nan_r]
end
t0_B = mag_data[:JulUTtime][idx_r[1]]
# 统一设置为相同物理时间范围
xlims!(ax2, jd_lo - t0_B, jd_hi - t0_B)
xlims!(ax3, jd_lo, jd_hi)

# ---------- 底部位置信息 ----------
tick_idx = [argmin(abs.(mag_data[:epoch] .- t)) for t in time_range]
xpos = mag_data[:position][tick_idx, 1]
ypos = mag_data[:position][tick_idx, 2]
zpos = mag_data[:position][tick_idx, 3]

ax3.xticks = (full_jd, [
    string(
        Dates.format(time_range[i], "HH:MM"),
        "\n",
        round(xpos[i]; digits=3),
        "\n",
        round(ypos[i]; digits=3),
        "\n",
        round(zpos[i]; digits=3)
    )
    for i in eachindex(time_range)
])
# 表头：放在 ax3 左侧边栏（ylabel 同竖线），底部与 xtick 平齐
# Label(fig[3, 1, Left()], L"UT\n$X_{MSO}$\n$Y_{MSO}$\n$Z_{MSO}$",
#     fontsize = 22, halign = :left, valign = :bottom,
#     tellheight = false, tellwidth = false)

name = "MAG_Overview_" * Dates.format(time_range[1], "yyyymmddTHHMM") * 
       "_to_" * Dates.format(time_range[end], "HHMM") * ".pdf"
save(joinpath(Overview_path, name), fig, pt_per_unit = 1)
println("Saved figure: ", joinpath(Overview_path, name))