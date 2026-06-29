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
shock_time = DateTime(date, Time(7, 27, 20))
wave_time = DateTime(date, Time(7, 30, 18))
shockJUL = datetime2julian.(shock_time)
waveJUL = datetime2julian.(wave_time)

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
MAG_plot.plot_spacecraft_orbit_xr(ax1,mag_data,[time_range[1],shock_time,time_range[end]];
    xlim=xlim1,ylim=(0,5))

idx = findall(time_range[1].<mag_data[:epoch].<time_range[end])
ax2 = Axis(fig[2, 1:3], ylabel = L"Magnetic $\mathbf{B}$ (nT)")
MAG_plot.plot_B_total(ax2, mag_data, time_range)
MAG_plot.plot_B(ax2, mag_data, time_range)
ax2.xticklabelsvisible=false
ax2_ymax = max(maximum(mag_data[:B][idx,:]),maximum(mag_data[:B_total][idx]))+1.0
ylims!(min(minimum(mag_data[:B][idx,:]),minimum(mag_data[:B_total][idx]))-1.0, ax2_ymax)
Legend(fig[2, 4], ax2)

ax3 = Axis(fig[3, 1:3], ylabel = L"Period T(s)", yscale = log2)
hm = MAG_plot.plot_wavelet(ax3, mag_data, time_range, 1/32.0; ylim=(1/8.0, 8.0))
Colorbar(fig[3, 4], hm, label = L"Wavelet Power $P_{\mathrm{B}}$")

full_jd = datetime2julian.(time_range)

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
Label(fig[3, 1:3, BottomLeft()], rich("UT\n") * rich(rich("X"), subscript("MSO"), "\n") *
    rich(rich("Y"), subscript("MSO"), "\n") * rich(rich("Z"), subscript("MSO")),
    fontsize = 22, halign = :left,  valign = :bottom, 
    # offset = (10, -5),
    justification = :left, lineheight = 1.5,
    padding = (-10, 0, 0, 0),
    tellheight = true, tellwidth = true
)

if !isempty(shockJUL)
    vlines!(ax2, shockJUL; color=:black, linewidth=2, linestyle=:dash)
    vlines!(ax3, shockJUL; color=:black, linewidth=2, linestyle=:dash)
    text!(ax2, shockJUL, ax2_ymax; text=" Bow Shock",
        color=:black, fontsize=20, align=(:left, :top))
    text!(ax3, shockJUL, 8.0; text=" Bow Shock",
        color=:black, fontsize=20, align=(:left, :top))
end
if !isempty(waveJUL)
    vlines!(ax2, waveJUL; color=:red, linewidth=2, linestyle=:dash)
    vlines!(ax3, waveJUL; color=:red, linewidth=2, linestyle=:dash)
    text!(ax2, waveJUL, ax2_ymax; text=" Wave",
        color=:red, fontsize=20, align=(:left, :top))
    text!(ax3, waveJUL, 8.0; text=" Wave",
        color=:red, fontsize=20, align=(:left, :top))
end

# 子图序号标注
text!(ax1, 0.02, 0.98, text="(a)", space=:relative, fontsize=30, align=(:left, :top))
text!(ax2, 0.02, 0.98, text="(b)", space=:relative, fontsize=30, align=(:left, :top))
text!(ax3, 0.02, 0.98, text="(c)", space=:relative, fontsize=30, align=(:left, :top))

name = "MAG_Overview_" * Dates.format(time_range[1], "yyyymmddTHHMM") * 
       "_to_" * Dates.format(time_range[end], "HHMM") * ".pdf"
save(joinpath(Overview_path, name), fig, pt_per_unit = 1)
println("Saved figure: ", joinpath(Overview_path, name))
