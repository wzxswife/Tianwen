using Glob
using CommonDataFormat
using Dates
using Statistics
using CairoMakie
include("../../src/scripts/MAG_plot.jl")
include("../../src/scripts/TW_load.jl")
using .MAG_plot
using .TW_load
const Rm = 3389.5

date = Date(2022, 8, 24)
datestr = Dates.format.(date, "yyyymmdd")
year_num = Dates.format(date, "yyyy")
month_num = Dates.format(date, "mm")

data_dir = "E:\\Data"
output_dir = "E:\\Output"
TW_path = joinpath(data_dir, "Tianwen-1")
MOMAG_path = joinpath(TW_path, "MOMAG", "2C.V03", year_num, month_num)
TW_output_path = joinpath(output_dir, "Tianwen-1")
MOMAG_output_path = joinpath(TW_output_path, "MOMAG")
Overview_path = joinpath(MOMAG_output_path, "Overview", year_num, month_num)
mkpath(Overview_path)
MOMAG_file = glob("TW1_MOMAG_MSO*$datestr*.dat", MOMAG_path)
println("MAG file: ", MOMAG_file)
println("Output path: ", Overview_path)

data01Hz_raw = TW_load.load_mag_2c_bydlm(MOMAG_file[1])
data32Hz_raw = TW_load.load_mag_2c_bydlm(MOMAG_file[2])
println("Data Loaded.")
# data01Hz[:position] = data01Hz[:position]./Rm
data32Hz_raw[:position] = data32Hz_raw[:position]./Rm

time_range = DateTime(date, Time(7, 30, 00)).+Dates.Minute(6) .* range(0, 6)
shock_time = DateTime(date, Time(7, 55, 20))
wave_time = DateTime(date, Time(7, 40, 00))
shockJUL = datetime2julian.(shock_time)
waveJUL = datetime2julian.(wave_time)

key_vals = [:epoch, :JulUTtime, :B_total, :B, :position]
data01Hz = TW_load.find_avail_data(data01Hz_raw, time_range, key_vals)
data32Hz = TW_load.find_avail_data(data32Hz_raw, time_range, key_vals)

CairoMakie.activate!()
size_inches = (20, 18)
size_pt = 72 .* size_inches
fig = Figure(size = size_pt, fontsize = 24, figure_padding = (30, 10, 20, 20), rowgap = 20)
set_theme!(;
    Axis=(backgroundcolor = :white,
        bordercolor=:black, bordersize=15, spinewidth = 2, 
        xlabelsize=26, ylabelsize=26, xticklabelsize=22, yticklabelsize=22, 
        xticksize=8, xminorticksize=4, yticksize=8, yminorticksize=4,
        xtickwidth = 2, ytickwidth = 2, xminortickwidth = 1, yminortickwidth = 1,
        xminorticksvisible = true, yminorticksvisible = true,
        xminorticks = IntervalsBetween(10), yminorticks = IntervalsBetween(10))
)

xlim1 = (6,-12)
ax1 = Axis(fig[1,1:3], aspect=DataAspect(), 
    xlabel=L"$x$ ($R_{\mathrm{M}}$)", 
    ylabel=L"$\sqrt{y^2+z^2}$ ($R_{\mathrm{M}}$)")
MAG_plot.bowshock_plot_xr(ax1;xlim=xlim1,ylim=(0,5))
MAG_plot.plot_spacecraft_orbit_xr(ax1,data01Hz_raw,[time_range[1],shock_time,time_range[end]];
    xlim=xlim1,ylim=(0,5))

idx = findall(time_range[1].<data01Hz[:epoch].<time_range[end])
ax2 = Axis(fig[2, 1:3], ylabel = L"Magnetic $\mathbf{B}$ (nT)")
MAG_plot.plot_B_total(ax2, data01Hz, time_range; lw=1.5)
MAG_plot.plot_B(ax2, data01Hz, time_range; lw=1.5)
ax2.xticklabelsvisible=false
ax2_ymax = max(maximum(data01Hz[:B][idx,:]),maximum(data01Hz[:B_total][idx]))+1.0
ylims!(min(minimum(data01Hz[:B][idx,:]),minimum(data01Hz[:B_total][idx]))-1.0, ax2_ymax)
Legend(fig[2, 4], ax2)

ax3 = Axis(fig[3, 1:3], ylabel = L"Period T(s)", yscale = log2)
hm = MAG_plot.plot_wavelet(ax3, data32Hz, time_range, 1/32.0; ylim=(1/8.0, 8.0))
Colorbar(fig[3, 4], hm, label = L"Wavelet Power $P_{\mathrm{B}}$")

full_jd = datetime2julian.(time_range)

# dx = Dates.Second(10)

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

# ---------- 底部位置信息 ----------
tick_idx = [argmin(abs.(data01Hz[:epoch] .- t)) for t in time_range]
xpos = data01Hz[:position][tick_idx, 1]
ypos = data01Hz[:position][tick_idx, 2]
zpos = data01Hz[:position][tick_idx, 3]

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
    justification = :left, lineheight = 1.5,
    padding = (-15, 0, 0, 0),
    tellheight = true, tellwidth = true
)

# 子图序号标注
text!(ax1, 0.02, 0.98, text="(a)", space=:relative, fontsize=30, align=(:left, :top))
text!(ax2, 0.02, 0.98, text="(b)", space=:relative, fontsize=30, align=(:left, :top))
# text!(ax2, 0.02, 0.98, text="(b) 1Hz Sampling Rate", space=:relative, fontsize=30, align=(:left, :top))
text!(ax3, 0.02, 0.98, text="(c)", space=:relative, fontsize=30, align=(:left, :top))

name = "MOMAG_Overview_" * Dates.format(time_range[1], "yyyymmddTHHMM") * 
       "_to_" * Dates.format(time_range[end], "HHMM") * ".pdf"
save(joinpath(Overview_path, name), fig, pt_per_unit = 1)
println("Saved figure: ", joinpath(Overview_path, name))