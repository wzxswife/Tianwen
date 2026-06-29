using Glob
using CommonDataFormat
using Dates
using Statistics
using CairoMakie
include("../../src/scripts/TW_load.jl")
include("../../src/scripts/MAG_plot.jl")
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
Wave_path = joinpath(MOMAG_output_path, "Wave", year_num, month_num)
mkpath(Wave_path)
MOMAG_file = glob("TW1_MOMAG_MSO*$datestr*.dat", MOMAG_path)
println("MAG file: ", MOMAG_file)
println("Output path: ", Wave_path)

# data01Hz_raw = TW_load.load_mag_2c_bydlm(MOMAG_file[1])
data32Hz_raw = TW_load.load_mag_2c_bydlm(MOMAG_file[2])
println("Data Loaded.")
# data01Hz[:position] = data01Hz[:position]./Rm
data32Hz_raw[:position] = data32Hz_raw[:position]./Rm

time_range = DateTime(date, Time(7, 40, 00)).+Dates.Second(20) .* range(0, 6)
wave1 = DateTime(date, Time(7, 40, 10)).+Dates.Second(5) .* range(0, 4)
wave2 = DateTime(date, Time(7, 40, 30)).+Dates.Second(5) .* range(0, 4)
wave3 = DateTime(date, Time(7, 41, 20)).+Dates.Second(5) .* range(0, 4)

key_vals = [:epoch, :JulUTtime, :B_total, :B, :position]
# data01Hz = TW_load.find_avail_data(data01Hz_raw, time_range, key_vals)
data32Hz = TW_load.find_avail_data(data32Hz_raw, time_range, key_vals)
datawave = TW_load.find_avail_data(data32Hz, time_range, key_vals)
datawave1 = TW_load.find_avail_data(data32Hz, wave1, key_vals)
datawave2 = TW_load.find_avail_data(data32Hz, wave2, key_vals)
datawave3 = TW_load.find_avail_data(data32Hz, wave3, key_vals)

Bwave = datawave[:B].-mean(datawave[:B], dims=1)
Bwave1 = datawave1[:B].-mean(datawave1[:B], dims=1)
Bwave2 = datawave1[:B].-mean(datawave2[:B], dims=1)
Bwave3 = datawave1[:B].-mean(datawave3[:B], dims=1)

CairoMakie.activate!()
size_inches = (30, 40)
size_pt = 72 .* size_inches
fig = Figure(size = size_pt, fontsize = 24, figure_padding = (20, 80, 20, 20), rowgap = 20)
set_theme!(;
    Axis=(backgroundcolor = :white,
        bordercolor=:black, bordersize=0, spinewidth = 2, 
        xlabelsize=26, ylabelsize=26, xticklabelsize=22, yticklabelsize=22, 
        xticksize=8, xminorticksize=4, yticksize=8, yminorticksize=4,
        xtickwidth = 2, ytickwidth = 2, xminortickwidth = 1, yminortickwidth = 1,
        xminorticksvisible = true, yminorticksvisible = true,
        xminorticks = IntervalsBetween(5), yminorticks = IntervalsBetween(5))
)

ax1 =  Axis(fig[1,1:3], xlabel="Time (UT)", ylabel = "Magnetic B (nT)")
ax2 =  Axis(fig[2,1:3], xlabel="Time (UT)", ylabel = "Magnetic B (nT)")
ax3 =  Axis(fig[3,1:3], xlabel="Time (UT)", ylabel = "Magnetic B (nT)")
ax_Phase1 =  Axis(fig[4,1], aspect=1)
ax_Phase2 =  Axis(fig[4,2], aspect=1)
ax_Phase3 =  Axis(fig[4,3], aspect=1)
# ax_Phase1 =  Axis(fig[2,1], aspect=1)
# ax_Phase2 =  Axis(fig[2,2], aspect=1)
# ax_Phase3 =  Axis(fig[2,3], aspect=1)
# MAG_plot.plot_wave(ax1, data32Hz, time_range;lw=2)
MAG_plot.plot_wave(ax1, data32Hz, wave1;lw=2)
MAG_plot.plot_wave(ax2, data32Hz, wave2;lw=2)
MAG_plot.plot_wave(ax3, data32Hz, wave3;lw=2)
MAG_plot.plot_MVA(ax_Phase1,ax_Phase2,ax_Phase3,Bwave2; lw=4)
# wavestr = Dates.format(time_range[1], "yyyy-mm-ddTHH:MM:SS")*"-"*Dates.format(time_range[end], "HH:MM:SS")
wavestr1 = Dates.format(wave1[1], "yyyy-mm-ddTHH:MM:SS")*"-"*Dates.format(wave1[end], "HH:MM:SS")
wavestr2 = Dates.format(wave2[1], "yyyy-mm-ddTHH:MM:SS")*"-"*Dates.format(wave2[end], "HH:MM:SS")
wavestr3 = Dates.format(wave3[1], "yyyy-mm-ddTHH:MM:SS")*"-"*Dates.format(wave3[end], "HH:MM:SS")
timestr = Dates.format(wave2[1], "HH:MM:SS")*"-"*Dates.format(wave2[end], "HH:MM:SS")
# text!(ax1, 0.02, 0.98, text="(a) $wavestr", space=:relative, fontsize=30, align=(:left, :top))
text!(ax1, 0.02, 0.98, text="(a) $wavestr1", space=:relative, fontsize=30, align=(:left, :top))
text!(ax2, 0.02, 0.98, text="(b) $wavestr2", space=:relative, fontsize=30, align=(:left, :top))
text!(ax3, 0.02, 0.98, text="(c) $wavestr3", space=:relative, fontsize=30, align=(:left, :top))
text!(ax_Phase1, 0.02, 0.98, text="(d) $timestr", space=:relative, fontsize=30, align=(:left, :top))
text!(ax_Phase2, 0.02, 0.98, text="(e) $timestr", space=:relative, fontsize=30, align=(:left, :top))
text!(ax_Phase3, 0.02, 0.98, text="(f) $timestr", space=:relative, fontsize=30, align=(:left, :top))

name = "MOMAG_Wave_" * Dates.format(time_range[1], "yyyymmddTHHMM") * 
       "_to_" * Dates.format(time_range[end], "HHMM") * ".pdf"
save(joinpath(Wave_path, name), fig, pt_per_unit = 1)
println("Saved figure: ", joinpath(Wave_path, name))
