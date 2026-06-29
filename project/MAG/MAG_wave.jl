using Glob
using CommonDataFormat
using Dates
using Statistics
using CairoMakie
include("../../src/scripts/TW_load.jl")
include("../../src/scripts/MAG_plot.jl")
include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_plot.jl")
using .TW_load
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
Wave_path = joinpath(MAG_output_path, "Wave", year_num, month_num)
mkpath(Wave_path)
MAG_file = glob("mvn_mag_l2*$datestr*.sts", MAG_path)
KP_file = glob("mvn_kp*$datestr*.jld2", KP_path)
SWIA_file = glob("mvn_swi_l2_coarsesvy3d*$datestr*.cdf", SWIA_path)
STA_file = glob("mvn_sta_l2*32e4d16a8m*$datestr*.cdf", STATIC_path)
println("MAG file: ", MAG_file)
println("KP file: ", KP_file)
println("SWIA file: ", SWIA_file)
println("STATIC file: ", STA_file)
println("Output path: ", Wave_path)

mag_data = load_mag_l2(MAG_file[1])
println("Data Loaded.")
const Rm = 3389.5
mag_data[:position] = mag_data[:position]./Rm

time_range = DateTime(date, Time(7, 30, 10)).+Dates.Second(5) .* range(0, 4)
wave = DateTime(date, Time(7, 30, 18)).+Dates.Second(1) .* range(0, 4)

key_vals = [:epoch, :JulUTtime, :B_total, :B, :position]
wave_data = TW_load.find_avail_data(mag_data, time_range, key_vals)

Bwave = Float64.(wave_data[:B].-mean(wave_data[:B], dims=1))

CairoMakie.activate!()
size_inches = (30, 30)
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
ax_Phase1 =  Axis(fig[3,1], aspect=1)
ax_Phase2 =  Axis(fig[3,2], aspect=1)
ax_Phase3 =  Axis(fig[3,3], aspect=1)
MAG_plot.plot_wave(ax1, wave_data, time_range;lw=2)
MAG_plot.plot_wave(ax2, wave_data, wave;lw=2)
MAG_plot.plot_MVA(ax_Phase1,ax_Phase2,ax_Phase3,Bwave;lw=4)
wavestr = Dates.format(time_range[1], "yyyy-mm-ddTHH:MM:SS")*"-"*Dates.format(time_range[end], "HH:MM:SS")
wavestr1 = Dates.format(wave[1], "yyyy-mm-ddTHH:MM:SS")*"-"*Dates.format(wave[end], "HH:MM:SS")
timestr = Dates.format(time_range[1], "HH:MM:SS")*"-"*Dates.format(time_range[end], "HH:MM:SS")
text!(ax1, 0.02, 0.98, text="(a) $wavestr", space=:relative, fontsize=30, align=(:left, :top))
text!(ax2, 0.02, 0.98, text="(b) $wavestr1", space=:relative, fontsize=30, align=(:left, :top))
text!(ax_Phase1, 0.02, 0.98, text="(c) $timestr", space=:relative, fontsize=30, align=(:left, :top))
text!(ax_Phase2, 0.02, 0.98, text="(d) $timestr", space=:relative, fontsize=30, align=(:left, :top))
text!(ax_Phase3, 0.02, 0.98, text="(e) $timestr", space=:relative, fontsize=30, align=(:left, :top))


name = "MAG_Wave_" * Dates.format(time_range[1], "yyyymmddTHHMM") * 
       "_to_" * Dates.format(time_range[end], "HHMM") * ".pdf"
save(joinpath(Wave_path, name), fig, pt_per_unit = 1)
println("Saved figure: ", joinpath(Wave_path, name))
