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
Phase_path = joinpath(MOMAG_output_path, "Phase", year_num, month_num)
mkpath(Phase_path)
MOMAG_file = glob("TW1_MOMAG_MSO*$datestr*.dat", MOMAG_path)
println("MAG file: ", MOMAG_file)
println("Output path: ", Phase_path)

# data01Hz_raw = load_mag_2c_bydlm(MOMAG_file[1])
data32Hz_raw = TW_load.load_mag_2c_bydlm(MOMAG_file[2])
println("Data Loaded.")
# data01Hz[:position] = data01Hz[:position]./Rm
data32Hz_raw[:position] = data32Hz_raw[:position]./Rm

time_range = DateTime(date, Time(7, 40, 00)).+Dates.Minute(1) .* range(0, 2)
wave1 = DateTime(date, Time(7, 40, 20)).+Dates.Second(4) .* range(0, 4)

key_vals = [:epoch, :JulUTtime, :B_total, :B, :position]
# data01Hz = TW_load.find_avail_data(data01Hz_raw, time_range, key_vals)
data32Hz = TW_load.find_avail_data(data32Hz_raw, time_range, key_vals)
datawave1 = TW_load.find_avail_data(data32Hz, wave1, key_vals)

Bwave = datawave1[:B].-mean(datawave1[:B], dims=1)

CairoMakie.activate!()
size_inches = (30, 10)
size_pt = 72 .* size_inches
fig = Figure(size = size_pt, fontsize = 24, figure_padding = (30, 20, 20, 20), rowgap = 20)
set_theme!(;
    Axis=(backgroundcolor = :white,
        bordercolor=:black, bordersize=0, spinewidth = 2, 
        xlabelsize=26, ylabelsize=26, xticklabelsize=22, yticklabelsize=22, 
        xticksize=8, xminorticksize=4, yticksize=8, yminorticksize=4,
        xtickwidth = 2, ytickwidth = 2, xminortickwidth = 1, yminortickwidth = 1,
        xminorticksvisible = true, yminorticksvisible = true,
        xminorticks = IntervalsBetween(5), yminorticks = IntervalsBetween(5))
)

ax_phase1 =  Axis(fig[1,1])
ax_phase2 =  Axis(fig[1,2])
ax_phase3 =  Axis(fig[1,3])
for ax in [ax_phase1, ax_phase2, ax_phase3]
    ax.xtickwidth = 2
    ax.ytickwidth = 2
    ax.xminortickwidth = 1
    ax.yminortickwidth = 1
    ax.spinewidth = 2
end
plot_MVA(ax_phase1,ax_phase2,ax_phase3,Bwave)

name = "MOMAG_Phase_" * Dates.format(wave1[1], "yyyymmddTHHMM") * 
       "_to_" * Dates.format(wave1[end], "HHMM") * ".pdf"
# name = "MOMAG_Phase_" * Dates.format(time_range[1], "yyyymmddTHHMM") * 
#        "_to_" * Dates.format(time_range[end], "HHMM") * ".pdf"
save(joinpath(Phase_path, name), fig, pt_per_unit = 1)
println("Saved figure: ", joinpath(Phase_path, name))