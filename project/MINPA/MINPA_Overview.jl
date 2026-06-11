using Dates, Glob
using GLMakie, CairoMakie
include("../../src/scripts/TW_MINPA.jl")
include("../../src/scripts/TW_load.jl")
include("../../src/scripts/MINPA_plot.jl")
using .TW_MINPA
using .TW_load
using .MINPA_plot

date = DateTime(2022, 08, 24, 00, 00, 00)
datejul = datetime2julian(date)
datestr = Dates.format.(date, "yyyymmdd")
ts = DateTime(2022, 08, 24, 8, 10, 00) 
te = DateTime(2022, 08, 24, 14, 00, 00)
time_range = ts .+ Dates.Minute(70) .* range(0, 5)
ts_str = Dates.format(ts, "yyyymmddTHH:MM")

# read TW data 
TW_path = "E:/Data/Tianwen-1/"
MAG_path = joinpath(TW_path, "MOMAG/2C.V03/2022/08/")
MINPA_path = joinpath(TW_path, "MINPA/2B/2022/08/")
MAG_file = joinpath(MAG_path, "TW1_MOMAG_MSO_01hz_20220824_2C_v03.dat")
MINPA_file_list = glob("*20220824*.2B", MINPA_path)

println("MINPA files found: ", length(MINPA_file_list))
for file in MINPA_file_list
    println("  ", file)
end

minpa_file = MINPA_file_list[2]
mod = 1
save_path = "E:/Output/Tianwen-1/MINPA/Overview"
mkpath(save_path)

mag_data = TW_load.load_mag_2c_bydlm(MAG_file)
minpa_data = TW_load.load_minpa(minpa_file, mod) 

im = 1
pic_name1 = "MINPA_H+eflux_Overview"* ts_str * ".pdf"
TW_MINPA_eflux_Overview = TW_MINPA_eflux_Overview_plot(minpa_data, time_range, im)
save(joinpath(save_path, pic_name1), TW_MINPA_eflux_Overview, pt_per_unit = 1)
println("Saved: ", pic_name1, " at ", save_path)
im = 2
pic_name2 = "MINPA_He++eflux_Overview"* ts_str * ".pdf"
TW_MINPA_eflux_Overview = TW_MINPA_eflux_Overview_plot(minpa_data, time_range, im)
save(joinpath(save_path, pic_name2), TW_MINPA_eflux_Overview, pt_per_unit = 1)
println("Saved: ", pic_name2, " at ", save_path)
im = 4
pic_name3 = "MINPA_O+eflux_Overview"* ts_str * ".pdf"
TW_MINPA_eflux_Overview = TW_MINPA_eflux_Overview_plot(minpa_data, time_range, im)
save(joinpath(save_path, pic_name3), TW_MINPA_eflux_Overview, pt_per_unit = 1)
println("Saved: ", pic_name3, " at ", save_path)
im = 6
pic_name6 = "MINPA_O2+eflux_Overview"* ts_str * ".pdf"
TW_MINPA_eflux_Overview = TW_MINPA_eflux_Overview_plot(minpa_data, time_range, im)
save(joinpath(save_path, pic_name6), TW_MINPA_eflux_Overview, pt_per_unit = 1)
println("Saved: ", pic_name6, " at ", save_path)
im = 7
pic_name7 = "MINPA_CO2+eflux_Overview"* ts_str * ".pdf"
TW_MINPA_eflux_Overview = TW_MINPA_eflux_Overview_plot(minpa_data, time_range, im)
save(joinpath(save_path, pic_name7), TW_MINPA_eflux_Overview, pt_per_unit = 1)
println("Saved: ", pic_name7, " at ", save_path)