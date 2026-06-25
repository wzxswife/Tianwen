using Glob
using CommonDataFormat
using Dates
using Statistics
using CairoMakie
include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_plot.jl")
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
STATIC_path = joinpath(MAVEN_path, "STATIC", "l2", year_num, month_num)
MAVEN_output_path = joinpath(output_dir, "MAVEN")
STATIC_output_path = joinpath(MAVEN_output_path, "STATIC")
VDF_path = joinpath(STATIC_output_path, "VDF", year_num, month_num)
mkpath(VDF_path)
MAG_file = glob("mvn_mag_l2*$datestr*.sts", MAG_path)
STA_file = glob("mvn_sta_l2*32e4d16a8m*$datestr*.cdf", STATIC_path)
println("MAG file: ", MAG_file)
println("STATIC file: ", STA_file)
println("Output path: ", VDF_path)

mag_data = MAVEN_load.load_mag_l2(MAG_file[1])
sta_data = MAVEN_load.load_STATIC(STA_file[1])
println("Data loaded.")

