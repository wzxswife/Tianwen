using Glob
using CommonDataFormat
using Dates, CSV, DataFrames
using Statistics
using CairoMakie
include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_plot.jl")
include("../../src/scripts/MAVEN_SWIA.jl")
using .MAVEN_load
using .MAVEN_plot
using .MAVEN_SWIA

date = Date(2022, 8, 24)
datestr = Dates.format.(date, "yyyymmdd")
year_num = Dates.format(date, "yyyy")
month_num = Dates.format(date, "mm")

# Define paths
data_dir = "E:\\Data"       # Windows path
output_dir = "E:\\Output"   # Windows path
MAVEN_path = joinpath(data_dir, "MAVEN")
MAG_path = joinpath(MAVEN_path, "MAG", "l2", year_num, month_num)
SWIA_path = joinpath(MAVEN_path, "SWIA", "l2", year_num, month_num)
quat_path = joinpath(MAVEN_path, "QUAT", year_num, month_num)
MAVEN_output_path = joinpath(output_dir, "MAVEN")
SWIA_output_path = joinpath(MAVEN_output_path, "SWIA")
Overview_path = joinpath(SWIA_output_path, "Overview", year_num, month_num)
mkpath(Overview_path)
MAG_file = glob("mvn_mag_l2*$datestr*.sts", MAG_path)
SWIA_file = glob("mvn_swi_l2_coarsesvy3d*$datestr*.cdf", SWIA_path)
quat_file = glob("mvn_spice_swia_qu_$(datestr).csv", quat_path)
println("MAG file: ", MAG_file)
println("SWIA file: ", SWIA_file)
println("Quaternion file: ", quat_file)
println("Output path: ", Overview_path)

# Load data
mag_data = MAVEN_load.load_mag_l2(MAG_file[1])
swi_data = MAVEN_load.load_cdf(SWIA_file[1])
quat_data = MAVEN_load.load_quat(quat_file[1])
quat_data[:data_load_flag] = true
swi_data = MAVEN_SWIA.get_3dc!(swi_data, quat_data = quat_data)
println("Data loaded.")

# Define time range for plotting
h_list = 0:3:24
title = rich("MAVEN SWIA ") * rich(rich("Ion Differential Energy Flux (cm"), 
    superscript("-2"), rich("s"), superscript("-1"), rich("sr"), 
    superscript("-1"), rich(")"))
for hstart in h_list
    xd = DateTime.(date + Dates.Hour(hstart) .+ Dates.Minute(30) .* range(0, 6))
    name = "MVN_SWIA_" * "eflux_Overview_" * Dates.format(xd[1], "yyyymmddTHH") * ".pdf"
    fig = MAVEN_plot.MVN_SWIA_eflux_Overview_plot(swi_data, xd)
    Label(fig[0, :], title, fontsize=20, tellwidth=false)
    save(joinpath(Overview_path, name), fig, pt_per_unit = 1)
    println("Saved figure: ", joinpath(Overview_path, name))
end
xd_total = DateTime.(date + Dates.Hour(0) .+ Dates.Hour(4) .* range(0, 6)) 
name = "MVN_SWIA_" * "eflux_Overview_" * Dates.format(xd_total[1], "yyyymmdd") * "_total.pdf"
fig = MAVEN_plot.MVN_SWIA_eflux_Overview_plot(swi_data, xd_total)
Label(fig[0, :], title, fontsize=20, tellwidth=false)
save(joinpath(Overview_path, name), fig, pt_per_unit = 1)
println("Saved figure: ", joinpath(Overview_path, name))