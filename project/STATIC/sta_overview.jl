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
SWIA_path = joinpath(MAVEN_path, "SWIA", "l2", year_num, month_num)
STATIC_path = joinpath(MAVEN_path, "STATIC", "l2", year_num, month_num)
MAVEN_output_path = joinpath(output_dir, "MAVEN")
STATIC_output_path = joinpath(MAVEN_output_path, "STATIC")
Overview_path = joinpath(STATIC_output_path, "Overview", year_num, month_num)
mkpath(Overview_path)
MAG_file = glob("mvn_mag_l2*$datestr*.sts", MAG_path)
SWIA_file = glob("mvn_swi_l2_coarsesvy3d*$datestr*.cdf", SWIA_path)
STA_file = glob("mvn_sta_l2*32e4d16a8m*$datestr*.cdf", STATIC_path)
println("MAG file: ", MAG_file)
println("SWIA file: ", SWIA_file)
println("STATIC file: ", STA_file)
println("Output path: ", Overview_path)

sta_data = MAVEN_load.load_STATIC(STA_file[1])

img = 1
sp_title = ["H+", "He++", "m/q=4.573966", "m/q=9.22163", 
    "O+", "O2+", "CO2+", "m/q=69.26344"]
sp = [rich(rich("H"), superscript("+")), rich(rich("He"), superscript("++")), 
    rich("m/q=4.573966"), rich("m/q=9.22163"), rich(rich("O"), superscript("+")), 
    rich(rich("O"), subscript("2"), superscript("+")), rich(rich("CO"), 
    subscript("2"), superscript("+")), rich("m/q=69.26344")]

h_list = 0:3:24
img_list = [1, 2, 5, 6, 7]

for img in img_list
    title = rich("MAVEN STATIC ") * sp[img] * rich(rich(" Differential Energy Flux (cm"), 
        superscript("-2"), rich("s"), superscript("-1"), rich("sr"), 
        superscript("-1"), rich(")"))
    for hstart in h_list
        xd = DateTime.(date + Dates.Hour(hstart) .+ Dates.Minute(30) .* range(0, 6))
        name = "MVN_STATIC_" * sp_title[img] * "eflux_Overview_" * Dates.format(xd[1], "yyyymmddTHH") * ".pdf"
        fig = MAVEN_plot.MVN_STATIC_eflux_Overview_plot(img, sta_data, xd)
        Label(fig[0, :], title, fontsize=20, tellwidth=false)
        save(joinpath(Overview_path, name), fig, pt_per_unit = 1)
        println("Saved figure: ", joinpath(Overview_path, name))
    end
    xd_total = DateTime.(date + Dates.Hour(0) .+ Dates.Hour(4) .* range(0, 6)) 
    name = "MVN_STATIC_" * sp_title[img] * "eflux_Overview_" * Dates.format(xd_total[1], "yyyymmdd") * "_total.pdf"
    fig = MAVEN_plot.MVN_STATIC_eflux_Overview_plot(img, sta_data, xd_total)
    Label(fig[0, :], title, fontsize=20, tellwidth=false)
    save(joinpath(Overview_path, name), fig, pt_per_unit = 1)
    println("Saved figure: ", joinpath(Overview_path, name))
end