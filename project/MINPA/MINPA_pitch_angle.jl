using Glob
using CommonDataFormat
using Dates, CSV, DataFrames
using LinearAlgebra
using Statistics
using CairoMakie
include("../../src/scripts/TW_load.jl")
include("../../src/scripts/MINPA_plot.jl")
include("../../src/scripts/MAVEN_SWIA.jl")
using .TW_load
using .MINPA_plot
using .MAVEN_SWIA

date = Date(2022, 8, 24)
datestr = Dates.format.(date, "yyyymmdd")
year_num = Dates.format(date, "yyyy")
month_num = Dates.format(date, "mm")

data_dir = "E:\\Data"
output_dir = "E:\\Output"
TW_path = joinpath(data_dir, "Tianwen-1")
MINPA_path = joinpath(TW_path, "MINPA", "2B", year_num, month_num)
MOMAG_path = joinpath(TW_path, "MOMAG", "2C.V03", year_num, month_num)
TW_output_path = joinpath(output_dir, "Tianwen-1")
MINPA_output_path = joinpath(TW_output_path, "MINPA")
pitch_angle_path = joinpath(MINPA_output_path, "pitch_angle", year_num, month_num)
mkpath(pitch_angle_path)
MINPA_file_list = glob("*$datestr*.2B", MINPA_path)
MOMAG_file = glob("*01Hz*$datestr*.dat", MOMAG_path)

println("MOMAG file: ", MOMAG_file)
println("MINPA files found: ", length(MINPA_file_list))
for file in MINPA_file_list
    println("  ", file)
end
println("Output path: ", pitch_angle_path)

mag_data = TW_load.load_mag_2c_bydlm(MOMAG_file[1])
mod = 4
minpa_data = TW_load.load_minpa(MINPA_file_list[6],mod)

h_list = [7]
img_list = [1, 2, 5, 6, 7]

sp_title = ["H+", "He++", "m/q=4.573966", "m/q=9.22163", 
    "O+", "O2+", "CO2+", "m/q=69.26344"]
sp = [rich(rich("H"), superscript("+")), rich(rich("He"), superscript("++")), 
    rich("m/q=4.573966"), rich("m/q=9.22163"), rich(rich("O"), superscript("+")), 
    rich(rich("O"), subscript("2"), superscript("+")), rich(rich("CO"), 
    subscript("2"), superscript("+")), rich("m/q=69.26344")]

using CairoMakie
CairoMakie.activate!()
set_theme!(;
    Axis=(; bordercolor=:black, bordersize=15, spinewidth = 3, 
        xlabelsize=18, ylabelsize=18, xticklabelsize=16, yticklabelsize=16)
)
for img in img_list
    title = rich("MINPA ") * sp[img] * rich(rich(" Differential Energy Flux (cm"), 
        superscript("-2"), rich("s"), superscript("-1"), rich("sr"), 
        superscript("-1"), rich(")"))
    for h in h_list
        time_range = DateTime.(date+Dates.Hour(h).+Dates.Minute(10).*range(0, 6))
        fig = Figure(size = (1600, 1200))
        yticks = 0:30:180
        # ax1 = Axis(fig[1, 1:3], ylabel=sp[img] * " 4-10 keV\nφ angle (°)", title=title,
        #     titlesize=20)
        ax2 = Axis(fig[1, 1:3], ylabel=sp[img] * " <4 keV\npitch angle (°)")
        ax3 = Axis(fig[2, 1:3], ylabel=sp[img] * " 4-10 keV\npitch angle (°)")
        ax4 = Axis(fig[3, 1:3], xlabel="Time (UT)", ylabel=sp[img] * " 10-25 keV\npitch angle (°)")
        # hm1 = MINPA_plot.TW_MINPA_phi_angle(ax1, minpa_data, time_range, [4.0, 10.0]; 
        #     c_range=(1e3, 1e7))
        hm2 = MINPA_plot.TW_MINPA_pitch_angle(ax2, minpa_data, mag_data, img, time_range, 
            [0.0, 4.0]; c_range=(1e5, 1e8))
        hm3 = MINPA_plot.TW_MINPA_pitch_angle(ax3, minpa_data, mag_data, img, time_range, 
            [4.0, 10.0]; c_range=(1e5, 1e8))
        hm4 = MINPA_plot.TW_MINPA_pitch_angle(ax4, minpa_data, mag_data, img, time_range, 
            [10.0, 25.0]; c_range=(1e5, 1e8))
        Colorbar(fig[1:3, 4], hm2, label="Differential Energy Flux (eV/cm²/s/sr/eV)", 
            width=15, labelsize=14, ticklabelsize=12)
        # ax1.xticklabelsvisible=false
        ax2.xticklabelsvisible=false
        ax3.xticklabelsvisible=false
        name = "MINPA_" * sp_title[img] * " pitch_angle_" * Dates.format(time_range[1], "yyyymmddTHH") * ".pdf"
        # Label(fig[0, :], title, fontsize=20, tellwidth=false)
        save(joinpath(pitch_angle_path, name), fig, pt_per_unit = 1)
        println("Saved figure: ", joinpath(pitch_angle_path, name))
    end
    time_range = DateTime.(date+Dates.Hour(h_list[1]).+Dates.Minute(10).*range(0, 6))
    fig = Figure(size = (1600, 1200))
    yticks = 0:30:180
    # ax1 = Axis(fig[1, 1:3], ylabel=sp[img] * " 4-10 keV\nφ angle (°)", title=title,
    #     titlesize=20)
    ax2 = Axis(fig[1, 1:3], ylabel=sp[img] * " <4 keV\npitch angle (°)")
    ax3 = Axis(fig[2, 1:3], ylabel=sp[img] * " 4-10 keV\npitch angle (°)")
    ax4 = Axis(fig[3, 1:3], xlabel="Time (UT)", ylabel=sp[img] * " 10-25 keV\npitch angle (°)")
    # hm1 = MINPA_plot.TW_MINPA_phi_angle(ax1, minpa_data, time_range, [4.0, 10.0]; 
    #     c_range=(1e3, 1e7))
    hm2 = MINPA_plot.TW_MINPA_pitch_angle(ax2, minpa_data, mag_data, img, time_range, 
        [0.0, 4.0]; c_range=(1e5, 1e8))
    hm3 = MINPA_plot.TW_MINPA_pitch_angle(ax3, minpa_data, mag_data, img, time_range, 
        [4.0, 10.0]; c_range=(1e5, 1e8))
    hm4 = MINPA_plot.TW_MINPA_pitch_angle(ax4, minpa_data, mag_data, img, time_range, 
        [10.0, 25.0]; c_range=(1e5, 1e8))
    Colorbar(fig[1:3, 4], hm2, label="Differential Energy Flux (eV/cm²/s/sr/eV)", 
        width=15, labelsize=14, ticklabelsize=12)
    # ax1.xticklabelsvisible=false
    ax2.xticklabelsvisible=false
    ax3.xticklabelsvisible=false
    name = "MINPA_" * sp_title[img] * " pitch_angle_" * Dates.format(time_range[1], "yyyymmddTHH") * ".pdf"
    # Label(fig[0, :], title, fontsize=20, tellwidth=false)
    save(joinpath(pitch_angle_path, name), fig, pt_per_unit = 1)
    println("Saved figure: ", joinpath(pitch_angle_path, name))
end