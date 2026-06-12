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
pitch_angle_path = joinpath(STATIC_output_path, "PitchAngle", year_num, month_num)
mkpath(pitch_angle_path)
MAG_file = glob("mvn_mag_l2*$datestr*.sts", MAG_path)
STA_file = glob("mvn_sta_l2*32e4d16a8m*$datestr*.cdf", STATIC_path)
println("MAG file: ", MAG_file)
println("STATIC file: ", STA_file)
println("Output path: ", pitch_angle_path)

mag_data = MAVEN_load.load_mag_l2(MAG_file[1])
sta_data = MAVEN_load.load_STATIC(STA_file[1])
println("Data loaded.")

img = 1
sp_title = ["H+", "He++", "m/q=4.573966", "m/q=9.22163", 
    "O+", "O2+", "CO2+", "m/q=69.26344"]
sp = [rich(rich("H"), superscript("+")), rich(rich("He"), superscript("++")), 
    rich("m/q=4.573966"), rich("m/q=9.22163"), rich(rich("O"), superscript("+")), 
    rich(rich("O"), subscript("2"), superscript("+")), rich(rich("CO"), 
    subscript("2"), superscript("+")), rich("m/q=69.26344")]

h_list = 0:3:24
img_list = [1, 2, 5, 6, 7]

using CairoMakie
CairoMakie.activate!()
set_theme!(;
    Axis=(; bordercolor=:black, bordersize=15, spinewidth = 3, 
        xlabelsize=18, ylabelsize=18, xticklabelsize=16, yticklabelsize=16)
)

for img in img_list
    title = rich("MAVEN STATIC ") * sp[img] * rich(rich(" Differential Energy Flux (cm"), 
        superscript("-2"), rich("s"), superscript("-1"), rich("sr"), 
        superscript("-1"), rich(")"))
    for hstart in h_list
        time_range = DateTime.(date + Dates.Hour(hstart) .+ Dates.Minute(30) .* range(0, 6))
        fig = Figure(size = (1600, 1200))
        yticks = 0:30:180
        # ax1 = Axis(fig[1, 1:3], ylabel=sp[img] * " 4-10 keV\nφ angle (°)", title=title,
        #     titlesize=20)
        ax2 = Axis(fig[1, 1:3], ylabel=sp[img] * " <4 keV\npitch angle (°)")
        ax3 = Axis(fig[2, 1:3], ylabel=sp[img] * " 4-10 keV\npitch angle (°)")
        ax4 = Axis(fig[3, 1:3], xlabel="Time (UT)", ylabel=sp[img] * " 10-25 keV\npitch angle (°)")
        # hm1 = MAVEN_plot.sta_phi_angle(ax1, sta_data, time_range, [4.0, 10.0]; 
        #     c_range=(1e3, 1e7))
        hm2 = MAVEN_plot.sta_pitch_angle(ax2, sta_data, mag_data, img, time_range, 
            [0.0, 4.0]; c_range=(1e5, 1e8))
        hm3 = MAVEN_plot.sta_pitch_angle(ax3, sta_data, mag_data, img, time_range, 
            [4.0, 10.0]; c_range=(1e5, 1e8))
        hm4 = MAVEN_plot.sta_pitch_angle(ax4, sta_data, mag_data, img, time_range, 
            [10.0, 25.0]; c_range=(1e5, 1e8))
        Colorbar(fig[1:3, 4], hm2, label="Differential Energy Flux (eV/cm²/s/sr/eV)", 
            width=15, labelsize=14, ticklabelsize=12)
        # ax1.xticklabelsvisible=false
        ax2.xticklabelsvisible=false
        ax3.xticklabelsvisible=false
        name = "MVN_STATIC_" * sp_title[img] * " pitch_angle_" * Dates.format(time_range[1], "yyyymmddTHH") * ".pdf"
        # Label(fig[0, :], title, fontsize=20, tellwidth=false)
        save(joinpath(pitch_angle_path, name), fig, pt_per_unit = 1)
        println("Saved figure: ", joinpath(pitch_angle_path, name))
    end
end
