using CommonDataFormat
using Dates, CSV, DataFrames
using Statistics
using CairoMakie
include("../../src/scripts/TW_load.jl")
include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_plot.jl")
include("../../src/scripts/MAVEN_SWIA.jl")
using .MAVEN_load
using .MAVEN_load
using .MAVEN_plot
using .MAVEN_SWIA

date = Date(2022, 8, 24)
datestr = Dates.format.(date, "yyyymmdd")
time_range = DateTime(date, Time(6, 30, 00)) .+ Dates.Minute(20) .* range(0, 6)

date_str = Dates.format.(time_range[1], "yyyymmdd")
date_str_dye = "$(year(time_range[1]))$(lpad(dayofyear(time_range[1]), 3, '0'))"

save_name = "MAVEN_SWIA_Ion_VDF_2D_$(replace(string(time_range[1]), ":" => "-")).png"
# dir = pwd()
dir = joinpath(@__DIR__, "..", "..")
input_path = joinpath(dir, "Data", "MAVEN")
output_path = joinpath(dir, "Results", "MAVEN", "PitchAngle")
save_path = joinpath(dir, "results", "IonBeam", save_name)
data_path = joinpath(input_path, "mvn_swi_l2_coarsesvy3d_$(date_str)_v02_r01.cdf")
quat_path = joinpath(input_path, "mvn_spice_swia_qu_$(date_str).csv")
mag_path = joinpath(input_path, "mvn_mag_l3_$(date_str_dye)ss1s_$(date_str)_v01_r01.f77_unformatted")

data = load_cdf(data_path)
quat_data = load_quat(quat_path)
quat_data[:data_load_flag] = true
swi_data = get_3dc!(data, quat_data = quat_data)
mag_data = load_mag_l3(mag_path)
println("Data loaded.")

energy_range = [4.0, 10.0]
title = "MAVEN Ion Pitch Angle Distribution\n
    $(Dates.format(time_range[1], "yyyy-mm-dd HH:MM:SS")) - $(Dates.format(time_range[end], "yyyy-mm-dd HH:MM:SS"))UT"
pic_name = "MAVEN_Ion_Pitch_Angle_$(Dates.format(time_range[1], "yyyy-mm-dd_HH-MM")).png"

using CairoMakie
CairoMakie.activate!()
set_theme!(;
    Axis=(; bordercolor=:black, bordersize=15, spinewidth = 3, 
        xlabelsize=18, ylabelsize=18, xticklabelsize=16, yticklabelsize=16)
)
fig = Figure(size = (1600, 1200))
yticks = 0:30:180
ax1 = Axis(fig[1, 1:3], ylabel="ions 4-10 keV\nφ angle (°)", title=title,
    titlesize=20)
ax2 = Axis(fig[2, 1:3], ylabel="ions <4 keV\npitch angle (°)")
ax3 = Axis(fig[3, 1:3], ylabel="ions 4-10 keV\npitch angle (°)")
ax4 = Axis(fig[4, 1:3], xlabel="Time (UT)", ylabel="ions 10-25 keV\npitch angle (°)")

hm1 = swi_phi_angle(ax1, swi_data, time_range, energy_range; 
    c_range=(1e3, 1e7))
hm2 = swi_pitch_angle(ax2, swi_data, mag_data, time_range, [0.0, 4.0]; 
     c_range=(1e5, 1e8))
hm3 = swi_pitch_angle(ax3, swi_data, mag_data, time_range, [4.0, 10.0]; 
    c_range=(1e5, 1e8))
hm4 = swi_pitch_angle(ax4, swi_data, mag_data, time_range, [10.0, 25.0]; 
    c_range=(1e5, 1e8))
Colorbar(fig[1, 4], hm1, label="Differential Energy Flux (eV/cm²/s/sr/eV)", 
    width=15, labelsize=14, ticklabelsize=12)
Colorbar(fig[2, 4], hm2, label="Differential Energy Flux (eV/cm²/s/sr/eV)", 
    width=15, labelsize=14, ticklabelsize=12)
Colorbar(fig[3, 4], hm3, label="Differential Energy Flux (eV/cm²/s/sr/eV)", 
    width=15, labelsize=14, ticklabelsize=12)
Colorbar(fig[4, 4], hm4, label="Differential Energy Flux (eV/cm²/s/sr/eV)", 
    width=15, labelsize=14, ticklabelsize=12)
ax1.xticklabelsvisible=false
ax2.xticklabelsvisible=false
ax3.xticklabelsvisible=false
save(joinpath(output_path, pic_name), fig)
println("Phi angle plot saved to $(joinpath(output_path, pic_name)).")