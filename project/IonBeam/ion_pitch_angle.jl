using CommonDataFormat
using Dates, CSV, DataFrames
using Statistics
using CairoMakie
include("../../src/scripts/TW_load.jl")
include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_plot.jl")
include("../../src/scripts/MAVEN_SWIA.jl")
include("../../src/scripts/MAVEN_STATIC.jl")
using .MAVEN_load
using .MAVEN_load
using .MAVEN_plot
using .MAVEN_SWIA
using .MAVEN_STATIC

date = Date(2022, 8, 24)
datestr = Dates.format.(date, "yyyymmdd")
time_stemp = DateTime(date, Time(7, 00, 00)) .+ Dates.Minute(20) .* range(0, 6)

date_str = Dates.format.(time_stemp[1], "yyyymmdd")
date_str_dye = "$(year(time_stemp[1]))$(lpad(dayofyear(time_stemp[1]), 3, '0'))"

save_name = "MAVEN_SWIA_Ion_VDF_2D_$(replace(string(time_stemp[1]), ":" => "-")).png"
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

energy_range = [4.0, 10.0]

using CairoMakie
CairoMakie.activate!()
fig = Figure(size = (1600, 400))
ax = Axis(fig[1, 1:3], xlabel="Time (UT)", ylabel=L"pitch angle (°)")
hm = swi_pitch_angle(ax, swi_data, time_stemp, energy_range; 
    c_range=(1e4, 1e7))
Colorbar(fig[1, 4], hm, label="Differential Energy Flux (eV/cm²/s/sr/eV)", width=15)
pic_name = "MAVEN.png"
save(joinpath(output_path, pic_name), fig)