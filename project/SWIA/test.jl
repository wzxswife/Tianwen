using CommonDataFormat
using Dates, CSV, DataFrames
using Statistics
using CairoMakie
include("../../src/scripts/TW_load.jl")
include("../../src/scripts/wave_caculate.jl")
include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_plot.jl")
include("../../src/scripts/MAVEN_SWIA.jl")
# include("../../src/scripts/MAVEN_STATIC.jl")
using .TW_load
using .WaveCaculate
using .MAVEN_load
using .MAVEN_load
using .MAVEN_plot
using .MAVEN_SWIA
# using .MAVEN_STATIC

date = Date(2022, 8, 24)
time_range = DateTime(date, Time(7, 00, 00)) .+ Dates.Minute(20) .* range(0, 6)
date_str = Dates.format.(date, "yyyymmdd")
date_str_dye = "$(year(date))$(lpad(dayofyear(date), 3, '0'))"
# dir = pwd()
dir = joinpath(@__DIR__, "..", "..")
input_path = joinpath(dir, "Data", "MAVEN")
# save_path = joinpath(dir, "results", "MAVEN", save_name)
data_path = joinpath(input_path, "mvn_swi_l2_coarsesvy3d_$(date_str)_v02_r01.cdf")
quat_path = joinpath(input_path, "mvn_spice_swia_qu_$(date_str).csv")
# mag_path = joinpath(input_path, "mvn_mag_l3_$(date_str_dye)ss1s_$(date_str)_v01_r01.f77_unformatted")

data = MAVEN_load.load_cdf(data_path)
quat_data = MAVEN_load.load_quat(quat_path)
quat_data[:data_load_flag] = true
swi_data = MAVEN_SWIA.get_3dc!(data, quat_data = quat_data)
# mag_data = MAVEN_load.load_mag_l3(mag_path)
println("Data loaded.")

time_slice = DateTime(date, Time(7, 30, 00))
swi_idx = argmin(abs.(swi_data[:epoch] .- time_slice))
flux_4d = swi_data[:diff_en_fluxes][swi_idx, :, :, :]
flux_phi = dropdims(mean(flux_4d, dims = 2), dims = 2)'
flux_theta = dropdims(mean(flux_4d, dims = 1), dims = 1)'

