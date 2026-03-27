using Dates
using DataFrames
using Statistics
using LinearAlgebra

include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_plot.jl")
include("../../src/scripts/TW_load.jl")
include("../../src/scripts/TW_plot.jl")
include("../../src/scripts/wave_caculate.jl")
using .MAVEN_load
using .MAVEN_plot
using .TW_load
using .TW_plot
using .WaveCaculate

const Rm = 3390.0
date = Date(2022, 8, 24)
datestr = Dates.format.(date, "yyyymmdd")
dir = joinpath(@__DIR__, "..", "..")
data_path = joinpath(dir, "data", "MAVEN")
mag_file = "mvn_mag_l2_2022236ss_" * datestr * "_v01_r01.sts"
swi_file = "mvn_swi_l2_coarsesvy3d_" * datestr * "_v02_r01.cdf"

mag_data = load_mag_l2(joinpath(data_path, mag_file))
mag_data[:position] .= mag_data[:position] ./ Rm
swi_data = load_cdf(joinpath(data_path, swi_file))
println("Data loaded.")

start_time = DateTime(date, Time(6, 0, 0, 0))
time_range = start_time .+ Dates.Minute(30) .* range(0, 6)

using CairoMakie
CairoMakie.activate!()

fig = Figure(resolution = (1600, 1200))
ax1 = Axis(fig[1,1], aspect=DataAspect(), 
    xlabel=L"$x$ ($R_{\mathrm{M}}$)", 
    ylabel=L"$\sqrt{y^2+z^2}$ ($R_{\mathrm{M}}$)", 
    xticks=range(-10, 10), yticks=range(-10, 10))
ax2 = Axis(fig[1,2], aspect=DataAspect(), 
    xlabel=L"$x$ ($R_{\mathrm{M}}$)", ylabel=L"$y$ ($R_{\mathrm{M}}$)", 
    xticks=range(-10, 10, 11), yticks=range(-10, 10, 11))
ax3 = Axis(fig[1,3], aspect=DataAspect(), 
    xlabel=L"$x$ ($R_{\mathrm{M}}$)", ylabel=L"$z$ ($R_{\mathrm{M}}$)", 
    xticks=range(-10, 10, 11), yticks=range(-10, 10, 11))

bowshock_plot(ax1, ax2, ax3)
plot_spacecraft_orbit(ax1, ax2, ax3, mag_data, time_range)

ax4 = Axis(fig[2,1:3], xlabel="Time (UT)", ylabel="B_total (nT)")
plot_B_total(ax4, mag_data, time_range)
Legend(fig[2,4], ax4)

ax5 = Axis(fig[3,1:3], xlabel = "UT", ylabel = L"Magnetic $\mathbf{B}$ (nT)")
plot_B(ax5, mag_data, time_range)
Legend(fig[3,4], ax5)

mag_data[:B] = convert.(Float64, mag_data[:B])
avail_data = find_avail_data(mag_data, time_range, [:epoch, :JulUTtime, :B])
dt = 1.0/32
Bpower, period = caculate_wavelet(avail_data[:B], dt)
time_data = avail_data[:JulUTtime]
plot_wavelet(fig, 4, time_data, time_range, Bpower, period, dt)

swi_heatmap(fig, 5, swi_data, time_range)

save(joinpath(@__DIR__, "bowshock.png"), fig)
println("Plot saved.")