using Dates
using DataFrames
using Statistics
using LinearAlgebra

include("../../src/scripts/TW_load.jl")
include("../../src/scripts/MAG_plot.jl")
include("../../src/scripts/wave_caculate.jl")
using .TW_load
using .MAG_plot
using .WaveCaculate

const Rm = 3389.5 # 火星半径，单位km

date = Date(2022, 8, 24)
datestr = Dates.format.(date, "yyyymmdd")

data_path = joinpath(@__DIR__, "..", "..", "data", "32Hz")
file_name = "TW1_MOMAG_MSO_32Hz_" * datestr * "_2C_v03.dat"
mag2c32hz = load_mag_2c_bydlm(joinpath(data_path, file_name))
mag2c32hz[:position] .= mag2c32hz[:position] ./ Rm
println("Data loaded.")

start_time = DateTime(date, Time(7, 10, 0))
time_range = start_time .+ Dates.Minute(10) .* range(0, 5)

mag_data = find_avail_data(mag2c32hz, time_range, [:epoch, :JulUTtime, :B])
dt = 1.0/32
Bpower, period = caculate_wavelet(mag_data[:B], dt)
time_data = mag_data[:JulUTtime]

using CairoMakie, GeometryBasics
CairoMakie.activate!()
println("Drawing...")
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
plot_spacecraft_orbit(ax1, ax2, ax3, mag2c32hz, [time_range[1], time_range[end]])

ax4 = Axis(fig[2, 1:3], xlabel = "Time (UT)", 
    ylabel = L"Magnetic ${\mathrm B_{total}}$ (nT)")
plot_B_total(ax4, mag2c32hz, time_range)
Legend(fig[2, 4], ax4)

ax5 = Axis(fig[3, 1:3], xlabel = "Time (UT)", 
    ylabel = L"Magnetic $\mathbf{B}$ (nT)")
plot_B(ax5, mag2c32hz, time_range)
Legend(fig[3, 4], ax5)

plot_wavelet(fig, 4, time_data, time_range, Bpower, period, dt)
save(joinpath(@__DIR__, "TW1_32Hz_wavelet.png"), fig)