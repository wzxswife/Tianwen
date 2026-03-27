using Dates
using DataFrames
using Statistics
using LinearAlgebra

include("./scripts/TW_load.jl")
include("./scripts/TW_plot.jl")
include("./scripts/wave_caculate.jl")
using .TW_load
using .TW_plot
using .WaveCaculate

data_path = joinpath(@__DIR__, "..", "data", "32Hz")
file_name = "TW1_MOMAG_MSO_32Hz_20220824_2C_v03.dat"
mag2c32hz = load_mag_2c_bydlm(joinpath(data_path, file_name))
mag2c32hz[:position] .= mag2c32hz[:position] ./ 3390.0
println("Data loaded.")

data = mag2c32hz
shock_time = DateTime("2022-08-24T07:27:20")
time_range = shock_time - Dates.Minute(1) .+ Dates.Minute(1) .* range(0, 3)

mag_data = find_avail_data(data, time_range, [:epoch, :JulUTtime, :B])
dt = 1.0/32
Bpower, period = caculate_wavelet(mag_data[:B], dt)
time_data = mag_data[:JulUTtime]

using CairoMakie, GeometryBasics
CairoMakie.activate!()
println("Drawing...")
fig1 = Figure(resolution = (800, 400))
plot_wavelet(fig1, 1, time_data, time_range, Bpower, period, dt)
save(joinpath(@__DIR__, "TW1_32Hz_wavelet.png"), fig1)

println("Drawing...")
fig2 = Figure()
ax2 = Axis(fig2[1, 1:3], title = "TW1 32Hz Magnetic Field", xlabel = "Time (UTC)", ylabel = "B (nT)")
plot_B(ax2, data, time_range)
Legend(fig2[1, 4], ax2)
save(joinpath(@__DIR__, "TW1_32Hz_B.png"), fig2)

println("Drawing...")
fig3 = Figure()
ax3 = Axis(fig3[1, 1], title = "TW1 32Hz Magnetic Field", xlabel = "Time (UTC)", ylabel = "B (nT)")
plot_B_total(ax3, data, time_range)
save(joinpath(@__DIR__, "TW1_32Hz_B_total.png"), fig3)

println("Drawing...")
fig4 = Figure(resolution=(800, 1200))
ax41 = Axis(fig4[1,1], aspect=DataAspect(), xlabel=L"$x$ ($R_{\mathrm{M}}$)", 
    ylabel=L"$\sqrt{y^2+z^2}$ ($R_{\mathrm{M}}$)", 
    xticks=range(-10, 10), yticks=range(-10, 10))
ax42 = Axis(fig4[2,1], aspect=DataAspect(), xlabel=L"$x$ ($R_{\mathrm{M}}$)", 
    ylabel=L"$\sqrt{y^2+z^2}$ ($R_{\mathrm{M}}$)", 
    xticks=range(-10, 10), yticks=range(-10, 10))
ax43 = Axis(fig4[3,1], aspect=DataAspect(), xlabel=L"$x$ ($R_{\mathrm{M}}$)", 
    ylabel=L"$\sqrt{y^2+z^2}$ ($R_{\mathrm{M}}$)", 
    xticks=range(-10, 10), yticks=range(-10, 10))
bowshock_plot(ax41, ax42, ax43)
plot_spacecraft_orbit(ax41, ax42, ax43, data, time_range)
save("bowshock.png", fig4)