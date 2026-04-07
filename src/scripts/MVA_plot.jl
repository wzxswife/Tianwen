using Dates
include("MAVEN_load.jl")
include("TW_load.jl")
include("MAG_plot.jl")
using .MAVEN_load
using .TW_load
using .MAG_plot

# read data 
const Rm = 3390.0
date = Date(2022, 8, 24)
datestr = Dates.format.(date, "yyyymmdd")
dir = joinpath(@__DIR__, "..", "..")
data_path = joinpath(dir, "data", "MAVEN")
out_path = joinpath(dir, "Results", "MAVEN", "OWWPI")
mag_file = "mvn_mag_l2_2022236ss_" * datestr * "_v01_r01.sts"
mag_data = load_mag_l2(joinpath(data_path, mag_file))
mag_data[:position] .= mag_data[:position] ./ Rm
println("Data loaded.")

time_range = DateTime(date, Time(7, 30, 00)) .+ Dates.Second(10) .* range(0, 4)
bwave = mag_wave(mag_data, time_range)

using CairoMakie, GeometryBasics, LaTeXStrings
CairoMakie.activate!()
println("Drawing...")

size_inches = (15, 5)
size_pt = 72 .* size_inches
fig = Figure(size = size_pt, fontsize = 25, font = "Times New Roman")
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])
ax3 = Axis(fig[1, 3])
plot_MVA(ax1, ax2, ax3, bwave; smooth_window=4)
save("MVA_plot.png", fig)

size_inches = (30, 14)
size_pt = 72 .* size_inches
fig2 = Figure(size = size_pt, fontsize = 25, font = "Times New Roman")
ax21 = Axis(fig2[1, 1])
ax22 = Axis(fig2[2, 1])
ax23 = Axis(fig2[3, 1])
plot_MVA_time(ax21, ax22, ax23, mag_data, time_range)
Legend(fig2[1,2],ax21)
Legend(fig2[2,2],ax22)
Legend(fig2[3,2],ax23)
save("MVA_plot_time.png", fig2)
