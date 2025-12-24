using DelimitedFiles, DataFrames, EzXML #read txt  and xml and then convert them to DataFrames
using Dates
using LaTeXStrings
using LinearAlgebra
using Interpolations

ds = 60

date = DateTime(2021, 12, 15, 15, 10, 0)
# time = date .+ Dates.Hour(0) .+ Dates.Minute(0) .+ Dates.Second(0) .+ Dates.Second(ds) .* range(0,12)
time = date .+ Dates.Second(ds) .* range(0,10)
titles = Dates.format.(date, "yyyy-mm-dd") * " " * Dates.format.(time[1], "HH:MM:SS") * "-" * 
    Dates.format.(time[end], "HH:MM:SS") * "    " * "Martian 1Hz Waves"
# println(time)

# read data 
datapath2c32hz = "E:/Acode/JuliaCode/Tianwen/Data/32Hz/"
datestr = Dates.format.(date, "yyyymmdd")
file2c32hz = datapath2c32hz * "TW1_MOMAG_MSO_32Hz_" * datestr * "_2C_v03.dat"
println("Reading file: ")
println(file2c32hz)
global mag2c32hz = identity.(DataFrame(readdlm(file2c32hz, skipstart=19), :auto))
name = ["Time", "Sampling_Rate", "X_MSO", "Y_MSO", "Z_MSO", "Probe_Position_X_MSO", "Probe_Position_Y_MSO", 
    "Probe_Position_Z_MSO", "Roll", "Pitch", "Yaw",  "Quality_Flags"]
rename!(mag2c32hz, name)
mag2c32hz[!, :Time] = map(x->DateTime(x[begin:end-4], DateFormat("y-m-dTH:M:S.s")), mag2c32hz[!, :Time])
mag2c32hz[!, :JulUT] = datetime2julian.(mag2c32hz[!, :Time])
unique!(mag2c32hz) #remove depulicate rows
sort!(mag2c32hz) #sorting
magut32hz = mag2c32hz[:, :JulUT].-mag2c32hz[1, :JulUT]
global BMSO32hz = mag2c32hz[:, [:X_MSO, :Y_MSO, :Z_MSO]]
ind0 = findall(x-> !isnan(x), BMSO32hz[!, 1] ) 
ind1 = findall(x-> isnan(x), BMSO32hz[!, 1] ) 
if length(ind1) >= 1 && length(ind0) > 1
    for ib in 1:3 
        interp_linear = linear_interpolation(magut32hz[ind0], BMSO32hz[ind0, ib]; extrapolation_bc=Line())
        BMSO32hz[ind1, ib] = interp_linear(magut32hz[ind1])
    end
end
println("Data read finished!")

ind = findall(x-> x<=maximum(time)  && x>=minimum(time), mag2c32hz[!, :Time])
ReB = BMSO32hz[:, [:X_MSO, :Y_MSO, :Z_MSO]]
ReB[!, :X_MSO] = ReB[!, :X_MSO]
ReB[!, :Y_MSO] = ReB[!, :Y_MSO]
ReB[!, :Z_MSO] = ReB[!, :Z_MSO]
#ReB[!, :X_MSO] = ReB[!, :X_MSO] .- round(Int, mean(mag2c32hz[ind, :X_MSO]))
#ReB[!, :Y_MSO] = ReB[!, :Y_MSO] .- round(Int, mean(mag2c32hz[ind, :Y_MSO])) 
#ReB[!, :Z_MSO] = ReB[!, :Z_MSO] .- round(Int, mean(mag2c32hz[ind, :Z_MSO]))

#overview
using CairoMakie
CairoMakie.activate!()
println("Drawing...")
# light_latexfonts = merge(theme_minimal(), theme_latexfonts())
# set_theme!(light_latexfonts)

# size_inches = (30, 17)
# size_pt = 72 .* size_inches
# fig = Figure(size = size_pt, fontsize = 25)\
size_inches = (30, 17)
size_pt = 72 .* size_inches
fig = Figure(size = size_pt, fontsize = 25)
xtk = datetime2julian.(time)
xtk = Float32.(xtk .- mag2c32hz[1, :JulUT])
ax = Axis(fig[1,1], xlabel = L"Time(MM:SS)", ylabel = L"Magnetic $\mathbf{B}$ (nT)", xticks=xtk, title=titles, 
    xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(5), 
    yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(4))
ax.xticks = (xtk, Dates.format.(time, "MM:SS"))
# ax.xticks = (range(1, step=ds, length=10))
xlims!(ax, minimum(xtk), maximum(xtk))
ylims!(ax, minimum(Matrix(ReB[ind, :]))-1.1, maximum(Matrix(ReB[ind, :]))+1.1)

lines!(ax, mag2c32hz[!, :JulUT].-mag2c32hz[1, :JulUT], ReB[!, :X_MSO], label = L"$B_{\mathrm{x}}$", 
    overdraw = true, linewidth=1.5)
lines!(ax, mag2c32hz[!, :JulUT].-mag2c32hz[1, :JulUT], ReB[!, :Y_MSO], label = L"$B_{\mathrm{y}}$", 
    overdraw = true, linewidth=1.5)
lines!(ax, mag2c32hz[!, :JulUT].-mag2c32hz[1, :JulUT], ReB[!, :Z_MSO], label = L"$B_{\mathrm{z}}$", 
    overdraw = true, linewidth=1.5)
axislegend(ax)
println("Drawing finished!")

save("./Picture/1HzWaves/" * datestr * ".png", fig)