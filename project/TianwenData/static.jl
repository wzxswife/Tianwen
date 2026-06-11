#using Glob
using DelimitedFiles, DataFrames, EzXML #read txt  and xml and then convert them to DataFrames
using Dates
#using LaTeXStrings
#using ContinuousWavelets  #wavelet power
#using Interpolations
using Statistics
using CSV

df = CSV.read("/Users/admin/WorkSpace/JuliaCode/Tianwen/Picture/1HzWaves.csv", DataFrame)

data = copy(df)
ind1 = findall(x -> x>90, data[:, :ProAngle])
data[ind1,:ProAngle] = 180 .- data[ind1,:ProAngle]
ind2 = findall(x -> x>90, data[:, :ShockAngle])
data[ind2,:ShockAngle] = 180 .- data[ind2,:ShockAngle]
data[:,:ShockAngle] = 90 .- data[:,:ShockAngle]

@show count((i -> i <= 45), data[:, :ShockAngle]) #count the number of ShockAngle <= 45

using CairoMakie
CairoMakie.activate!()

fig1 = Figure(resolution = (800, 600), fontsize = 20, font = "Times New Roman")
ax1 = Axis(fig1[1, 1], xlabel = "Magnetic" * " " * "Field(nT)", ylabel = "Frequency(Hz)")
scatter!(ax1, data[:, :magback], data[:, :Frequence], color = "turquoise3", markersize = 18, label = "Data")
x = [0.0: 0.5: 6;]
#y = b .+ k .* x
y= 0.196 .* x .+ 0.33
str1 = "f(Hz) =0.33 + 0.197B(nT)"
lines!(ax1, x, y, color = "coral3", linewidth = 3, label = "Linear Fit")
axislegend()
text!(ax1, 0.8, 0.9-0.1, text = str1, font = :bold, align = (:center, :center), 
    space = :relative, fontsize = 20, color=:black, colorrange=(1, 7), colormap=:tab10)
save("./Picture/static/MtoF" * ".png", fig1)

ind2 = findall(data.Frequence.<0.5 .|| data.magback.>7.5)
tempdata = select(data, :magback, :Frequence)
delete!(tempdata,ind2)
using GLM
lmr = lm(@formula(Frequence ~ magback), tempdata)
sqrt(r2(lmr))
using Statistics, StatsBase
StatsBase.cor(data.Frequence, data.magback)
nullmod = lm(@formula(Frequence ~ 1), tempdata)
str = "f(Hz) = " * string(round(coef(lmr)[1],digits=3)) * " + " * string(round(coef(lmr)[2],digits=3)) * "B(nT)"
mod1_predict_confidence = predict(lmr, 
    hcat(ones(nrow(tempdata)),collect(tempdata[:,1])),
    :confint)
x1 = collect(tempdata[:,1])
y1 = collect(tempdata[:,2])
y2 = mod1_predict_confidence.lower
y3 = mod1_predict_confidence.upper
yp = coef(lmr)[1] .+ coef(lmr)[2] .* x1

fig = Figure(resolution = (800, 600), fontsize = 20, font = "Times New Roman")
ax = Axis(fig[1, 1], xlabel = "Magnetic" * " " * "Field(nT)", ylabel = "Frequency(Hz)")
scatter!(ax, x1, y1, color = "turquoise3", markersize = 18, label = "Data")
lines!(ax, x1, yp, color = "coral3", linewidth = 3, label = "Linear Fit")
scatter!(ax, x1, y2, color = "purple" , markersize = 15, label = "Confidence Interval")
scatter!(ax, x1, y3, color = "purple", markersize = 15)
axislegend()
str = "f(Hz) = " * string(round(coef(lmr)[1],digits=3)) * " + " * string(round(coef(lmr)[2],digits=3)) * "B(nT)"
text!(ax, 0.8, 0.75, text = str, font = :bold, align = (:center, :center), 
    space = :relative, fontsize = 20, color=:black, colorrange=(1, 7), colormap=:tab10)
save("./Picture/static/MtoF2" * ".png", fig)
using LsqFit

@.model(m, p) = p[1] + p[2] * m

mdata = tempdata[:, :magback]
fdata = tempdata[:, :Frequence]
p0=[0.5,0.5]
fit=curve_fit(model, mdata, fdata, p0)
p = fit.param
fig2 = Figure(resolution = (800, 600), fontsize = 20, font = "Times New Roman")
ax2 = Axis(fig2[1, 1], xlabel = "Propagation" * " " * "Angle(deg)", ylabel = "Amplitude(Hz)")
scatter!(ax2, data[:, :ProAngle], data[:, :Amplitude])
save("./Picture/static/AntoAm" * ".png", fig2)

fig3 = Figure(resolution = (800, 600), fontsize = 20, font = "Times New Roman")
ax3 = Axis(fig3[1, 1], xlabel = "ShockAngle(deg)", ylabel = "Density")
#hist!(ax3, data[:, :ShockAngle], bins = 20, color = :gray74, strokewidth = 1, strokecolor = :black)
density!(ax3, data[:, :ShockAngle], color = (:gray60, 0.3), strokecolor = :gray60, strokewidth = 3, strokearound = true)
xlims!(ax3, 0, 90)
save("./Picture/static/ShockAngle" * ".png", fig3)

fig4 = Figure(resolution = (800, 600), fontsize = 20, font = "Times New Roman")
ax4 = Axis(fig4[1, 1], xlabel = "Frequency(Hz)", ylabel = "Density")
#hist!(ax4, data[:, :Frequence], bins = [0.5:0.1:1.6;], color = :gray60, strokewidth = 1, strokecolor = :black)
density!(ax4, data[:, :Frequence], color = (:gray, 0.3), strokecolor = :gray, strokewidth = 3, strokearound = true)
xlims!(ax4, 0, 2.0)
ylims!(ax4, 0, 1.5)
save("./Picture/static/Fre1" * ".png", fig4)

using Glob  # find files with a specific pattern
using DelimitedFiles, DataFrames, EzXML #read txt  and xml and then convert them to DataFrames
using Dates
using LaTeXStrings
#using ContinuousWavelets  #wavelet power
using Interpolations # interpolation
#using StaticArrays #faster array
#using VMD # VMD decomposition
using DSP
using LinearAlgebra
using Wavelets
using FFTW
using Statistics
using CSV
include("TCWavelet.jl")

date = DateTime(2021, 12, 1, 11, 39, 0)
ds = 10

# read data
const Rm = 3390.0  #km 
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
mag2c32hz[!,[ :Probe_Position_X_MSO, :Probe_Position_Y_MSO, :Probe_Position_Z_MSO]] = 
    mag2c32hz[!,[ :Probe_Position_X_MSO, :Probe_Position_Y_MSO, :Probe_Position_Z_MSO]]./Rm
println("Data read finished!")
timeWave = date .+ Dates.Second(ds) .* range(0,6)

ind = findall(
    (minimum(timeWave) .<= mag2c32hz.Time .<= maximum(timeWave)) .&
    .!isnan.(mag2c32hz.X_MSO) .&
    .!isnan.(mag2c32hz.Y_MSO) .&
    .!isnan.(mag2c32hz.Z_MSO)
)
ReB = BMSO32hz[ind, [:X_MSO, :Y_MSO, :Z_MSO]]
ReB[:, :MSO] = sqrt.(ReB[:, 1].^2 .+ ReB[:, 2].^2 .+ ReB[:, 3].^2)

#wavelet power
ns32hz = size(ReB)[1]
dt32hz = 1.0/32.0

println("Wavelet power calculation...")
mother = "MORLET"
wave32hz, period32hz, scale32hz, coi32hz = wavelet(reshape(ReB[:, 1], ns32hz), dt32hz; pad=1, mother=mother)
xpower32hz = abs.(wave32hz).^2
wave32hz, period32hz, scale32hz, coi32hz = wavelet(reshape(ReB[:, 2], ns32hz), dt32hz; pad=1, mother=mother)
ypower32hz = abs.(wave32hz).^2
wave32hz, period32hz, scale32hz, coi32hz = wavelet(reshape(ReB[:, 3], ns32hz), dt32hz; pad=1, mother=mother)
zpower32hz = abs.(wave32hz).^2
Bpower32hz =xpower32hz+ypower32hz+zpower32hz
println("Wavelet power calculation finished!")

# Bpower = Bpower32hz[:, size(Bpower32hz)[2]]
Bpower = sum(Bpower32hz, dims=2)./size(Bpower32hz)[2]
fig = Figure(resolution = (600, 600), fontsize = 18, font = "Times New Roman")
ax = Axis(fig[1, 1], xscale = log2, yscale = log10, 
    xminorticksvisible = true, xminorgridvisible = true, xminorticks = IntervalsBetween(5), 
    yminorticksvisible = true, yminorgridvisible = true, yminorticks = IntervalsBetween(5), 
    xlabel = "Frequency(Hz)", ylabel = "Wavelet Power", title = "Wavelet Power of Magnetic Field")
lines!(ax, 1 ./period32hz, Array(Bpower)[:, 1])
save("./Picture/static/Bpower" * ".png", fig)
fig4 = Figure(resolution = (800, 600), fontsize = 20, font = "Times New Roman")
ax4 = Axis(fig4[1, 1], xlabel = "Amplitude(nT)", ylabel = "Number")
hist!(ax4, data[:, :Amplitude], bins = [0: 1: 10;], color = :gray60, strokewidth = 1, strokecolor = :black)
#density!(ax4, data[:, :Amplitude], color = (:gray, 0.3), strokecolor = :gray, strokewidth = 3, strokearound = true)
xlims!(ax4, 0, 10)
ylims!(ax4, 0, 20)
save("./Picture/static/Amp" * ".png", fig4)

using CairoMakie
CairoMakie.activate!()

fig5 = Figure(resolution = (800, 600), fontsize = 20, font = "Times New Roman")
ax5 = Axis(fig5[1, 1],  xlabel = "Frequency(Hz)", ylabel = "Number")
hist!(ax5, data[:, :Frequence], bins = [0.5:0.1:1.6;], color = :gray60, strokewidth = 1, strokecolor = :black)
xlims!(ax5, 0, 10)
ylims!(ax4, 0, 20)
ax6 = Axis(fig5[1, 2],  xlabel = "Frequency(Hz)", ylabel = "Density")
density!(ax6, data[:, :Frequence], color = (:gray, 0.3), strokecolor = :gray, strokewidth = 3, strokearound = true)
xlims!(ax4, 0, 10)
ylims!(ax4, 0, 20)
ax7 = Axis(fig5[2, 1], xlabel = "Amplitude(nT)", ylabel = "Number")
hist!(ax7, data[:, :Amplitude], bins = [0: 1: 10;], color = :gray60, strokewidth = 1, strokecolor = :black)
xlims!(ax4, 0, 10)
ylims!(ax4, 0, 20)
ax8 = Axis(fig5[2, 2], xlabel = "Amplitude(nT)", ylabel = "Density")
density!(ax8, data[:, :Amplitude], color = (:gray, 0.3), strokecolor = :gray, strokewidth = 3, strokearound = true)
xlims!(ax4, 0, 10)
ylims!(ax4, 0, 20)
save("./Picture/static/ALL" * ".png", fig5)

using CairoMakie
CairoMakie.activate!()
fig6 = Figure(resolution = (800, 600), fontsize = 20, font = "Times New Roman")
ax10 = Axis(fig6[1, 1], xlabel = "Time", ylabel = "Amplitude(nT)")
scatter!(ax10, data[:, :DateTime], data[:, :Amplitude], color = "turquoise3", markersize = 10, label = "Data")
save("./Picture/static/TimeAmp" * ".png", fig6)