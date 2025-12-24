# 用来画极化图的程序
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
using Plots
include("TCWavelet.jl")

const Rm = 3390.0  #km 
df = CSV.read("E:/work/Tianwen/Picture/1Hz Waves.csv", DataFrame)

function xyz_to_new_basis(point::Matrix{Float64}, basis_vectors::AbstractVector{<:AbstractVector})::Matrix{Float64}
    basis_matrix = hcat(basis_vectors...)
    return (basis_matrix' * point')'
end

for i = 1:size(df)[1]
    date = df[i, :DateTime]
    datestr = Dates.format.(date, "yyyymmdd")
    datapath2c32hz = "E:/work/Tianwen/Data/32Hz/"
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

    # 分离背景磁场和波动磁场
    ds = df[i, :Duration]/6
    timeWave = date .+ Dates.Second(ds) .* range(0,6)
    timeShock = date .- Dates.Minute(5) .+ Dates.Second(100+ds) .* range(0,6)
    timeSolar = date .+ Dates.Second(ds) .* range(0,6)
    indWave = findall(
    (minimum(timeWave) .<= mag2c32hz.Time .<= maximum(timeWave)) .&
    .!isnan.(mag2c32hz.X_MSO) .&
    .!isnan.(mag2c32hz.Y_MSO) .&
    .!isnan.(mag2c32hz.Z_MSO)
    )
    indShock = findall(
    (minimum(timeShock) .<= mag2c32hz.Time .<= maximum(timeShock)) .&
    .!isnan.(mag2c32hz.X_MSO) .&
    .!isnan.(mag2c32hz.Y_MSO) .&
    .!isnan.(mag2c32hz.Z_MSO)
    )
    indSolar = findall(
    (minimum(timeSolar) .<= mag2c32hz.Time .<= maximum(timeSolar)) .&
    .!isnan.(mag2c32hz.X_MSO) .&
    .!isnan.(mag2c32hz.Y_MSO) .&
    .!isnan.(mag2c32hz.Z_MSO)
    )
    magWave = mag2c32hz[indWave, [:X_MSO, :Y_MSO, :Z_MSO]]
    magShock = mag2c32hz[indShock, [:X_MSO, :Y_MSO, :Z_MSO]]
    magSolar = mag2c32hz[indSolar, [:X_MSO, :Y_MSO, :Z_MSO]]
    magWaveMean = zeros(1, 3)
    magSolarMean = zeros(1, 3)
    for i in 1:3
        magWaveMean[1, i] = mean(magWave[:, i])
        magSolarMean[1, i] = mean(magSolar[:, i])
    end
    magWave = magWave .- magWaveMean
    bWave = Array(magWave)
    Bmean = Array(magWaveMean)

    # 极化角度PolarAngle
    e1 = (Bmean./sqrt(sum(Bmean.* Bmean)))[1, :]
    b1 = bWave[1, :] ./ sqrt(sum(bWave[1, :].^ 2))
    e3 = cross(e1, b1)
    e3 = e3./sqrt(sum(e3.^2))
    e2 = cross(e3, e1)
    e2 = e2./sqrt(sum(e2.^2))
    bwave_back = xyz_to_new_basis(bWave, [e1, e2, e3])

    # 画极化图（有箭头）
    len = size(magWave)[1]
    # 绘制带箭头的折线图
    fig1 = Plots.plot(magWave[1:len-100, 2], magWave[1:len-100, 3], arrow=(Plots.arrow(:closed,:head,5.,2.)), linewidth = 0.5, color=:black, 
        xlabel=L"$B_M$", ylabel=L"$B_L$", dpi = 300, legend=:none)
    fig1 = Plots.plot!(magWave[len-100+1:len, 2], magWave[len-100+1:len, 3], arrow=(Plots.arrow(:closed,:head,5.,2.)), linewidth = 0.5, color=:red, 
        dpi = 300, legend=:none)

    savefig(fig1, "./Picture/phase/" * Dates.format(timeWave[1], "yyyymmdd HH-MM-SS") * ".png")
end