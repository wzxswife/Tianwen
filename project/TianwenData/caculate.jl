using DelimitedFiles, DataFrames, EzXML #read txt  and xml and then convert them to DataFrames
using Dates
using LaTeXStrings
#using ContinuousWavelets  #wavelet power
using DSP
using Interpolations
using LinearAlgebra
using Wavelets
using FFTW
using Statistics
using CSV
include("TCWavelet.jl")

df = CSV.read("E:/work/Tianwen/Picture/1Hz Waves.csv", DataFrame)

const Rm = 3390.0  #km 
datapath2c32hz = "E:/work/Tianwen/Data/32Hz/"
for i=1:(size(df)[1])
    # read data
    date = df[i, :DateTime]
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
    ds = round(Int, df[1, :Duration]/6)
    timeWave = date .+ Dates.Second(ds) .* range(0,6)
    timeForeShock = date .- Dates.Minute(5) .+ Dates.Second(50) .* range(0,6)
    timeShock = date .- Dates.Minute(3) .+ Dates.Second(60 + ds) .* range(0,6)
    timeSolar = date .- Dates.Minute(5) .+ Dates.Second(50 + ds) .* range(0,6)
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
    indForeShock = findall(
    (minimum(timeForeShock) .<= mag2c32hz.Time .<= maximum(timeForeShock)) .&
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
    magForeShock = mag2c32hz[indForeShock, [:X_MSO, :Y_MSO, :Z_MSO]]
    magSolar = mag2c32hz[indSolar, [:X_MSO, :Y_MSO, :Z_MSO]]
    magWaveMean = zeros(1, 3)
    magSolarMean = zeros(1, 3)
    magForeShockMean = zeros(1, 3)
    for i in 1:3
        magWaveMean[1, i] = mean(magWave[:, i])
        magSolarMean[1, i] = mean(magSolar[:, i])
        magForeShockMean[1, i] = mean(magForeShock[:, i])
    end
    magWave = magWave .- magWaveMean
    #FFT计算波动频率
    # 计算单分量信号的功率谱密度（PSD）
    function compute_psd(signal::Vector{Float64}, fs::Float64)
        n = length(signal)
        # 1. 去均值（消除直流分量）
        signal_centered = signal .- mean(signal)
        # 2. 加汉宁窗（减少频谱泄漏）
        window = hanning(n)
        signal_windowed = signal_centered .* window
        # 3. 计算FFT（实数信号，只需正频率部分）
        fft_result = rfft(signal_windowed)
        freqs = rfftfreq(n, fs)  # 正频率轴
        # 4. 计算功率谱密度（PSD）
        psd = abs2.(fft_result) ./ (fs * sum(window.^2))  # 归一化
        return (frequencies=freqs, psd=psd)
    end

    # 主函数：分析磁场波动的主频率（针对三分量）
    function analyze_dominant_frequency(mag2c32hz::DataFrame, fs::Float64)
        
        mag_event = mag2c32hz[:, [:X_MSO, :Y_MSO, :Z_MSO]]
        
        # 为每个分量计算PSD和主频率
        dominant_freqs = Vector{Float64}(undef, 3)
        psd_results = Vector{Any}(undef, 3)
        
        for i in 1:3
            component = mag_event[:, i]
            freqs, psd = compute_psd(component, fs)
            
            max_idx = argmax(psd)
            dominant_freq = freqs[max_idx]
            dominant_freqs[i] = dominant_freq
            psd_results[i] = (frequencies=freqs, psd=psd)
            
        end
        
        return (
            maxFreq=dominant_freqs[argmax(dominant_freqs)],  # 最大主频率
            dominant=dominant_freqs,  # 各分量的主频率 [Bx, By, Bz]
            psd=psd_results          # 各分量的PSD数据（绘图用）
        )
    end

    fs = 32.0  
    freq = analyze_dominant_frequency(magWave, fs)
    temp = sqrt.(sum(Array(magWave).^2, dims=2))
    freq_main,psd_main = compute_psd(temp[:], fs)
    freq_main = freq_main[argmax(psd_main)]

    # 波动性质和激波性质
    function MVA(magbc)
        nb = length(magbc[:, 1])
        bm = mean(magbc, dims=1)
        muv = magbc' * magbc ./ nb - bm' * bm
        return eigen(muv)
    end

    function xyz_to_new_basis(point::Matrix{Float64}, basis_vectors::AbstractVector{<:AbstractVector})::Matrix{Float64}
        basis_matrix = hcat(basis_vectors...)
        return (basis_matrix' * point')'
    end

    function angle_between(v1, v2)
        cosθ = dot(v1, v2) / (norm(v1) * norm(v2))
        return acosd(clamp(cosθ, -1.0, 1.0))  # 限制范围防止数值误差
    end

    # 波动性质
    bWave = Array(magWave)
    Bmean = Array(magWaveMean)
    bm = MVA(bWave)
    # 波矢量WaveNormal
    kWave = bm.vectors[:, argmin(bm.values)]
    Bmva = bWave*bm.vectors
    Bmvamean= dropdims(mean(Bmva, dims=1), dims=1)
    # 振幅Amplitude
    Amplitude = sum(bm.values[1:3])

    # 激波性质
    bShock = Array(magShock[:, 1:3])
    bSolar = Array(magSolar[:, 1:3])
    BmeanShock = dropdims(mean(bShock, dims=1), dims=1)
    bmShock = MVA(bShock)
    # 激波法向ShockNormal
    kShock = bmShock.vectors[:, argmin(bmShock.values)]
    # kShock = bmShock.vectors[:, 1]
    # 传播角PropagationAngle
    ProAngle = angle_between(kWave, magWaveMean)    # 波传播方向与背景磁场的夹角
    ShockAngle = angle_between(kShock, magForeShockMean)    # 激波法向与太阳风的夹角
    magback = sqrt(sum(Bmean.^2))
    maxAm = maximum([maximum(bWave[:,1])-minimum(bWave[:,1]), maximum(bWave[:,2])-minimum(bWave[:,2]), 
        maximum(bWave[:,3])-minimum(bWave[:,3])])

    # 记录数据
    println("写入数据...")
    dfStatic = DataFrame(
        DateTime=timeWave[begin],
        Duration=ds*6,
        Amplitude=round(Amplitude,digits = 3),
        Frequence=round(freq.maxFreq,digits = 3),
        ProAngle=round(ProAngle,digits = 3), 
        ShockAngle=round(ShockAngle,digits = 3),
        Fx=round(freq.dominant[1],digits = 3),
        Fy=round(freq.dominant[2],digits = 3),
        Fz=round(freq.dominant[3],digits = 3),
        magback=round(magback,digits = 3),
        maxAmplitude=round(maxAm,digits = 3)
    )
    CSV.write("E:/work/Tianwen/Picture/test.csv", dfStatic, append=true, header=false)
end