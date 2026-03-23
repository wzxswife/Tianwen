"""
波动性质计算模块
"""
module WaveCaculate
using Dates
using DSP
using LinearAlgebra
using Wavelets
using FFTW
include("TCWavelet.jl")
using Wavelets

export MVA, analyze_dominant_frequency, xyz_to_new_basis, angle_between

# MVA分析，计算特征值和特征向量
function MVA(magbc)
    nb = length(magbc[:, 1])
    bm = mean(magbc, dims=1)
    muv = magbc' * magbc ./ nb - bm' * bm
    return eigen(muv)
end

# 计算单分量信号的功率谱密度（PSD）
function compute_psd(signal::Vector, fs::Float64)
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
function analyze_dominant_frequency(mag::Matrix, fs::Float64)
    mag_event = mag
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

# 简单的坐标变换函数（示例：从MSO到新基底）
function xyz_to_new_basis(point::Matrix{Float64}, basis_vectors::AbstractVector{<:AbstractVector})::Matrix{Float64}
    basis_matrix = hcat(basis_vectors...)
    return (basis_matrix' * point')'
end
# 计算两个向量之间的夹角（单位：度）
function angle_between(v1, v2)
    cosθ = dot(v1, v2) / (norm(v1) * norm(v2))
    return acosd(clamp(cosθ, -1.0, 1.0))  # 限制范围防止数值误差
end


end