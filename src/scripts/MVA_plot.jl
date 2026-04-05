module MVA_plot
using Dates
using Statistics
using LinearAlgebra
include("MAVEN_load.jl")
include("TW_load.jl")
using .MAVEN_load
using .TW_load

export MVA, plot_MVA, mag_wave, smooth_data

function mag_wave(mag_data, time_range)
    mag_wave = find_avail_data(mag_data, time_range, [:B])
    mag_wave = convert.(Float64, mag_wave[:B])
    bmean = mean(mag_wave, dims=1)
    bwave = mag_wave .- bmean
    return bwave
end

function smooth_data(data::Matrix{Float64}, window::Int)
    n = size(data, 1)
    smoothed = similar(data)
    half_window = window ÷ 2
    for i in 1:n
        start_idx = max(1, i - half_window)
        end_idx = min(n, i + half_window)
        smoothed[i, :] = mean(data[start_idx:end_idx, :], dims=1)
    end
    return smoothed
end

function MVA(magbc)
    nb = length(magbc[:, 1])
    bm = mean(magbc, dims=1)
    muv = magbc' * magbc ./ nb - bm' * bm
    return eigen(muv)
end

function plot_MVA(ax1, ax2, ax3, mag_wave::Matrix{Float64}; 
    color=(0.2, 0.3, 0.7), lw=1.5, ms=20, smooth_window=10)
    bm = MVA(mag_wave)
    min_hat = bm.vectors[:, 1]
    mid_hat = bm.vectors[:, 2]
    max_hat = bm.vectors[:, 3]
    R = [max_hat mid_hat min_hat]
    bWaveMVA = mag_wave * R 
    if smooth_window > 1
        bWaveMVA = smooth_data(bWaveMVA, smooth_window)
    end
    lims_max = ceil(maximum(abs.(bWaveMVA)))
    lims_min = -lims_max
    ax1.xlabel = "B min"
    ax1.ylabel = "B mid"
    ax2.xlabel = "B min"
    ax2.ylabel = "B max"
    ax3.xlabel = "B mid"
    ax3.ylabel = "B max"

    limits!(ax1, lims_min, lims_max, lims_min, lims_max)
    limits!(ax2, lims_min, lims_max, lims_min, lims_max)
    limits!(ax3, lims_min, lims_max, lims_min, lims_max)
    lines!(ax1, bWaveMVA[:, 3], bWaveMVA[:, 2], color=color, linewidth=lw)
    lines!(ax2, bWaveMVA[:, 3], bWaveMVA[:, 1], color=color, linewidth=lw)
    lines!(ax3, bWaveMVA[:, 2], bWaveMVA[:, 1], color=color, linewidth=lw)
    scatter!(ax1, [bWaveMVA[1, 3]], [bWaveMVA[1, 2]], 
        color=:lightgreen, markersize=ms, marker=:circle, strokewidth=0)
    scatter!(ax1, [bWaveMVA[end, 3]], [bWaveMVA[end, 2]], 
        color=:darkred, markersize=ms, marker='x')
    scatter!(ax2, [bWaveMVA[1, 3]], [bWaveMVA[1, 1]], 
        color=:lightgreen, markersize=ms, marker=:circle, strokewidth=0)
    scatter!(ax2, [bWaveMVA[end, 3]], [bWaveMVA[end, 1]], 
        color=:darkred, markersize=ms, marker='x')
    scatter!(ax3, [bWaveMVA[1, 2]], [bWaveMVA[1, 1]], 
        color=:lightgreen, markersize=ms, marker=:circle, strokewidth=0)
    scatter!(ax3, [bWaveMVA[end, 2]], [bWaveMVA[end, 1]], 
        color=:darkred, markersize=ms, marker='x')
end

# -----------以下为测试脚本部分----------------
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
plot_MVA(ax1, ax2, ax3, bwave)

save("MVA_plot.png", fig)



end