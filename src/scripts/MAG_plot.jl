# 绘制Tianwen相关图像

module MAG_plot
using Dates
using DelimitedFiles
using Interpolations
using Statistics
using LinearAlgebra
using DSP
using Wavelets
using CairoMakie
using GeometryBasics
using LaTeXStrings
include("TW_load.jl")
using .TW_load

export bowshock, magnetopause
export bowshock_plot_xr, plot_spacecraft_orbit_xr, 
    bowshock_plot, plot_spacecraft_orbit, 
    plot_B_total, plot_B, plot_wavelet
export mag_wave, smooth_data
export MVA, plot_MVA, plot_MVA_time


const Rm = 3389.5 # 火星半径，单位km

"""
    bowshock(x)
火星弓激波模型
参数:
    x: 火星固连坐标系X坐标 (Rm)
返回值:
    弓激波半径 (Rm), Inf64 表示超出模型范围
"""
function bowshock(xshock)
    xF = 0.55
    ϵ = 1.05
    L = 2.10
    temp = (ϵ^2-1.0)*(xshock-xF)^2 - 2ϵ*L*(xshock-xF) + L^2
    return temp >= 0 ? sqrt(temp) : Inf64
end

"""
    magnetopause(x)
火星磁层顶模型
参数:
    x: 火星固连坐标系X坐标 (Rm)
返回值:
    磁层顶半径 (Rm), Inf64 表示超出模型范围
"""
function magnetopause(xmp)
    rSD = 1.33
    xF = 0.86
    ϵ = 0.92
    L = 0.90
    temp = (ϵ^2-1.0)*(xmp-xF)^2 - 2ϵ*L*(xmp-xF) + L^2
    return temp >= 0 ? sqrt(temp) : Inf64
end

function mag_wave(mag_data::Dict, time_range::Vector{DateTime})
    local mag_wave = find_avail_data(mag_data, time_range, [:B])
    mag_wave = convert.(Float64, mag_wave[:B])
    local bmean = mean(mag_wave, dims=1)
    local bwave = mag_wave .- bmean
    return bwave
end

function smooth_data(data::Matrix{Float64}, window::Int)
    local n = size(data, 1)
    local smoothed = similar(data)
    local half_window = window ÷ 2
    for i in 1:n
        start_idx = max(1, i - half_window)
        end_idx = min(n, i + half_window)
        smoothed[i, :] = mean(data[start_idx:end_idx, :], dims=1)
    end
    return smoothed
end

function MVA(magbc)
    local nb = length(magbc[:, 1])
    local bm = mean(magbc, dims=1)
    local muv = magbc' * magbc ./ nb - bm' * bm
    return eigen(muv)
end

"""
画火星弓激波、磁层顶和火星位置（测试通过）
"""
# 只有x_r图
function bowshock_plot_xr(ax)
    np = 5001
    xshock = range(-10, 3, length=np)
    ryzshock = bowshock.(xshock)
    ind1 = findall(isfinite, ryzshock)
    xshock = [xshock[ind1]; reverse(xshock[ind1])]
    ryzshock = [ryzshock[ind1]; -reverse(ryzshock[ind1])]
    xmp = range(-10, 3, length=np)
    ryzmp = magnetopause.(xmp)
    ind1 = findall(isfinite, ryzmp)
    xmp = [xmp[ind1]; reverse(xmp[ind1])]
    ryzmp = [ryzmp[ind1]; -reverse(ryzmp[ind1])]

    xlims!(ax, -5, 5) 
    ylims!(ax, 0, 5) 
    arc!(ax, Point2f(0), 1, 0, π; color=:black, linewidth=2)
    phi = range(0, 0.5π; length=180)
    x = [cos.(phi); 0.0]
    y = [sin.(phi); 0.0]
    pn = Polygon(Point2f[(xi, yi) for (xi, yi) in zip(x, y)])
    poly!(ax, pn)

    lines!(ax, xshock, ryzshock; linewidth=2, overdraw = true)
    lines!(ax, xmp, ryzmp; linewidth=2, overdraw = true)
end
# x_r, x_y, x_z图
function bowshock_plot(ax1, ax2, ax3)
    np = 5001
    xshock = range(-10, 3, length=np)
    ryzshock = bowshock.(xshock)
    ind1 = findall(isfinite, ryzshock)
    xshock = [xshock[ind1]; reverse(xshock[ind1])]
    ryzshock = [ryzshock[ind1]; -reverse(ryzshock[ind1])]
    xmp = range(-10, 3, length=np)
    ryzmp = magnetopause.(xmp)
    ind1 = findall(isfinite, ryzmp)
    xmp = [xmp[ind1]; reverse(xmp[ind1])]
    ryzmp = [ryzmp[ind1]; -reverse(ryzmp[ind1])]

    xlims!(ax1, -5, 5) 
    ylims!(ax1, 0, 5) 
    a1 = arc!(ax1, Point2f(0), 1, 0, π; color=:black, linewidth=2)
    phi = range(0, 0.5π; length=180)
    x = [cos.(phi); 0.0]
    y = [sin.(phi); 0.0]
    pn = Polygon(Point2f[(xi, yi) for (xi, yi) in zip(x, y)])
    poly!(ax1, pn)
    l1 = lines!(ax1, xshock, ryzshock; linewidth=2, overdraw = true)
    l2 = lines!(ax1, xmp, ryzmp; linewidth=2, overdraw = true)

    xlims!(ax2, -10, 10) 
    ylims!(ax2, -5, 5) 
    a2 = arc!(ax2, Point2f(0), 1, 0, 2π; color=:black, linewidth=2)
    phi = range(-0.5π, 0.5π; length=180)
    x = [cos.(phi); 0.0]
    y = [sin.(phi); 0.0]
    pn = Polygon(Point2f[(xi, yi) for (xi, yi) in zip(x, y)])
    poly!(ax2, pn, color = :black)
    #l1 = lines!(ax2, xshock, ryzshock; linewidth=2)
    #l2 = lines!(ax2, xmp, ryzmp; linewidth=2)

    xlims!(ax3, -10, 10) 
    ylims!(ax3, -5, 5) 
    a3 = arc!(ax3, Point2f(0), 1, 0, 2π; color=:black, linewidth=2)
    phi = range(-0.5π, 0.5π; length=180)
    x = [cos.(phi); 0.0]
    y = [sin.(phi); 0.0]
    pn = Polygon(Point2f[(xi, yi) for (xi, yi) in zip(x, y)])
    poly!(ax3, pn, color = :black)
    #l1 = lines!(ax3, xshock, ryzshock; linewidth=2)
    #l2 = lines!(ax3, xmp, ryzmp; linewidth=2)
end

"""
画航天器轨道，和航天器位置（测试通过）
"""
function plot_spacecraft_orbit_xr(ax, data, time_range)
    xlims!(ax, -5, 5) 
    ylims!(ax, 0, 5)
    l3 = lines!(ax, data[:position][:, 1], 
        sqrt.(data[:position][:, 2].^2+data[:position][:, 3].^2); 
        color=:black, overdraw = true)
    len = length(time_range)
    for ti in 1:len
        local time = time_range[ti]
        local ind = findfirst(x-> x<time && x>time-Dates.Second(1), 
            data[:epoch])
        local str = Dates.format(time, "HH:MM")
        # println("Time: $time, Index: $ind /n")
        if ind !== nothing
            scatter!(ax, data[:position][ind, 1], 
                sqrt.(data[:position][ind, 2].^2+data[:position][ind, 3].^2); 
                color=ti, colorrange=(1, len), colormap=:tab10, 
                markersize=10)
            text!(ax, 0.9, 0.9-ti*0.1, text = str, font = :bold, 
                align = (:center, :center), space = :relative, fontsize = 25, 
                color=ti, colorrange=(1, len), colormap=:tab10)
        end
    end
end
function plot_spacecraft_orbit(ax1, ax2, ax3, data, time_range)
    xlims!(ax1, -5, 5) 
    ylims!(ax1, 0, 5)
    xlims!(ax2, -10, 10) 
    ylims!(ax2, -5, 5) 
    xlims!(ax3, -10, 10) 
    ylims!(ax3, -5, 5)
    lines!(ax1, data[:position][:, 1], 
        sqrt.(data[:position][:, 2].^2+data[:position][:, 3].^2); 
        color=:black, overdraw = true)
    lines!(ax2, data[:position][:, 1], data[:position][:, 2]; 
        color=:black, overdraw = true)
    lines!(ax3, data[:position][:, 1], data[:position][:, 3]; 
        color=:black, overdraw = true)
    len = length(time_range)
    for ti in 1:len
        local time = time_range[ti]
        local ind = findfirst(x-> x<time && x>time-Dates.Second(1), 
            data[:epoch])
        local str = Dates.format(time, "HH:MM")
        # println("Time: $time, Index: $ind /n")
        if ind !== nothing
            scatter!(ax1, data[:position][ind, 1], 
                sqrt.(data[:position][ind, 2].^2+data[:position][ind, 3].^2); 
                color=ti, colorrange=(1, len), colormap=:tab10, 
                markersize=10)
            text!(ax1, 0.9, 0.9-ti*0.1, text = str, font = :bold, 
                align = (:center, :center), space = :relative, fontsize = 24, 
                color=ti, colorrange=(1, len), colormap=:tab10)
            scatter!(ax2, data[:position][ind, 1], data[:position][ind, 2]; 
                color=ti, colorrange=(1, len), colormap=:tab10, 
                markersize=10)
            text!(ax2, 0.9, 0.9-ti*0.1, text = str, font = :bold, 
                align = (:center, :center), space = :relative, fontsize = 24, 
                color=ti, colorrange=(1, len), colormap=:tab10)
            scatter!(ax3, data[:position][ind, 1], data[:position][ind, 3]; 
                color=ti, colorrange=(1, len), colormap=:tab10, 
                markersize=10)
            text!(ax3, 0.9, 0.9-ti*0.1, text = str, font = :bold, 
                align = (:center, :center), space = :relative, fontsize = 24, 
                color=ti, colorrange=(1, len), colormap=:tab10)
        end
    end
end

"""
画磁场数据（测试通过）
"""
# 总磁场大小
function plot_B_total(ax, data, time_range)
    local lw = 0.8
    local dataB = find_avail_data(data, time_range, 
        [:epoch, :B_total, :JulUTtime])
    local xtk = datetime2julian.(time_range)
    local t0 = dataB[:JulUTtime][1]
    xtk = xtk .- t0
    xlims!(ax, xtk[1], xtk[end])
    ylims!(ax, minimum(dataB[:B_total])-1.0, maximum(dataB[:B_total])+1.0)
    lines!(ax, dataB[:JulUTtime].-t0, dataB[:B_total]; 
        label = L"$B_{\mathrm{total}}$", 
        color=:black, linewidth=lw, overdraw = true)
    ax.xticks = (xtk, Dates.format.(time_range, "HH:MM:SS"))
end
# 磁场三分量
function plot_B(ax, data, time_range)
    local lw = 0.8
    local dataB = find_avail_data(data, time_range, 
        [:epoch, :B, :JulUTtime])
    local xtk = datetime2julian.(time_range)
    local t0 = dataB[:JulUTtime][1]
    xtk = xtk .- t0
    xlims!(ax, xtk[1], xtk[end])
    ylims!(ax, minimum(dataB[:B])-1.0, maximum(dataB[:B])+1.0)
    lines!(ax, dataB[:JulUTtime].-t0, dataB[:B][:, 1]; 
        label = L"$B_{\mathrm{x}}$", linewidth=lw, overdraw = true)
    lines!(ax, dataB[:JulUTtime].-t0, dataB[:B][:, 2]; 
        label = L"$B_{\mathrm{y}}$", linewidth=lw, overdraw = true)
    lines!(ax, dataB[:JulUTtime].-t0, dataB[:B][:, 3]; 
        label = L"$B_{\mathrm{z}}$", linewidth=lw, overdraw = true)
    ax.xticks = (xtk, Dates.format.(time_range, "HH:MM:SS"))
end

"""
磁场wavelet频谱图（测试通过）
"""
function plot_wavelet(ax, data, time_range, dt;
    colorscale=log10, colormap=:gist_earth, colorrange=(2e-2, 3e2))
    local avail_data = find_avail_data(data, time_range, [:epoch, :JulUTtime, :B])
    local Bpower, period = caculate_wavelet(avail_data[:B], dt)
    local time_data = avail_data[:JulUTtime]
    local xtk = datetime2julian.(time_range)
    ax.xticks = (xtk, Dates.format.(time_range, "HH:MM:SS"))
    xlims!(ax, minimum(xtk), maximum(xtk)) 
    ylims!(ax, 2*dt, 256*dt)
    ax.yreversed = true
    local hp = heatmap!(ax, time_data, period, Bpower', 
        colorscale=colorscale, colormap=colormap, colorrange=colorrange)
    tightlimits!(ax)
    return hp
end

function plot_phase(ax, data::Dict, time_range; 
    lw=2, ms=20,smooth_window=10)
    local avail_data = find_avail_data(data, time_range, [:epoch, :JulUTtime, :B])
    if smooth_window > 1
        b_data = smooth_data(avail_data[:B], smooth_window)
    end
end

function plot_MVA(ax1, ax2, ax3, mag_wave::Matrix{Float64}; 
    color=(0.2, 0.3, 0.7), lw=2, ms=20, smooth_window=10)
    local bm = MVA(mag_wave)
    local min_hat = bm.vectors[:, 1]
    local mid_hat = bm.vectors[:, 2]
    local max_hat = bm.vectors[:, 3]
    local R = [max_hat mid_hat min_hat]
    local bWaveMVA = mag_wave * R 
    if smooth_window > 1
        bWaveMVA = smooth_data(bWaveMVA, smooth_window)
    end

    local lims_max = ceil(maximum(abs.(bWaveMVA)))
    local lims_min = -lims_max
    ax1.xlabel = "B min"
    ax1.ylabel = "B mid"
    ax2.xlabel = "B min"
    ax2.ylabel = "B max"
    ax3.xlabel = "B mid"
    ax3.ylabel = "B max"

    limits!(ax1, lims_min, lims_max, lims_min, lims_max)
    limits!(ax2, lims_min, lims_max, lims_min, lims_max)
    limits!(ax3, lims_min, lims_max, lims_min, lims_max)
    lines!(ax1, bWaveMVA[:, 3], bWaveMVA[:, 1], color=color, linewidth=lw)
    lines!(ax2, bWaveMVA[:, 3], bWaveMVA[:, 2], color=color, linewidth=lw)
    lines!(ax3, bWaveMVA[:, 1], bWaveMVA[:, 2], color=color, linewidth=lw)
    scatter!(ax1, [bWaveMVA[1, 3]], [bWaveMVA[1, 1]], 
        color=:lightgreen, markersize=ms, marker=:circle, strokewidth=0)
    scatter!(ax1, [bWaveMVA[end, 3]], [bWaveMVA[end, 1]], 
        color=:darkred, markersize=ms, marker='x')
    scatter!(ax2, [bWaveMVA[1, 3]], [bWaveMVA[1, 2]], 
        color=:lightgreen, markersize=ms, marker=:circle, strokewidth=0)
    scatter!(ax2, [bWaveMVA[end, 3]], [bWaveMVA[end, 1]], 
        color=:darkred, markersize=ms, marker='x')
    scatter!(ax3, [bWaveMVA[1, 1]], [bWaveMVA[1, 2]], 
        color=:lightgreen, markersize=ms, marker=:circle, strokewidth=0)
    scatter!(ax3, [bWaveMVA[end, 1]], [bWaveMVA[end, 2]], 
        color=:darkred, markersize=ms, marker='x')
end

function plot_MVA_time(ax1, ax2, ax3, data::Dict, time_range::Vector{DateTime};
    color=(0.2, 0.3, 0.7), lw=1.5)
    local mag_data = find_avail_data(data, time_range, [:epoch, :JulUTtime, :B])
    local bwave = mag_wave(mag_data, time_range)
    local bm = MVA(bwave)
    local min_hat = bm.vectors[:, 1]
    local mid_hat = bm.vectors[:, 2]
    local max_hat = bm.vectors[:, 3]
    local R = [max_hat mid_hat min_hat]
    local bWaveMVA = bwave * R

    local lims_max = ceil(maximum(abs.(bWaveMVA)))
    local lims_min = -lims_max
    local xtk = datetime2julian.(time_range)
    local t0 = mag_data[:JulUTtime][1]
    xtk = xtk .- t0
    limits!(ax1, xtk[1], xtk[end], lims_min, lims_max)
    limits!(ax2, xtk[1], xtk[end], lims_min, lims_max)
    limits!(ax3, xtk[1], xtk[end], lims_min, lims_max)
    # ylims!(ax1, lims_min, lims_max)
    # ylims!(ax2, lims_min, lims_max)
    # ylims!(ax3, lims_min, lims_max)
    lines!(ax1, mag_data[:JulUTtime].-t0, bWaveMVA[:, 1], 
        label = L"$B_{\mathrm{max}}$", color=color, linewidth=lw)
    lines!(ax2, mag_data[:JulUTtime].-t0, bWaveMVA[:, 2], 
        label = L"$B_{\mathrm{mid}}$", color=color, linewidth=lw)
    lines!(ax3, mag_data[:JulUTtime].-t0, bWaveMVA[:, 3], 
        label = L"$B_{\mathrm{min}}$", color=color, linewidth=lw)
    ax1.xticks = (xtk, Dates.format.(time_range, "HH:MM:SS"))
    ax2.xticks = (xtk, Dates.format.(time_range, "HH:MM:SS"))
    ax3.xticks = (xtk, Dates.format.(time_range, "HH:MM:SS"))
end

end