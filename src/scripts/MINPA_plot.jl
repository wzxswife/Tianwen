module MINPA_plot
using Dates, Statistics
using GLMakie, LinearAlgebra, CoordinateTransformations

include("MAVEN_SWIA.jl")
include("SphericalSplines.jl")
using .MAVEN_SWIA
using .SphericalSplines

export plot_vdf_slices!, plot_vdf_slice_vperp1!, plot_arrow!
export TW_MINPA_eflux_Overview_plot

global colorrange = (1e4, 1e8)
global cmap = :gist_earth


""" 绘制 v∥ 切片画廊 (含透明度支持) """
function plot_vdf_slices!(scene, result; 
    colorrange=colorrange, colormap=cmap, alpha=1.0)
    vpara  = result.vpara
    vperp1 = result.vperp1
    vperp2 = result.vperp2
    rect   = result.rect

    for ipara in 1:length(vpara)
        vparai = vpara[ipara]
        heatmap!(scene, vperp1, vperp2, rect[:, ipara, :],
                 transformation = (:xz, vparai),
                 colormap = colormap, colorrange = colorrange,
                 colorscale = log10, alpha = alpha)
    end
end

""" 绘制 v⊥1≈0 剖面 (含透明度支持) """
function plot_vdf_slice_vperp1!(scene, result; 
    colorrange=colorrange, colormap=cmap, alpha=1.0)
    vpara  = result.vpara
    vperp1 = result.vperp1
    vperp2 = result.vperp2
    rect   = result.rect

    iperp1 = findmin(abs.(vperp1))[2]
    vperp1i = vperp1[iperp1]
    heatmap!(scene, vpara, vperp2, rect[iperp1, :, :],
             transformation = (:yz, vperp1i),
             colormap = colormap, colorrange = colorrange,
             colorscale = log10, alpha = alpha)
end

""" 添加箭头 (标注峰值方向) """
function plot_arrow!(scene, result; color=:red)
    efxmean = result.efxmean
    paphi   = result.paphi
    patheta = result.patheta
    pvel    = result.pvelocity
    perp1   = result.perpminpa1
    para    = result.paraminpa
    perp2   = result.perpminpa2

    i, j = 1, 13
    ind = argmax(efxmean[1, j, i, :])
    v1 = CartesianFromSpherical()(CoordinateTransformations.Spherical(
        pvel[1, 22], deg2rad(paphi[j]), deg2rad(patheta[i]+22.5)))
    ps = Point3f(0.0, 0.0, 0.0)
    pe = Point3f(dot(v1, perp1), dot(v1, para), dot(v1, perp2))
    arrows3d!(scene, ps, pe, argmode=:endpoint, color=color,
              shaftradius=0.02, tipradius=0.05, tailradius=0.02)
end

function TW_MINPA_eflux_Overview_plot(minpa_data::Dict, time_range::Vector{DateTime}, im::Int)

    size_inches = (20, 30)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 20, figure_padding = (5,20, 9, 10))
    cn = cgrad(:seaborn_dark, 10, categorical = true)
    lcn = Makie.wong_colors()

    xd = time_range
    xtk = datetime2julian.(xd)
    xtk = Float32.(xtk .- minpa_data[:JulUTtime][1])
    idx = findall(time_range[begin] .<= minpa_data[:epoch] .<= time_range[end])
    axes = [Axis(fig[j, i], yscale=log10, xticks = (xtk, Dates.format.(xd, "HH:MM")), backgroundcolor=:gray80, xticklabelrotation=-π/3) for j in 1:16, i in 1:4]

    for j in 1:16
        for i in 1:4
            ax = axes[j, i]
            heatmap!(ax, Float32.(minpa_data[:JulUTtime].-minpa_data[:JulUTtime][1]), 
                minpa_data[:energy], minpa_data[:eflux][:,im,j,i,:], 
                colormap=:gist_earth, colorscale=log10, colorrange = (1e5, 1e9))
            text!(ax, 0.0, 1.0, text = rich("θ=$i ϕ=$j"), 
                align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 20, color = :white)
            xlims!(ax, xtk[begin], xtk[end]) 
            ylims!(ax, 1e1, 1e4)

            # 在 time_range 内沿时间平均，找出平均通量最大的能道，用红色横线标注能量值（keV）
            if !isempty(idx)
                eflux_slice = minpa_data[:eflux][idx, im, j, i, :]  # [nt, nenergy]
                mean_flux = vec(mean(eflux_slice, dims=1))           # 每个能道的平均通量
                max_flux, max_energy_idx = findmax(mean_flux)
                if !isnan(max_flux) && isfinite(max_flux)
                    max_energy_eV = minpa_data[:energy][max_energy_idx]
                    hlines!(ax, [max_energy_eV], color=:red, linewidth=2)
                    text!(ax, xtk[end], max_energy_eV, text=rich("Max: ", rich("$(round(max_energy_eV / 1000, digits=2)) keV", color=:red, fontweight=:bold)),
                        align=(:right, :bottom), fontsize=14, strokewidth=3, strokecolor=:black, color=:white)
                end
            end
            
            if im == 1 && i==1 && j==12
                # annotation!(ax, -80, -40, datetime2julian(DateTime(2022, 06, 12, 10, 00, 00))-minpajulianmod1[1], 600, path = Ann.Paths.Corner(), text = "SW", color=cn[2], fontsize=20)
                # annotation!(ax, -120, -60, datetime2julian(DateTime(2022, 06, 12, 10, 00, 00))-minpajulianmod1[1], 2000, path = Ann.Paths.Corner(), text = "PU", color=cn[2], fontsize=20)
                annotation!(ax, -80, -40, datetime2julian(xd[1])-minpa_data[:JulUTtime][1], 600, path = Ann.Paths.Corner(), text = "SW", color=cn[2], fontsize=20)
                annotation!(ax, -120, -60, datetime2julian(xd[1])-minpa_data[:JulUTtime][1], 2000, path = Ann.Paths.Corner(), text = "PU", color=cn[2], fontsize=20)
            end           
        end
    end
    linkaxes!(axes...)
    hidexdecorations!.(axes[1:15,:], grid=false)
    hideydecorations!.(axes[:, 2:end], grid=false)
    colgap!(fig.layout, 30)

    sp = [rich(rich("H"), superscript("+")), rich(rich("He"), superscript("++")), rich("m/q=4"), rich(rich("O"), superscript("+")), rich("m/q=28"), rich(rich("O"), subscript("2"), superscript("+")), rich(rich("CO"), subscript("2"), superscript("+")), rich("m/q=64")]

    Colorbar(fig[1:16, 5], colormap=:gist_earth, scale=log10, colorrange = (1e5, 1e9), label = rich("Tianwen-1 MINPA ") * sp[im] * rich(rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")))
        Label(fig[1:16, 0], rich("Energy (eV)"), fontsize = 25, halign = :center, rotation = π/2)
        Label(fig[17, 1:5], rich("Universal Time"), fontsize = 20, halign = :center)

    return fig
end
# TODO:完善MINPA pitch angle函数
function TW_MINPA_pitch_angle(ax, minpa_data::Dict, mag_data::Dict, im::Int,
    time_range::Vector{DateTime}, energy_range::Vector{Float64}; 
    c_range=(1e3, 1e8), colormap=cgrad(:RdYlBu, rev=true), colorscale=log10, 
    threshold=1e0)
    
    local time_idx_mag = findall(time_range[begin].<mag_data[:epoch].<time_range[end])
    local mag_epochs = mag_data[:epoch][time_idx_mag]
    local energies = minpa_data[:energy][:] ./ 1000 
    local time_idx = findall(time_range[begin].<minpa_data[:epoch].<time_range[end])
    local energy_idx = findall(energy_range[begin].<energies.<energy_range[end]) 
    local ntime = length(time_idx)
    local nenergy = length(energy_idx)
    local ntheta = 4
    local nphi = 16
    local minpa_epochs = minpa_data[:epoch][time_idx]
    local time_data = datetime2julian.(minpa_epochs) 
    local flux_4d = minpa_data[:eflux][time_idx, im, :, :, energy_idx]
    local pitch_angle_4d = Array{Float64, 4}(undef, ntime, nphi, ntheta, nenergy)
    for i in 1:ntime
        local energy = reshape(minpa_data[:energy][energy_idx],1,1,nenergy)
        local phi = reshape(minpa_data[:phi],nphi,1,1)
        local theta = reshape(minpa_data[:theta],1,ntheta,1)
        local v0 = MAVEN_SWIA.ion_energy2v.(energy,1)
        local vv = MAVEN_SWIA.sphere2xyz_for_SWIA.(v0,theta,phi)
        local mag_idx = argmin(abs.(mag_epochs .- minpa_epochs[i]))
        local B = convert(Vector{Float64},mag_data[:B][mag_idx, :])
        pitch_angle_4d[i, :, :, :] = angle_between.(vv, Ref(B))
    end

    local pitch_angle_2d = reshape(pitch_angle_4d, ntime, nphi*ntheta*nenergy)
    local flux_2d = reshape(flux_4d, ntime, nphi*ntheta*nenergy)

    # 创建 pitch angle bins (0-180度)
    local n_bins = 60
    local bin_edges = range(0, 180, length=n_bins+1)
    local bin_centers = (bin_edges[1:end-1] + bin_edges[2:end]) / 2
    # 构建热力图矩阵
    local heatmap_data = zeros(Float64, ntime, n_bins)
    for t in 1:ntime
        pitch_angle_t = pitch_angle_2d[t, :]
        flux_t = flux_2d[t, :]
        
        for i in 1:(nphi*ntheta*nenergy)
            pa_val = pitch_angle_t[i]
            flux_val = flux_t[i]
            isnan(pa_val) && continue 
            # 找到对应的 bin
            bin_idx = floor(Int, pa_val / 180 * n_bins) + 1
            bin_idx = clamp(bin_idx, 1, n_bins)
            
            heatmap_data[t, bin_idx] += flux_val
        end
    end

    local xtk = datetime2julian.(time_range)

    # 将小于阈值的值替换为 NaN
    local heatmap_data_masked = ifelse.(heatmap_data .< threshold, NaN, heatmap_data)
    hm = heatmap!(ax, time_data, range(0, 180, length=n_bins), heatmap_data_masked, 
        colormap=colormap, colorscale=colorscale, colorrange=c_range)
    xlims!(ax, xtk[1], xtk[end])
    ax.xticks = (xtk, Dates.format.(time_range, "HH:MM:SS"))
    ax.xminorticksvisible=true
    # ax.xminorgridvisible=true
    ax.xminorticks=IntervalsBetween(10)
    ax.yticks = 0:30:180
    return hm
end

# 常用函数
function angle_between(v1, v2)
    cosθ = dot(v1, v2) / (norm(v1) * norm(v2))
    return acosd(clamp(cosθ, -1.0, 1.0))  # 限制范围防止数值误差
end

end