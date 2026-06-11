module MINPA_plot
using Dates
using GLMakie, LinearAlgebra, CoordinateTransformations

include("SphericalSplines.jl")
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

function TW_MINPA_eflux_Overview_plot(minpa_data::Dict, 
    time_range::Vector{DateTime}, im::Int)
    size_inches = (20, 30)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 20, figure_padding = (5,20, 9, 10))
    cn = cgrad(:seaborn_dark, 10, categorical = true)
    lcn = Makie.wong_colors()

    xd = time_range
    xtk = datetime2julian.(xd)
    xtk = Float32.(xtk .- minpa_data[:JulUTtime][1])

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

function TW_MINPA_pitch_angle(ax, minpa_data::Dict, mag_data::Dict, im::Int,
    time_range::Vector{DateTime}, energy_range::Vector{Float64}; 
    c_range=(1e3, 1e8), colormap=cgrad(:RdYlBu, rev=true), colorscale=log10)
    
end

end