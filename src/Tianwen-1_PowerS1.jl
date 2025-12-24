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
include("TCWavelet.jl")
const Rm = 3390.0e3  #m
#Edberg, N. J. T., M. Lester, S. W. H. Cowley, and A. I. Eriksson (2008), Statistical analysis of the location of the Martian magnetic pileup boundary and bow shock and the influence of crustal magnetic fields, J. Geophys. Res., 113, A08206, doi:10.1029/2008JA013096.
"bow-shock model"
function bowshock(xshock)
    xF = 0.55 # Rm
    ϵ = 1.05
    L = 2.10 # rm
    rSD = 1.58
    temp = (ϵ^2-1.0)*(xshock-xF)^2-2ϵ*L*(xshock-xF)+L^2
    if temp>=0 
        return sqrt(temp)
    else
        return Inf64
    end
end
"magnetopause model"
function magnetopause(xmp)
    rSD = 1.33
	xF = 0.86
	ϵ = 0.92
	L = 0.90
    temp = (ϵ^2-1.0)*(xmp-xF)^2-2ϵ*L*(xmp-xF)+L^2
    if temp>=0 
        return sqrt(temp)
    else
        return Inf64
    end
end

np = 5001
xshock = range(-10, 3, length=np)
ryzshock = bowshock.(xshock)
ind = findall(isfinite, ryzshock)
xshock = [xshock[ind]; reverse(xshock[ind])]
ryzshock = [ryzshock[ind]; -reverse(ryzshock[ind])]
xmp = range(-10, 3, length=np)
ryzmp = magnetopause.(xmp)
ind = findall(isfinite, ryzmp)
xmp = [xmp[ind]; reverse(xmp[ind])]
ryzmp = [ryzmp[ind]; -reverse(ryzmp[ind])]


# read data 
project_root = dirname(@__DIR__)
datapath2c1hz = joinpath(project_root, "Data", "1Hz")
# date = DateTime(2021, 11, 16, 00, 00, 00) .+ Dates.Day(1) .* range(0, 226)
date = DateTime(2021, 12, 01, 00, 00, 00) .+ Dates.Day(1) .* range(0, 1)
#date = DateTime(2021, 11, 16, 00, 00, 00) .+ Dates.Day(1) .* range(0, 0)
#20220630
datestr = Dates.format.(date, "yyyymmdd")

for (dts, dte) in zip(datestr, date)
    filename = "TW1_MOMAG_MSO_01Hz_" * dts * "_2C_v03.dat"
    file2c1hz = joinpath(datapath2c1hz, filename)
    println(file2c1hz)
    global mag2c1hz = identity.(DataFrame(readdlm(file2c1hz, skipstart=19), :auto))
    name = ["Time", "Sampling_Rate", "X_MSO", "Y_MSO", "Z_MSO", "Probe_Position_X_MSO", "Probe_Position_Y_MSO", "Probe_Position_Z_MSO", "Roll", "Pitch", "Yaw",  "Quality_Flags"]
    rename!(mag2c1hz, name)
    mag2c1hz[!, :Time] = map(x->DateTime(x[begin:end-4], DateFormat("y-m-dTH:M:S.s")), mag2c1hz[!, :Time])
    mag2c1hz[!, :JulUT] = datetime2julian.(mag2c1hz[!, :Time])
    unique!(mag2c1hz) #remove depulicate rows
    sort!(mag2c1hz) #sorting
    magut1hz = mag2c1hz[:, :JulUT].-mag2c1hz[1, :JulUT]
    global BMSO1hz = mag2c1hz[:, [:X_MSO, :Y_MSO, :Z_MSO]]
    ind0 = findall(x-> !isnan(x), BMSO1hz[!, 1] ) 
    ind1 = findall(x-> isnan(x), BMSO1hz[!, 1] ) 
    if length(ind1) >= 1 && length(ind0) > 1
        for ib in 1:3 
            interp_linear = linear_interpolation(magut1hz[ind0], BMSO1hz[ind0, ib]; extrapolation_bc=Line())
            BMSO1hz[ind1, ib] = interp_linear(magut1hz[ind1])
        end
    end
    
    #wavelet power
    ns1hz = length(magut1hz)
    dt1hz = 1.0

    mother = "MORLET"
    wave1hz, period1hz, scale1hz, coi1hz = wavelet(reshape(BMSO1hz[:, 1], ns1hz), dt1hz; pad=1, mother=mother)
    xpower1hz = abs.(wave1hz).^2
    wave1hz, period1hz, scale1hz, coi1hz = wavelet(reshape(BMSO1hz[:, 2], ns1hz), dt1hz; pad=1, mother=mother)
    ypower1hz = abs.(wave1hz).^2
    wave1hz, period1hz, scale1hz, coi1hz = wavelet(reshape(BMSO1hz[:, 3], ns1hz), dt1hz; pad=1, mother=mother)
    zpower1hz = abs.(wave1hz).^2
    Bpower1hz =xpower1hz+ypower1hz+zpower1hz

    
    #overview
    using CairoMakie, GeometryBasics
    CairoMakie.activate!()
    dark_latexfonts = merge(theme_dark(), theme_latexfonts())
    set_theme!(dark_latexfonts)
    #update_theme!(fontsize=300)  #, figure_padding = (1, 0, 1, 5)
    size_inches = (30, 17)
    size_pt = 72 .* size_inches
    #fig = Figure(resolution = size_pt, fontsize = 12)
    for hr in 1:6 
        fig = Figure(size = size_pt, fontsize = 25)
        #g1 = fig[1, 1] = GridLayout()
        #g2 = fig[2, 1] = GridLayout()

        ax1 = Axis(fig[1,1], aspect=DataAspect(), xlabel=L"$x$ ($R_{\mathrm{M}}$)", ylabel=L"$\sqrt{y^2+z^2}$ ($R_{\mathrm{M}}$)", xticks=range(-10, 10), yticks=range(-10, 10))
        xlims!(ax1, -5, 5) 
        ylims!(ax1, 0, 5) 
        a1 = arc!(ax1, Point2f(0), 1, 0, π; color=:white, linewidth=2)
        phi = range(0, 0.5π; length=180)
        x = [cos.(phi); 0.0]
        y = [sin.(phi); 0.0]
        pn = Polygon(Point2f[(xi, yi) for (xi, yi) in zip(x, y)])
        poly!(ax1, pn, color = :white)

        l1 = lines!(ax1, xshock, ryzshock; color=:white, linewidth=2, overdraw = true)
        l2 = lines!(ax1, xmp, ryzmp; color=:white, linewidth=2, overdraw = true)
        l3 = lines!(ax1, mag2c1hz[!, :Probe_Position_X_MSO], sqrt.(mag2c1hz[!, :Probe_Position_Y_MSO].^2+mag2c1hz[!, :Probe_Position_Z_MSO].^2); color=:white, overdraw = true)
        #Legend(fig[1,4], ax1, padding = (50, 0, 0, 0))
#        axislegend(ax1, framevisible = false)


        ax2 = Axis(fig[1,2], aspect=DataAspect(), xlabel=L"$x$ ($R_{\mathrm{M}}$)", ylabel=L"$y$ ($R_{\mathrm{M}}$)", xticks=range(-10, 10, 11), yticks=range(-10, 10, 11))
        xlims!(ax2, -10, 10) 
        ylims!(ax2, -5, 5) 
        a1 = arc!(ax2, Point2f(0), 1, 0, 2π; color=:white, linewidth=2)
        phi = range(-0.5π, 0.5π; length=180)
        x = [cos.(phi); 0.0]
        y = [sin.(phi); 0.0]
        pn = Polygon(Point2f[(xi, yi) for (xi, yi) in zip(x, y)])
        poly!(ax2, pn, color = :white)

        #l1 = lines!(ax2, xshock, ryzshock; color=:white, linewidth=2)
        #l2 = lines!(ax2, xmp, ryzmp; color=:white, linewidth=2)
        l3 = lines!(ax2, mag2c1hz[!, :Probe_Position_X_MSO], mag2c1hz[!, :Probe_Position_Y_MSO]; color=:white, overdraw = true)
        ax2.title = dts

        ax3 = Axis(fig[1,3], aspect=DataAspect(), xlabel=L"$x$ ($R_{\mathrm{M}}$)", ylabel=L"$z$ ($R_{\mathrm{M}}$)", xticks=range(-10, 10, 11), yticks=range(-10, 10, 11))
        xlims!(ax3, -10, 10) 
        ylims!(ax3, -5, 5) 
        a1 = arc!(ax3, Point2f(0), 1, 0, 2π; color=:white, linewidth=2)
        phi = range(-0.5π, 0.5π; length=180)
        x = [cos.(phi); 0.0]
        y = [sin.(phi); 0.0]
        pn = Polygon(Point2f[(xi, yi) for (xi, yi) in zip(x, y)])
        poly!(ax3, pn, color = :white)

        #l1 = lines!(ax3, xshock, ryzshock; color=:white, linewidth=2)
        #l2 = lines!(ax3, xmp, ryzmp; color=:white, linewidth=2)
        l3 = lines!(ax3, mag2c1hz[!, :Probe_Position_X_MSO], mag2c1hz[!, :Probe_Position_Z_MSO]; color=:white, overdraw = true)
        #=
        hidexdecorations!(ax1, grid=false)
        hidexdecorations!(ax2, grid=false)
        linkxaxes!(ax1, ax2, ax3)
        =#

        for ti in 0:8
            local ind = findall(x-> x<dte + Dates.Hour((hr-1)*4) + Dates.Minute(ti*30) + Dates.Second(1)  && x>dte + Dates.Hour((hr-1)*4) + Dates.Minute(ti*30) - Dates.Second(1), mag2c1hz[!, :Time]) 
            local str = Dates.format(dte + Dates.Hour((hr-1)*4) + Dates.Minute(ti*30), "HH:MM")
            scatter!(ax1, mag2c1hz[ind, :Probe_Position_X_MSO], sqrt.(mag2c1hz[ind, :Probe_Position_Y_MSO].^2+mag2c1hz[ind, :Probe_Position_Z_MSO].^2), markersize=20, label=str, color=ti+1, colorrange=(1, 9), colormap=:tab10)
            text!(ax1, 0.9, 0.9-ti*0.1, text = str, font = :bold, align = (:center, :center), space = :relative, fontsize = 25, color=ti+1, colorrange=(1, 9), colormap=:tab10)

            scatter!(ax2, mag2c1hz[ind, :Probe_Position_X_MSO], mag2c1hz[ind, :Probe_Position_Y_MSO], markersize=20, label=str, color=ti+1, colorrange=(1, 9), colormap=:tab10)
            text!(ax2, 0.9, 0.9-ti*0.1, text = str, font = :bold, align = (:center, :center), space = :relative, fontsize = 25, color=ti+1, colorrange=(1, 9), colormap=:tab10)

            scatter!(ax3, mag2c1hz[ind, :Probe_Position_X_MSO], mag2c1hz[ind, :Probe_Position_Z_MSO], markersize=20, label=str, color=ti+1, colorrange=(1, 9), colormap=:tab10)
            text!(ax3, 0.9, 0.9-ti*0.1, text = str, font = :bold, align = (:center, :center), space = :relative, fontsize = 25, color=ti+1, colorrange=(1, 9), colormap=:tab10)
        end


        xd = dte .+ Dates.Hour((hr-1)*4) .+ Dates.Minute(10) .* range(0, 12)
        #xd = DateTime(2022, 03, 23, 22, 50, 0) .+ Dates.Minute(1) .* range(0, 10)

        ind = findall(x-> x<=maximum(xd)  && x>=minimum(xd), mag2c1hz[!, :Time]) 
        xtk = datetime2julian.(xd)
        xtk = Float32.(xtk .- mag2c1hz[1, :JulUT])

        ax4 = Axis(fig[2, 1:3], xlabel = "UT", ylabel = L"Magnetic $\mathbf{B}$ (nT)", xticks=xtk)
        ax4.xticks = (xtk, Dates.format.(xd, "HH:MM"))
        xlims!(ax4, minimum(xtk), maximum(xtk))
        if !isnan(minimum(Matrix(BMSO1hz[ind, :]))) && !isnan(maximum(Matrix(BMSO1hz[ind, :])))
            ylims!(ax4, minimum(Matrix(BMSO1hz[ind, :]))-0.1, maximum(Matrix(BMSO1hz[ind, :]))+0.1)
        end
        #=
        println(minimum(minimum(eachcol(BMSO1hz[ind, :]))), "**", maximum(maximum(eachcol(BMSO1hz[ind, :]))))
        println(minimum(xd), maximum(xd))
        =#
        lines!(ax4, mag2c1hz[!, :JulUT].-mag2c1hz[1, :JulUT], mag2c1hz[!, :X_MSO], label = L"$B_{\mathrm{x}}$", overdraw = true, linewidth=0.8)
        lines!(ax4, mag2c1hz[!, :JulUT].-mag2c1hz[1, :JulUT], mag2c1hz[!, :Y_MSO], label = L"$B_{\mathrm{y}}$", overdraw = true, linewidth=0.8)
        lines!(ax4, mag2c1hz[!, :JulUT].-mag2c1hz[1, :JulUT], mag2c1hz[!, :Z_MSO], label = L"$B_{\mathrm{z}}$", overdraw = true, linewidth=0.8)
        #lines!(ax4, mag2c1hz[!, :JulUT].-mag2c1hz[1, :JulUT],  sqrt.(mag2c1hz[!, :X_MSO].^2+mag2c1hz[!, :Y_MSO].^2+mag2c1hz[!, :Z_MSO].^2), label = L"$B_{\mathrm{t}}$", overdraw = true, color=:white)

        #=
        lines!(ax4, mag2c1hz[!, :JulUT].-mag2c1hz[1, :JulUT], BMSO1hz[:,1], label = L"$B_{\mathrm{x}}$", overdraw = true, colormap=:darkrainbow)
        lines!(ax4, mag2c1hz[!, :JulUT].-mag2c1hz[1, :JulUT], BMSO1hz[:,2], label = L"$B_{\mathrm{y}}$", overdraw = true, colormap=:darkrainbow)
        lines!(ax4, mag2c1hz[!, :JulUT].-mag2c1hz[1, :JulUT], BMSO1hz[:,3], label = L"$B_{\mathrm{z}}$", overdraw = true, colormap=:darkrainbow)
        =#
        Legend(fig[2,4], ax4)


        ax5=Axis(fig[3,1:3], xlabel = "UT", ylabel = L"Period $T$ (s)",xticks=xtk, yscale=log2)
        ax5.xticks = (xtk, Dates.format.(xd, "HH:MM"))
        xlims!(ax5, minimum(xtk), maximum(xtk)) 
        ylims!(ax5, 2, 128)
        ax5.yreversed = true
        hp1 = heatmap!(ax5, magut1hz, period1hz, Bpower1hz', colorscale=log10, colormap=:gist_earth, colorrange=(2e-2, 3e2))
        #cr1 = contourf!(ax5, magut1hz, period1hz, log10.(Bcomppower1hz'), levels=range(-2, 2, 12), colormap=:viridis, extendlow = :auto, extendhigh = :auto)

        #lines!(ax5, mag2c1hz[!, :JulUT].-mag2c1hz[1, :JulUT],  1836/28 ./ sqrt.(mag2c1hz[!, :X_MSO].^2+mag2c1hz[!, :Y_MSO].^2+mag2c1hz[!, :Z_MSO].^2), label = L"$B_{\mathrm{y}}$", overdraw = true)
        tightlimits!(ax5)

        cb = Colorbar(fig[3, 4], hp1, label = L"Wavelet Power $P_{\mathrm{B}}$")


        xd = dte .+ Dates.Hour((hr-1)*4+2) .+ Dates.Minute(10) .* range(0, 12)
        ind = findall(x-> x<=maximum(xd)  && x>=minimum(xd), mag2c1hz[!, :Time]) 
        #xd = DateTime(2022, 03, 23, 22, 50, 0) .+ Dates.Minute(1) .* range(0, 10)
        xtk = datetime2julian.(xd)
        xtk = Float32.(xtk .- mag2c1hz[1, :JulUT])

        ax6 = Axis(fig[4, 1:3], xlabel = "UT", ylabel = L"Magnetic $\mathbf{B}$ (nT)", xticks=xtk)
        ax6.xticks = (xtk, Dates.format.(xd, "HH:MM"))
        xlims!(ax6, minimum(xtk), maximum(xtk)) 
        if !isnan(minimum(Matrix(BMSO1hz[ind, :]))) && !isnan(maximum(Matrix(BMSO1hz[ind, :])))
            ylims!(ax6, minimum(Matrix(BMSO1hz[ind, :]))-0.1, maximum(Matrix(BMSO1hz[ind, :]))+0.1)
        end
        #=
        println(minimum(minimum(eachcol(BMSO1hz[ind, :]))), "**", maximum(maximum(eachcol(BMSO1hz[ind, :]))))
        println(minimum(xd), maximum(xd))
        =#
        lines!(ax6, mag2c1hz[!, :JulUT].-mag2c1hz[1, :JulUT], mag2c1hz[!, :X_MSO], label = L"$B_{\mathrm{x}}$", overdraw = true, colormap=:darkrainbow, linewidth=0.8)
        lines!(ax6, mag2c1hz[!, :JulUT].-mag2c1hz[1, :JulUT], mag2c1hz[!, :Y_MSO], label = L"$B_{\mathrm{y}}$", overdraw = true, colormap=:darkrainbow, linewidth=0.8)
        lines!(ax6, mag2c1hz[!, :JulUT].-mag2c1hz[1, :JulUT], mag2c1hz[!, :Z_MSO], label = L"$B_{\mathrm{z}}$", overdraw = true, colormap=:darkrainbow, linewidth=0.8)
        #lines!(ax6, mag2c1hz[!, :JulUT].-mag2c1hz[1, :JulUT],  sqrt.(mag2c1hz[!, :X_MSO].^2+mag2c1hz[!, :Y_MSO].^2+mag2c1hz[!, :Z_MSO].^2), label = L"$B_{\mathrm{y}}$", overdraw = true, color=:white)
        #=
        lines!(ax6, mag2c1hz[!, :JulUT].-mag2c1hz[1, :JulUT], BMSO1hz[:,1], label = L"$B_{\mathrm{x}}$", overdraw = true, colormap=:darkrainbow)
        lines!(ax6, mag2c1hz[!, :JulUT].-mag2c1hz[1, :JulUT], BMSO1hz[:,2], label = L"$B_{\mathrm{y}}$", overdraw = true, colormap=:darkrainbow)
        lines!(ax6, mag2c1hz[!, :JulUT].-mag2c1hz[1, :JulUT], BMSO1hz[:,3], label = L"$B_{\mathrm{z}}$", overdraw = true, colormap=:darkrainbow)
        Legend(fig[4,4], ax6)
        =#

        ax7 = Axis(fig[5,1:3], xlabel = "UT", ylabel = L"Period $T$ (s)", xticks=xtk, yscale=log2)
        ax7.xticks = (xtk, Dates.format.(xd, "HH:MM"))
        xlims!(ax7, minimum(xtk), maximum(xtk)) 
        ylims!(ax7, 2, 128)
        ax7.yreversed = true
        hp1 = heatmap!(ax7, magut1hz, period1hz, Bpower1hz', colorscale=log10, colormap=:gist_earth, colorrange=(2e-2, 3e2))
        #cr1 = contourf!(ax7, magut1hz, period1hz, log10.(Bcomppower1hz'), levels=range(-2, 2, 12), colormap=:viridis, extendlow = :auto, extendhigh = :auto)
        tightlimits!(ax7)

        cb = Colorbar(fig[5, 4], hp1, label = L"Wavelet Power $P_{\mathrm{B}}$")

        colgap!(fig.layout, 10)
        rowgap!(fig.layout, 12)
        rowsize!(fig.layout, 1, Relative(1.5/5))
        colsize!(fig.layout, 4, Relative(1/25))
        #colsize!(fig.layout, 1, Relative(1/5))

        hidexdecorations!(ax4, grid=false)
        hidexdecorations!(ax6, grid=false)
        linkxaxes!(ax4, ax5)
        linkxaxes!(ax6, ax7)
        resize_to_layout!(fig)
        #save("overview-1Hz.pdf", fig, pt_per_unit = 1)
        save("./Picture/1Hz/" * dts * string(hr) * ".png", fig, pt_per_unit = 1)
    end
end
#=
B0MSO1hz = copy(BMSO1hz)
for ib in 1:3
    for jc in 5:length(magut1hz)-5
        B0MSO1hz[jc, ib] = mean(BMSO1hz[jc-4:jc+5, ib]) 
    end
end
B0MSO1hz = Array(B0MSO1hz)
=#