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


const Rm = 3390.0 #km
#Edberg, N. J. T., M. Lester, S. W. H. Cowley, and A. I. Eriksson (2008), Statistical analysis of the location of the Martian magnetic pileup boundary and bow shock and the influence of crustal magnetic fields, J. Geophys. Res., 113, A08206, doi:10.1029/2008JA013096.
#bow-shock model
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
#magnetopause model
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
datapath2c1hz = "E:/work/Tianwen/Data/MAVEN"
date = DateTime(2014, 10, 16, 00, 00, 00) .+ Dates.Day(1) .* range(0, 3650)
#date = DateTime(2014, 10, 11, 00, 00, 00) .+ Dates.Day(1) .* range(0, 1)
#20220630
datestr = Dates.format.(date, "yyyymmdd")

for (dts, dte) in zip(datestr, date)
	year = Dates.format(dte, "yyyy")
    file2c1hz = glob("*$dts*.sts", datapath2c1hz)
    println(file2c1hz)
	if isempty(file2c1hz) 
		continue
	end
	file2c1hz=file2c1hz[1]
	skip = 1 
	for line in eachline(file2c1hz)
		if line[1:6] != "  $year"
			skip = skip+1
		else
			break
		end
	end
    global mag2c1hz = identity.(DataFrame(readdlm(file2c1hz, skipstart=skip), :auto))

	global BMSO1hz = Array(mag2c1hz[:, 8:10])
	global  magut1hz = @. DateTime(mag2c1hz[:, 1]) + Day(mag2c1hz[:, 2]-1) + Hour(mag2c1hz[:, 3]) + Minute(mag2c1hz[:, 4]) + Second(mag2c1hz[:, 5]) + Millisecond(mag2c1hz[:, 6])
	global  magjlut1hz = datetime2julian.(magut1hz)
	magjlut1hz= magjlut1hz.-magjlut1hz[1]
	global  PosMSO1hz = Array(mag2c1hz[:, 12:14])./Rm


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
        a1 = arc!(ax1, Point2f(0), 1, 0, π; color=:white, lw=2)
        phi = range(0, 0.5π; length=180)
        x = [cos.(phi); 0.0]
        y = [sin.(phi); 0.0]
        pn = Polygon(Point2f[(xi, yi) for (xi, yi) in zip(x, y)])
        poly!(ax1, pn, color = :white)

        l1 = lines!(ax1, xshock, ryzshock; color=:white, lw=2, overdraw = true)
        l2 = lines!(ax1, xmp, ryzmp; color=:white, lw=2, overdraw = true)
        l3 = lines!(ax1, PosMSO1hz[:, 1], sqrt.(PosMSO1hz[:, 2].^2+PosMSO1hz[:, 3].^2); color=:white, overdraw = true)
        #Legend(fig[1,4], ax1, padding = (50, 0, 0, 0))
#        axislegend(ax1, framevisible = false)


        ax2 = Axis(fig[1,2], aspect=DataAspect(), xlabel=L"$x$ ($R_{\mathrm{M}}$)", ylabel=L"$y$ ($R_{\mathrm{M}}$)", xticks=range(-10, 10, 11), yticks=range(-10, 10, 11))
        xlims!(ax2, -10, 10) 
        ylims!(ax2, -5, 5) 
        a1 = arc!(ax2, Point2f(0), 1, 0, 2π; color=:white, lw=2)
        phi = range(-0.5π, 0.5π; length=180)
        x = [cos.(phi); 0.0]
        y = [sin.(phi); 0.0]
        pn = Polygon(Point2f[(xi, yi) for (xi, yi) in zip(x, y)])
        poly!(ax2, pn, color = :white)

        #l1 = lines!(ax2, xshock, ryzshock; color=:white, lw=2)
        #l2 = lines!(ax2, xmp, ryzmp; color=:white, lw=2)
        l3 = lines!(ax2, PosMSO1hz[:, 1], PosMSO1hz[:, 2]; color=:white, overdraw = true)
        ax2.title = dts

        ax3 = Axis(fig[1,3], aspect=DataAspect(), xlabel=L"$x$ ($R_{\mathrm{M}}$)", ylabel=L"$z$ ($R_{\mathrm{M}}$)", xticks=range(-10, 10, 11), yticks=range(-10, 10, 11))
        xlims!(ax3, -10, 10) 
        ylims!(ax3, -5, 5) 
        a1 = arc!(ax3, Point2f(0), 1, 0, 2π; color=:white, lw=2)
        phi = range(-0.5π, 0.5π; length=180)
        x = [cos.(phi); 0.0]
        y = [sin.(phi); 0.0]
        pn = Polygon(Point2f[(xi, yi) for (xi, yi) in zip(x, y)])
        poly!(ax3, pn, color = :white)

        #l1 = lines!(ax3, xshock, ryzshock; color=:white, lw=2)
        #l2 = lines!(ax3, xmp, ryzmp; color=:white, lw=2)
        l3 = lines!(ax3, PosMSO1hz[:, 1], PosMSO1hz[:, 3]; color=:white, overdraw = true)
        #=
        hidexdecorations!(ax1, grid=false)
        hidexdecorations!(ax2, grid=false)
        linkxaxes!(ax1, ax2, ax3)
        =#

        for ti in 0:8
            local ind = findall(x-> x<dte + Dates.Hour((hr-1)*4) + Dates.Minute(ti*30) + Dates.Second(1)  && x>dte + Dates.Hour((hr-1)*4) + Dates.Minute(ti*30) - Dates.Second(1), magut1hz) 
            local str = Dates.format(dte + Dates.Hour((hr-1)*4) + Dates.Minute(ti*30), "HH:MM")
            scatter!(ax1, PosMSO1hz[ind, 1], sqrt.(PosMSO1hz[ind, 2].^2+PosMSO1hz[ind, 3].^2), markersize=20, label=str, color=ti+1, colorrange=(1, 9), colormap=:tab10)
            text!(ax1, 0.9, 0.9-ti*0.1, text = str, font = :bold, align = (:center, :center), space = :relative, fontsize = 25, color=ti+1, colorrange=(1, 9), colormap=:tab10)

            scatter!(ax2, PosMSO1hz[ind, 1], PosMSO1hz[ind, 2], markersize=20, label=str, color=ti+1, colorrange=(1, 9), colormap=:tab10)
            text!(ax2, 0.9, 0.9-ti*0.1, text = str, font = :bold, align = (:center, :center), space = :relative, fontsize = 25, color=ti+1, colorrange=(1, 9), colormap=:tab10)

            scatter!(ax3, PosMSO1hz[ind, 1], PosMSO1hz[ind, 3], markersize=20, label=str, color=ti+1, colorrange=(1, 9), colormap=:tab10)
            text!(ax3, 0.9, 0.9-ti*0.1, text = str, font = :bold, align = (:center, :center), space = :relative, fontsize = 25, color=ti+1, colorrange=(1, 9), colormap=:tab10)
        end


        xd = dte .+ Dates.Hour((hr-1)*4) .+ Dates.Minute(10) .* range(0, 12)
        #xd = DateTime(2022, 03, 23, 22, 50, 0) .+ Dates.Minute(1) .* range(0, 10)

        ind = findall(x-> x<=maximum(xd)  && x>=minimum(xd), magut1hz) 
        xtk = datetime2julian.(xd)
        xtk = Float32.(xtk .- datetime2julian(magut1hz[1]))

        ax4 = Axis(fig[2, 1:3], xlabel = "UT", ylabel = L"Magnetic $\mathbf{B}$ (nT)", xticks=xtk)
        ax4.xticks = (xtk, Dates.format.(xd, "HH:MM"))
        xlims!(ax4, minimum(xtk), maximum(xtk))
		#=
        if !isnan(minimum(Matrix(BMSO1hz[ind, :]))) && !isnan(maximum(Matrix(BMSO1hz[ind, :])))
            ylims!(ax4, minimum(Matrix(BMSO1hz[ind, :]))-0.1, maximum(Matrix(BMSO1hz[ind, :]))+0.1)
        end
		=#
        #=
        println(minimum(minimum(eachcol(BMSO1hz[ind, :]))), "**", maximum(maximum(eachcol(BMSO1hz[ind, :]))))
        println(minimum(xd), maximum(xd))
        =#
		if !isempty(ind) 
			lines!(ax4, magjlut1hz[ind], BMSO1hz[ind, 1], label = L"$B_{\mathrm{x}}$", overdraw = true, linewidth=0.8)
			lines!(ax4, magjlut1hz[ind], BMSO1hz[ind, 2], label = L"$B_{\mathrm{y}}$", overdraw = true, linewidth=0.8)
			lines!(ax4, magjlut1hz[ind], BMSO1hz[ind, 3], label = L"$B_{\mathrm{z}}$", overdraw = true, linewidth=0.8)
            Legend(fig[2,4], ax4)
		end
        #lines!(ax4, magjlut1hz,  sqrt.(BMSO1hz[:, 1].^2+BMSO1hz[:, 2].^2+BMSO1hz[:, 3].^2), label = L"$B_{\mathrm{t}}$", overdraw = true, color=:white)

        #=
        lines!(ax4, magjlut1hz, BMSO1hz[:,1], label = L"$B_{\mathrm{x}}$", overdraw = true, colormap=:darkrainbow)
        lines!(ax4, magjlut1hz, BMSO1hz[:,2], label = L"$B_{\mathrm{y}}$", overdraw = true, colormap=:darkrainbow)
        lines!(ax4, magjlut1hz, BMSO1hz[:,3], label = L"$B_{\mathrm{z}}$", overdraw = true, colormap=:darkrainbow)
        =#
        


        ax5=Axis(fig[3,1:3], xlabel = "UT", ylabel = L"Period $T$ (s)",xticks=xtk, yscale=log2)
        ax5.xticks = (xtk, Dates.format.(xd, "HH:MM"))
        xlims!(ax5, minimum(xtk), maximum(xtk)) 
        ylims!(ax5, 2, 128)
        ax5.yreversed = true
		if !isempty(ind) 
			hp1 = heatmap!(ax5, magjlut1hz, period1hz, Bpower1hz', colorscale=log10, colormap=:gist_earth, colorrange=(2e-2, 3e2))
			cb = Colorbar(fig[3, 4], hp1, label = L"Wavelet Power $P_{\mathrm{B}}$")
		end
        tightlimits!(ax5)
        c4 = !isempty(ind) 

        #


        xd = dte .+ Dates.Hour((hr-1)*4+2) .+ Dates.Minute(10) .* range(0, 12)
        ind = findall(x-> x<=maximum(xd)  && x>=minimum(xd), magut1hz) 
        #xd = DateTime(2022, 03, 23, 22, 50, 0) .+ Dates.Minute(1) .* range(0, 10)
        xtk = datetime2julian.(xd)
        xtk = Float32.(xtk .- datetime2julian(magut1hz[1]))

        ax6 = Axis(fig[4, 1:3], xlabel = "UT", ylabel = L"Magnetic $\mathbf{B}$ (nT)", xticks=xtk)
        ax6.xticks = (xtk, Dates.format.(xd, "HH:MM"))
        xlims!(ax6, minimum(xtk), maximum(xtk)) 
		#=
        if !isnan(minimum(Matrix(BMSO1hz[ind, :]))) && !isnan(maximum(Matrix(BMSO1hz[ind, :])))
            ylims!(ax6, minimum(Matrix(BMSO1hz[ind, :]))-0.1, maximum(Matrix(BMSO1hz[ind, :]))+0.1)
        end
		=#
        #=
        println(minimum(minimum(eachcol(BMSO1hz[ind, :]))), "**", maximum(maximum(eachcol(BMSO1hz[ind, :]))))
        println(minimum(xd), maximum(xd))
        =#
        if !isempty(ind) 
			lines!(ax6, magjlut1hz[ind], BMSO1hz[ind, 1], label = L"$B_{\mathrm{x}}$", overdraw = true, linewidth=0.8)
			lines!(ax6, magjlut1hz[ind], BMSO1hz[ind, 2], label = L"$B_{\mathrm{y}}$", overdraw = true, linewidth=0.8)
			lines!(ax6, magjlut1hz[ind], BMSO1hz[ind, 3], label = L"$B_{\mathrm{z}}$", overdraw = true, linewidth=0.8)
		end
        #lines!(ax6, magjlut1hz,  sqrt.(BMSO1hz[:, 1].^2+BMSO1hz[:, 2].^2+BMSO1hz[:, 3].^2), label = L"$B_{\mathrm{y}}$", overdraw = true, color=:white)
        #=
        lines!(ax6, magjlut1hz, BMSO1hz[:,1], label = L"$B_{\mathrm{x}}$", overdraw = true, colormap=:darkrainbow)
        lines!(ax6, magjlut1hz, BMSO1hz[:,2], label = L"$B_{\mathrm{y}}$", overdraw = true, colormap=:darkrainbow)
        lines!(ax6, magjlut1hz, BMSO1hz[:,3], label = L"$B_{\mathrm{z}}$", overdraw = true, colormap=:darkrainbow)
        Legend(fig[4,4], ax6)
        =#

        ax7 = Axis(fig[5,1:3], xlabel = "UT", ylabel = L"Period $T$ (s)", xticks=xtk, yscale=log2)
        ax7.xticks = (xtk, Dates.format.(xd, "HH:MM"))
        xlims!(ax7, minimum(xtk), maximum(xtk)) 
        ylims!(ax7, 2, 128)
        ax7.yreversed = true
		if !isempty(ind)
			hp1 = heatmap!(ax7, magjlut1hz, period1hz, Bpower1hz', colorscale=log10, colormap=:gist_earth, colorrange=(2e-2, 3e2))
			cb = Colorbar(fig[5, 4], hp1, label = L"Wavelet Power $P_{\mathrm{B}}$")
		end
        #cr1 = contourf!(ax7, magut1hz, period1hz, log10.(Bcomppower1hz'), levels=range(-2, 2, 12), colormap=:viridis, extendlow = :auto, extendhigh = :auto)
        tightlimits!(ax7)

        c4 = c4 || !isempty(ind)
        

        colgap!(fig.layout, 10)
        rowgap!(fig.layout, 12)
        rowsize!(fig.layout, 1, Relative(1.5/5))
        if c4 
            colsize!(fig.layout, 4, Relative(1/25))
        end
        #colsize!(fig.layout, 1, Relative(1/5))

        hidexdecorations!(ax4, grid=false)
        hidexdecorations!(ax6, grid=false)
        linkxaxes!(ax4, ax5)
        linkxaxes!(ax6, ax7)
        resize_to_layout!(fig)
        #save("overview-1hz.pdf", fig, pt_per_unit = 1)
        save("E:/MAVEN/MAG/MagPower1hz/" * dts * string(hr) * ".png", fig, pt_per_unit = 1)
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





