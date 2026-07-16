using CairoMakie
CairoMakie.activate!()
#Mars Orbit
function Ephem_plot()
    size_inches = (9, 9)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5,10, 3, 5))
    cn = cgrad(:seaborn_dark, 10, categorical = true)
    ax1 = Axis(fig[1,1], aspect=DataAspect(), xlabel=rich(rich("x", font = :italic), " (AU)"), ylabel=rich(rich("y", font = :italic), " (AU)"), xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticksvisible = true, backgroundcolor=:gray80) # xticks=range(0, 256, step=32)

    xlims!(ax1, -2, 2) 
    ylims!(ax1, -2, 2) 
    #lines!(ax1, merpos[1,:], merpos[2,:], linewidth=2)
    #lines!(ax1, venpos[1,:], venpos[2,:], linewidth=2)
    #lines!(ax1, earpos[1,:], earpos[2,:], linewidth=2)
    lines!(ax1, marpos[1,:], marpos[2,:], linewidth=2, linestyle=:dash, color=cn[2])
    scatter!(ax1, 0, 0, color=:gold, markersize =40, label="Sun", marker=:star8)
    #scatter!(ax1, merpos[1, 1], merpos[2, 1], markersize =20, label="Mercury")
    #scatter!(ax1, venpos[1, 1], venpos[2, 1], markersize =20, label="Venus")
    #scatter!(ax1, earpos[1, 1], earpos[2, 1], markersize =20, label="Earth")

    scatter!(ax1, marpos[1, pind], marpos[2, pind], markersize =20, label="Mars Perihelion", marker=:star6, color=cn[3])
    scatter!(ax1, marpos[1, aind], marpos[2, aind], markersize =20, label="Mars Aphelion", marker=:star6, color=cn[5])

    date0 = DateTime(2022, 06, 12, 0, 0, 0)
    marp = MarPos(date0)
    scatter!(ax1, marp[1], marp[2], markersize =20, label="Mars", color=cn[2])
    text!(ax1, marp[1]-0.08, marp[2], text = Dates.format.(date0, "yyyy-mm-dd"), align = (:right, :center), space = :data, fontsize = 20, color=cn[2])
    text!(ax1, marp[1]+0.08, marp[2], text = "Mars", align = (:left, :center), space = :data, fontsize = 20, color=cn[2])
    
    
    earthp = EarthPos(date0)
    scatter!(ax1, earthp[1], earthp[2], markersize =20, label="Earth", color=cn[1])
    text!(ax1, earthp[1]-0.08, earthp[2], text = Dates.format.(date0, "yyyy-mm-dd"), align = (:right, :center), space = :data, fontsize = 20, color=cn[1])
    text!(ax1, earthp[1]+0.08, earthp[2], text = "Earth", align = (:left, :center), space = :data, fontsize = 20, color=cn[1])
    
    @show rad2deg(angle(marp[1]+marp[2]*im))-rad2deg(angle(earthp[1]+earthp[2]*im))
    
    date0 = DateTime(2022, 07, 15, 0, 0, 0)
    marp = MarPos(date0)
    scatter!(ax1, marp[1], marp[2], markersize =20, label="Mars", color=cn[2])
    text!(ax1, marp[1]-0.08, marp[2], text = Dates.format.(date0, "yyyy-mm-dd"), align = (:right, :center), space = :data, fontsize = 20, color=cn[2])

    earthp = EarthPos(date0)
    scatter!(ax1, earthp[1], earthp[2], markersize =20, label="Earth", color=cn[1])
    text!(ax1, earthp[1]-0.08, earthp[2], text = Dates.format.(date0, "yyyy-mm-dd"), align = (:right, :center), space = :data, fontsize = 20, color=cn[1])
    text!(ax1, earthp[1]+0.08, earthp[2], text = "Earth", align = (:left, :center), space = :data, fontsize = 20, color=cn[1])

    @show rad2deg(angle(marp[1]+marp[2]*im))-rad2deg(angle(earthp[1]+earthp[2]*im))
    
    date0 = DateTime(2022, 05, 21, 0, 0, 0)
    marp = MarPos(date0)
    scatter!(ax1, marp[1], marp[2], markersize =20, label="Mars", color=cn[2])
    text!(ax1, marp[1]-0.08, marp[2], text = Dates.format.(date0, "yyyy-mm-dd"), align = (:right, :center), space = :data, fontsize = 20, color=cn[2])

    earthp = EarthPos(date0)
    scatter!(ax1, earthp[1], earthp[2], markersize =20, label="Earth", color=cn[1])
    text!(ax1, earthp[1]-0.08, earthp[2], text = Dates.format.(date0, "yyyy-mm-dd"), align = (:right, :center), space = :data, fontsize = 20, color=cn[1])
    text!(ax1, earthp[1]+0.08, earthp[2], text = "Earth", align = (:left, :center), space = :data, fontsize = 20, color=cn[1])

    @show rad2deg(angle(marp[1]+marp[2]*im))-rad2deg(angle(earthp[1]+earthp[2]*im))
    

    text!(ax1, 0.0, -0.2, text = "Sun", align = (:center, :center), space = :data, fontsize = 20, color=:gold)
    #text!(ax1, 0.0, 1.7, text = "Mars orbit", align = (:center, :center), space = :data, fontsize = 20)
    text!(ax1, 1.3, 0.0, text = "Perihelion", align = (:right, :center), space = :data, fontsize = 20, color=cn[3])
    text!(ax1, -1.6, 0.0, text = "Aphelion", align = (:left, :center), space = :data, fontsize = 20, color=cn[5])
    fig    
end 

#Solar cycle
function SolarCycle_plot()
    size_inches = (9, 5)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5,30, 3, 8))
    cn = cgrad(:seaborn_dark, 10, categorical = true)
    lcn = Makie.wong_colors()
    

    mp = MarPos.(sc[!, :time])
    mpr = norm.(mp)

    xd = DateTime(1965, 1, 1, 0, 0, 0) .+ Year(10).*range(0, 6, step=1) 
    xtk = datetime2julian.(xd)
    xtk = Float32.(xtk .- sc[1, :julian])
    ax0 = Axis(fig[1, 1], backgroundcolor=:gray80, xminorticks=IntervalsBetween(10), xminorticksvisible = true, yminorticks=IntervalsBetween(5), yminorticksvisible = true, ylabel=rich("Solar Index F10.7 (sfu)")) #xticks=xtk, yminorticks=IntervalsBetween(5), yminorticksvisible = true, 
    #xlims!(ax1, 200, 5e3)
    ylims!(ax0, 60, 300)
    ax0.xticks = (xtk, Dates.format.(xd, "yyyy"))
    xlims!(ax0, minimum(xtk), maximum(xtk))
    bp=barplot!(ax0, sc[!, :julian] .- sc[1, :julian], sc[!, :f107], color=mpr, colormap=:gist_earth )

    ax0.title=rich(rich("Mariner 6/7/9 ", color = (lcn[2], 0.5)), rich("Mars Express ", color = (lcn[4], 0.5)), rich("Rosetta ", color = (lcn[5], 0.5)), rich("Hubble Space Telescope ", color = (lcn[6], 0.5)), rich("MAVEN ", color = (lcn[3], 0.5)), rich("Tianwen-1", color = (lcn[1], 0.5)))

    #MAVEN
    hspan!([300], [60], xmin=(datetime2julian(DateTime(2014, 11, 1, 0, 0, 0))-sc[1, :julian]-xtk[begin])/(xtk[end]-xtk[begin]), xmax=(datetime2julian(DateTime(2021, 3, 1, 0, 0, 0))-sc[1, :julian]-xtk[begin])/(xtk[end]-xtk[begin]), color= (lcn[3], 0.5),transparency = true) 

    #MEX
    hspan!([300], [60], xmin=(datetime2julian(DateTime(2006, 11, 1, 0, 0, 0))-sc[1, :julian]-xtk[begin])/(xtk[end]-xtk[begin]), xmax=(datetime2julian(DateTime(2009, 11, 1, 0, 0, 0))-sc[1, :julian]-xtk[begin])/(xtk[end]-xtk[begin]), color= (lcn[4], 0.5), transparency = true)

    #HST
    hspan!([300], [60], xmin=(datetime2julian(DateTime(2007, 10, 1, 0, 0, 0))-sc[1, :julian]-xtk[begin])/(xtk[end]-xtk[begin]), xmax=(datetime2julian(DateTime(2007, 11, 1, 0, 0, 0))-sc[1, :julian]-xtk[begin])/(xtk[end]-xtk[begin]), color= (lcn[6], 0.5), transparency = true)

    #HST
    hspan!([300], [60], xmin=(datetime2julian(DateTime(2014, 9, 1, 0, 0, 0))-sc[1, :julian]-xtk[begin])/(xtk[end]-xtk[begin]), xmax=(datetime2julian(DateTime(2014, 11, 1, 0, 0, 0))-sc[1, :julian]-xtk[begin])/(xtk[end]-xtk[begin]), color= (lcn[6], 0.5), transparency = true)

    #Rosetta
    hspan!([300], [60], xmin=(datetime2julian(DateTime(2007, 2, 1, 0, 0, 0))-sc[1, :julian]-xtk[begin])/(xtk[end]-xtk[begin]), xmax=(datetime2julian(DateTime(2007, 3, 1, 0, 0, 0))-sc[1, :julian]-xtk[begin])/(xtk[end]-xtk[begin]), color= (lcn[5], 0.5), transparency = true)

    #Mariner
    hspan!([300], [60], xmin=(datetime2julian(DateTime(1969, 7, 1, 0, 0, 0))-sc[1, :julian]-xtk[begin])/(xtk[end]-xtk[begin]), xmax=(datetime2julian(DateTime(1969, 8, 1, 0, 0, 0))-sc[1, :julian]-xtk[begin])/(xtk[end]-xtk[begin]), color= (lcn[2], 0.5), transparency = true)
    #Mariner
    hspan!([300], [60], xmin=(datetime2julian(DateTime(1971, 12, 1, 0, 0, 0))-sc[1, :julian]-xtk[begin])/(xtk[end]-xtk[begin]), xmax=(datetime2julian(DateTime(1972, 1, 1, 0, 0, 0))-sc[1, :julian]-xtk[begin])/(xtk[end]-xtk[begin]), color= (lcn[2], 0.5), transparency = true)

    #Tianwen-1
    hspan!([300], [60], xmin=(datetime2julian(DateTime(2022, 6, 1, 0, 0, 0))-sc[1, :julian]-xtk[begin])/(xtk[end]-xtk[begin]), xmax=(datetime2julian(DateTime(2022, 7, 1, 0, 0, 0))-sc[1, :julian]-xtk[begin])/(xtk[end]-xtk[begin]), color= (lcn[1], 0.5), transparency = true)

    


    Colorbar(fig[2, 1], vertical = false, bp, label=rich("Martian Heliospheric Distance (AU)"))



    fig

    
end

#solar wind Statistics: new version
function SWSta_plot()
    size_inches = (9, 11.5)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (0,5, 0, -480))
    cn = cgrad(:seaborn_deep, 10, categorical = true)
    lcn = Makie.wong_colors()

    ax1 = Axis(fig[1,1], aspect=DataAspect(), xaxisposition = :top, xlabel=rich(rich("x", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich("(", rich("y", font = :italic), superscript("2"), "+", rich("z", font = :italic), superscript("2"), ")", superscript("1/2"), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10), yticks=range(-10, 10), xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticksvisible = true, yautolimitmargin=(0.00f0, 0.05f0))
    xlims!(ax1, -5, 3) 
    ylims!(ax1, 0, 3.5) 
    hmN = heatmap!(ax1, X, Y, MVN_N, colorscale=identity, colormap=:gist_earth, nan_color=:grey80, colorrange=(1, 8))
    #Colorbar(fig[1,2], vertical = true, hmN, label=rich("Density (cm", superscript("-3"), ")"))
    cb1 = Colorbar(fig[0, 1], vertical = false, flipaxis = false, hmN, label=rich("Density", "\n", "(cm", superscript("-3"), ")"), ticks = range(1,13; step=1),  tellwidth=false, tellheight=false, width=200)
    translate!(cb1.blockscene, -160, -254, 0)


    ax2 = Axis(fig[2,1], aspect=DataAspect(), xaxisposition = :bottom, xlabel=rich(rich("x", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich("(", rich("y", font = :italic), superscript("2"), "+", rich("z", font = :italic), superscript("2"), ")", superscript("1/2"), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10), yticks=range(-10, 10), xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticksvisible = true, yautolimitmargin=(0.05f0, 0.00f0))
    xlims!(ax2, -5, 3) 
    ylims!(ax2, 0, 3.5) 
    hmB = heatmap!(ax2, X, Y, MVN_B, colormap=:gist_earth, nan_color=:grey80, colorrange=(1, 13))
    #Colorbar(fig[2, 2], vertical = true, hmB, label="Magnetic Field Strength (nT)", ticks = range(1,13; step=2))

    cb2 = Colorbar(fig[-1, 1], vertical = false, flipaxis = false, hmB, label="Magnetic Field Strength" * "\n" * "(nT)", ticks = range(1,13; step=2),  tellwidth=false, tellheight=false, width=200)
    translate!(cb2.blockscene, -160, -594, 0)
    



    ax3 = Axis(fig[3,1], aspect=DataAspect(), xlabel=rich(rich("x", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich("(", rich("y", font = :italic), superscript("2"), "+", rich("z", font = :italic), superscript("2"), ")", superscript("1/2"), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10), yticks=range(-10, 10), xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticksvisible = true, yautolimitmargin=(0.00f0, 0.05f0))
    xlims!(ax3, -5, 3) 
    ylims!(ax3, 0, 3.5) 
    hmV = heatmap!(ax3, X, Y, MVN_V,  colormap=:gist_earth, nan_color=:grey80, colorrange=(200, 500))
    #Colorbar(fig[3, 2], vertical = true, hmV, label=rich("Bulk Velocity Magnitude (km⋅s", superscript("-1"), ")"))
    cb3 = Colorbar(fig[-2, 1], vertical = false, flipaxis = false, hmV, label=rich("Bulk Velocity Magnitude" * "\n" * "(km⋅s", superscript("-1"), ")"), tellwidth=false, tellheight=false, width=200)
    translate!(cb3.blockscene, -160, -1058, 0)



    ax = [ax1, ax2, ax3]
    txt = ["A", "B", "C", "D", "E"]
    img = load("marsxz1.png")
    for (i, axi) in enumerate(ax)
        image!(axi, (-1, 1), (-1, 1), rotr90(img))
        lines!(axi, xshock, ryzshock; color=cn[8], linewidth=4, overdraw = true)
        lines!(axi, xmp, ryzmp; color=cn[8], linewidth=4,  overdraw = true)

        lines!(axi, TW_mag2c1hz[!, :Probe_Position_X_MSO], sqrt.(TW_mag2c1hz[!, :Probe_Position_Y_MSO].^2+TW_mag2c1hz[!, :Probe_Position_Z_MSO].^2); color=lcn[1], overdraw = true)
        lines!(axi, MVN_mag2c1hz[!, :Probe_Position_X_MSO], sqrt.(MVN_mag2c1hz[!, :Probe_Position_Y_MSO].^2+MVN_mag2c1hz[!, :Probe_Position_Z_MSO].^2); color=lcn[3], overdraw = true)
        text!(axi, -3.5, 1.5, text = rich(rich("Tianwen-1", color=lcn[1])), align = (:left, :bottom), fontsize = 15)
        text!(axi, -1.0, 1.0, text = rich(rich("MAVEN", color=lcn[3])), align = (:left, :bottom), fontsize = 15)

        for ri in 2:4
            @show ri
            lx, lz, va = stream(3.0, 0.75ri)
            lines!(axi, lx, lz; color=cn[9], overdraw = true, linewidth=2)
            lines!(axi, lx, -lz; color=cn[9], overdraw = true, linewidth=2)
            arrows!(axi, [3.0], [lz[end]], [-1.0], [0.0], arrowsize = 15, lengthscale = 0.5, color=cn[9])
        end

        for ti in 0:4
            local ind = findall(x-> x<date + Dates.Hour(9) + Dates.Minute(30+ti*15) + Dates.Second(1)  && x>date + Dates.Hour(9) + Dates.Minute(30+ti*15) - Dates.Second(1), MVN_mag2c1hz[!, :Time]) 
            local str = Dates.format(date + Dates.Hour(9) + Dates.Minute(30+ti*15), "HH:MM")
            ind = ind[1]
            color = cn[ ti+4 == 8 ? 2 : ti+4]
            scatter!(axi, MVN_mag2c1hz[ind, :Probe_Position_X_MSO], sqrt.(MVN_mag2c1hz[ind, :Probe_Position_Y_MSO].^2+MVN_mag2c1hz[ind, :Probe_Position_Z_MSO].^2), markersize=10, label=str, color=color, transparency=true) 
            if ti == 0 
                text!(axi, 0.2, 2.3, text=str, align = (:right, :center), space = :data, fontsize = 15, color=color, rotation=deg2rad(-40))
            elseif ti == 1
                text!(axi, 0.34, 2.33, text=str, align = (:left, :center), space = :data, fontsize = 15, color=color, rotation=deg2rad(70))
            elseif ti == 2 
                text!(axi, 0.40, 2.20, text=str, align = (:left, :center), space = :data, fontsize = 15, color=color, rotation=deg2rad(40))
            elseif ti == 3 
                text!(axi, 0.40, 2.12, text=str, align = (:left, :center), space = :data, fontsize = 15, color=color, rotation=deg2rad(-10))
            elseif ti == 4 
                text!(axi, 0.38, 1.92, text=str, align = (:left, :center), space = :data, fontsize = 15, color=color, rotation=deg2rad(-10))
            end
            
            local ind = findall(x-> x<date + Dates.Hour(9) + Dates.Minute(30+ti*15) + Dates.Second(1)  && x>date + Dates.Hour(9) + Dates.Minute(30+ti*15) - Dates.Second(1), TW_mag2c1hz[!, :Time]) 
            ind = ind[1]
            local str = Dates.format(date + Dates.Hour(9) + Dates.Minute(30+ti*15), "HH:MM")
            scatter!(axi, TW_mag2c1hz[ind, :Probe_Position_X_MSO], sqrt.(TW_mag2c1hz[ind, :Probe_Position_Y_MSO].^2+TW_mag2c1hz[ind, :Probe_Position_Z_MSO].^2), markersize=10, label=str, color=color, transparency=true)
            text!(axi, TW_mag2c1hz[ind, :Probe_Position_X_MSO]-0.12, sqrt.(TW_mag2c1hz[ind, :Probe_Position_Y_MSO].^2+TW_mag2c1hz[ind, :Probe_Position_Z_MSO].^2)+0.16, text=str, align = (:left, :top), space = :data, fontsize = 15, color=color, rotation=deg2rad(60))
            
        end
    end

    
    txt = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"]
    Label(fig[1, 0], txt[1],font = :bold, valign=:top, tellheight=false, fontsize = 20)
    Label(fig[2, 0], txt[2],font = :bold, valign=:top, tellheight=false, fontsize = 20)
    Label(fig[3, 0], txt[3],font = :bold, valign=:top, tellheight=false, fontsize = 20)


    hidexdecorations!.(ax[begin:end-1], grid=false)
    linkxaxes!(ax...)
    rowsize!(fig.layout, 0, Relative(0.01))
    rowsize!(fig.layout, -1, Relative(0.01))
    #colsize!(fig.layout, 1, Relative(0.9))
    rowgap!(fig.layout, 4)
    colgap!(fig.layout, 0)

    fig
end

function Mach_plot()
    size_inches = (9, 4.9)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5,4, 0, 1))
    cn = cgrad(:seaborn_deep, 10, categorical = true)
    lcn = Makie.wong_colors()
    ax1 = Axis(fig[1,1], aspect=DataAspect(), xlabel=rich(rich("x", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich("(", rich("y", font = :italic), superscript("2"), "+", rich("z", font = :italic), superscript("2"), ")", superscript("1/2"), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10), yticks=range(-10, 10), yautolimitmargin=(0.00f0, 0.05f0))
    xlims!(ax1, -5, 3) 
    ylims!(ax1, 0, 4) 

    hmMₐ = heatmap!(X, Y, MVN_Mₐ, colormap=:gist_earth, nan_color=:grey80, colorrange=(1, 11))
    Colorbar(fig, hmMₐ, label = "Alfven Mach Number", 
    labelsize = 15, ticklabelsize = 15, bbox=ax1.scene.viewport, vertical=false, flipaxis=false, width=200, 
    alignmode = Outside(3), halign = :left, valign = :bottom)


    img = load("marsxz1.png")
    image!(ax1, (-1, 1), (-1, 1), rotr90(img))
    lines!(ax1, xshock, ryzshock; color=cn[8], linewidth=4, overdraw = true)
    lines!(ax1, xmp, ryzmp; color=cn[8], linewidth=4,  overdraw = true)
    lines!(ax1, TW_mag2c1hz[!, :Probe_Position_X_MSO], sqrt.(TW_mag2c1hz[!, :Probe_Position_Y_MSO].^2+TW_mag2c1hz[!, :Probe_Position_Z_MSO].^2); color=cn[1], overdraw = true)
    lines!(ax1, MVN_mag2c1hz[!, :Probe_Position_X_MSO], sqrt.(MVN_mag2c1hz[!, :Probe_Position_Y_MSO].^2+MVN_mag2c1hz[!, :Probe_Position_Z_MSO].^2); color=cn[3], overdraw = true)

    for ri in 2:4
        @show ri
        lx, lz, va = stream(3.0, 0.75ri)
        lines!(ax1, lx, lz; color=cn[9], overdraw = true, linewidth=2)
        lines!(ax1, lx, -lz; color=cn[9], overdraw = true, linewidth=2)
        arrows!(ax1, [3.0], [lz[end]], [-1.0], [0.0], arrowsize = 15, lengthscale = 0.5, color=cn[9])
    end

    for ti in 0:4
        local ind = findall(x-> x<date + Dates.Hour(9) + Dates.Minute(30+ti*15) + Dates.Second(1)  && x>date + Dates.Hour(9) + Dates.Minute(30+ti*15) - Dates.Second(1), MVN_mag2c1hz[!, :Time]) 
        local str = Dates.format(date + Dates.Hour(9) + Dates.Minute(30+ti*15), "HH:MM")
        ind = ind[1]
        color = cn[ ti+4 == 8 ? 2 : ti+4]
        for axi in [ax1]
            scatter!(axi, MVN_mag2c1hz[ind, :Probe_Position_X_MSO], sqrt.(MVN_mag2c1hz[ind, :Probe_Position_Y_MSO].^2+MVN_mag2c1hz[ind, :Probe_Position_Z_MSO].^2), markersize=15, label=str, color=color, transparency=true)   
        end
        if ti == 0 
            text!(ax1, 0.2, 2.3, text=str, align = (:right, :center), space = :data, fontsize = 15, color=color, rotation=deg2rad(-40))
        elseif ti == 1
            text!(ax1, 0.34, 2.33, text=str, align = (:left, :center), space = :data, fontsize = 15, color=color, rotation=deg2rad(70))
        elseif ti == 2 
            text!(ax1, 0.40, 2.20, text=str, align = (:left, :center), space = :data, fontsize = 15, color=color, rotation=deg2rad(40))
        elseif ti == 3 
            text!(ax1, 0.40, 2.12, text=str, align = (:left, :center), space = :data, fontsize = 15, color=color, rotation=deg2rad(-10))
        elseif ti == 4 
            text!(ax1, 0.38, 1.92, text=str, align = (:left, :center), space = :data, fontsize = 15, color=color, rotation=deg2rad(-10))
        end
        local ind = findall(x-> x<date + Dates.Hour(9) + Dates.Minute(30+ti*15) + Dates.Second(1)  && x>date + Dates.Hour(9) + Dates.Minute(30+ti*15) - Dates.Second(1), TW_mag2c1hz[!, :Time]) 
        ind = ind[1]
        local str = Dates.format(date + Dates.Hour(9) + Dates.Minute(30+ti*15), "HH:MM")
        for axi in [ax1]
            scatter!(axi, TW_mag2c1hz[ind, :Probe_Position_X_MSO], sqrt.(TW_mag2c1hz[ind, :Probe_Position_Y_MSO].^2+TW_mag2c1hz[ind, :Probe_Position_Z_MSO].^2), markersize=15, label=str, color=color, transparency=true)
        end
        text!(ax1, TW_mag2c1hz[ind, :Probe_Position_X_MSO]-0.12, sqrt.(TW_mag2c1hz[ind, :Probe_Position_Y_MSO].^2+TW_mag2c1hz[ind, :Probe_Position_Z_MSO].^2)+0.16, text=str, align = (:left, :top), space = :data, fontsize = 15, color=color, rotation=deg2rad(60))
        
    end

    text!(ax1, -3.5, 1.5, text = rich(rich("Tianwen-1", color=lcn[1])), align = (:left, :bottom), fontsize = 15)
    text!(ax1, -1.0, 1.0, text = rich(rich("MAVEN", color=lcn[3])), align = (:left, :bottom), fontsize = 15)

    fig
end

#solar wind Statistics
function SWStaBackup_plot()
    size_inches = (9, 6.5)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5,16, 3, 8))
    cn = cgrad(:seaborn_deep, 10, categorical = true)

    ax1 = Axis(fig[2,1], aspect=DataAspect(), xaxisposition = :top, xlabel=rich(rich("x", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich("(", rich("y", font = :italic), superscript("2"), "+", rich("z", font = :italic), superscript("2"), ")", superscript("1/2"), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10), yticks=range(-10, 10), xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticksvisible = true, yautolimitmargin=(0.00f0, 0.05f0))
    xlims!(ax1, -3, 3) 
    ylims!(ax1, 0, 3) 
    hmN = heatmap!(ax1, X, Y, MVN_N, colorscale=identity, colormap=:gist_earth, nan_color=:grey80, colorrange=(1, 8))
    Colorbar(fig[1,1], vertical = false, hmN, label=rich("Density (cm", superscript("-3"), ")"))


    ax2 = Axis(fig[3,1], aspect=DataAspect(), xaxisposition = :bottom, xlabel=rich(rich("x", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich("(", rich("y", font = :italic), superscript("2"), "+", rich("z", font = :italic), superscript("2"), ")", superscript("1/2"), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10), yticks=range(-10, 10), xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticksvisible = true, yautolimitmargin=(0.05f0, 0.00f0))
    xlims!(ax2, -3, 3) 
    ylims!(ax2, 0, 3) 
    hmB = heatmap!(ax2, X, Y, MVN_B, colormap=:gist_earth, nan_color=:grey80, colorrange=(1, 13))
    #Colorbar(fig[4, 1], vertical = false, flipaxis = false, hmB, label="Magnetic Field Strength (nT)", ticks = range(1,13; step=2))
    cb = Colorbar(fig[4, 1], vertical = true, flipaxis = false, hmB, label="Magnetic Field Strength" * "\n" * "(nT)", ticks = range(1,13; step=2),  tellwidth=false, tellheight=false, height=70)
    translate!(cb.blockscene, -50, 60, 0)



    ax3 = Axis(fig[2,3], aspect=DataAspect(), xaxisposition = :top, xlabel=rich(rich("x", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich("(", rich("y", font = :italic), superscript("2"), "+", rich("z", font = :italic), superscript("2"), ")", superscript("1/2"), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10), yticks=range(-10, 10), xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticksvisible = true, yautolimitmargin=(0.00f0, 0.05f0))
    xlims!(ax3, -3, 3) 
    ylims!(ax3, 0, 3) 
    hmV = heatmap!(ax3, X, Y, MVN_V,  colormap=:gist_earth, nan_color=:grey80, colorrange=(200, 500))
    Colorbar(fig[1, 3], vertical = false, hmV, label=rich("Bulk Velocity Magnitude (km⋅s", superscript("-1"), ")"))



    ax4 = Axis(fig[3,3], aspect=DataAspect(), xaxisposition = :bottom, xlabel=rich(rich("x", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich("(", rich("y", font = :italic), superscript("2"), "+", rich("z", font = :italic), superscript("2"), ")", superscript("1/2"), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10), yticks=range(-10, 10), xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticksvisible = true, yautolimitmargin=(0.05f0, 0.00f0))
    xlims!(ax4, -3, 3) 
    ylims!(ax4, 0, 3) 
    hmMₐ = heatmap!(ax4, X, Y, MVN_Mₐ, colormap=:gist_earth, nan_color=:grey80, colorrange=(1, 11))
    Colorbar(fig[4,3], vertical = false, hmMₐ, flipaxis = false, label="Alfven Mach Number", ticks = range(1,11; step=2))
	#, width=Relative(6/7)

    ax = [ax1, ax2, ax3, ax4]
    txt = ["A", "B", "C", "D", "E"]
    img = load("marsxz1.png")
    for (i, axi) in enumerate(ax)
        image!(axi, (-1, 1), (-1, 1), rotr90(img))
        lines!(axi, xshock, ryzshock; color=cn[8], linewidth=4, overdraw = true)
        lines!(axi, xmp, ryzmp; color=cn[8], linewidth=4,  overdraw = true)
    end

    txt = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"]
    Label(fig[1:2, 0], txt[1],font = :bold, valign=:top, tellheight=false, fontsize = 20)
    Label(fig[3:4, 0], txt[2],font = :bold, valign=:top, tellheight=false, fontsize = 20)
    Label(fig[1:2, 2], txt[3],font = :bold, valign=:top, tellheight=false, fontsize = 20)
    Label(fig[3:4, 2], txt[4],font = :bold, valign=:top, tellheight=false, fontsize = 20)

    hidexdecorations!(ax1, grid=false)
    hidexdecorations!(ax3, grid=false)
    hideydecorations!(ax3, grid=false)
    hideydecorations!(ax4, grid=false)
    linkyaxes!(ax1, ax3)
    linkyaxes!(ax2, ax4)

    rowgap!(fig.layout, 1)
    #g = fig[2:3, 1:2] = GridLayout()
    #rowgap!(g, 0.5)
	#colgap!(fig.layout, 25.0)

    fig
end

#solar wind streamline
function streamline_plot()
    size_inches = (9, 6)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5,16, 3, 8))
    cn = cgrad(:seaborn_deep, 10, categorical = true)
    #=
    g1 = fig[1, 1] = GridLayout()
    g2 = fig[1, 2] = GridLayout()
    g3 = fig[1, 3] = GridLayout()
    =#
    ax1 = Axis(fig[2,1], aspect=DataAspect(), xlabel=rich(rich("x", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich("(", rich("y", font = :italic), superscript("2"), "+", rich("z", font = :italic), superscript("2"), ")", superscript("1/2"), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10), yticks=range(-10, 10), yautolimitmargin=(0.00f0, 0.05f0))
    xlims!(ax1, -5, 3) 
    ylims!(ax1, 0, 4) 

    

    hmVA = heatmap!(ax1, X, Y, MVN_VA, colormap=:gist_earth, nan_color=:grey80, colorrange=(140, 180))

    img = load("marsxz1.png")
    image!(ax1, (-1, 1), (-1, 1), rotr90(img))
    lines!(ax1, xshock, ryzshock; color=cn[8], linewidth=4, overdraw = true)
    lines!(ax1, xmp, ryzmp; color=cn[8], linewidth=4,  overdraw = true)
    lines!(ax1, TW_mag2c1hz[!, :Probe_Position_X_MSO], sqrt.(TW_mag2c1hz[!, :Probe_Position_Y_MSO].^2+TW_mag2c1hz[!, :Probe_Position_Z_MSO].^2); color=cn[1], overdraw = true)
    lines!(ax1, MVN_mag2c1hz[!, :Probe_Position_X_MSO], sqrt.(MVN_mag2c1hz[!, :Probe_Position_Y_MSO].^2+MVN_mag2c1hz[!, :Probe_Position_Z_MSO].^2); color=cn[3], overdraw = true)

    for ri in 2:4
        @show ri
        lx, lz, va = stream(3.0, 0.75ri)
        lines!(ax1, lx, lz; color=cn[9], overdraw = true, linewidth=2)
        lines!(ax1, lx, -lz; color=cn[9], overdraw = true, linewidth=2)
        arrows!(ax1, [3.0], [lz[end]], [-1.0], [0.0], arrowsize = 15, lengthscale = 0.5, color=cn[9])
    end

    for ti in 0:4
        local ind = findall(x-> x<date + Dates.Hour(9) + Dates.Minute(30+ti*15) + Dates.Second(1)  && x>date + Dates.Hour(9) + Dates.Minute(30+ti*15) - Dates.Second(1), MVN_mag2c1hz[!, :Time]) 
        local str = Dates.format(date + Dates.Hour(9) + Dates.Minute(30+ti*15), "HH:MM")
        ind = ind[1]
        color = cn[ ti+4 == 8 ? 2 : ti+4]
        for axi in [ax1]
            scatter!(axi, MVN_mag2c1hz[ind, :Probe_Position_X_MSO], sqrt.(MVN_mag2c1hz[ind, :Probe_Position_Y_MSO].^2+MVN_mag2c1hz[ind, :Probe_Position_Z_MSO].^2), markersize=15, label=str, color=color, transparency=true)   
        end
        if ti == 0 
            text!(ax1, 0.2, 2.4, text=str, align = (:left, :center), space = :data, fontsize = 20, color=color, rotation=deg2rad(60))
        elseif ti == 1
            text!(ax1, 0.35, 2.25, text=str, align = (:left, :center), space = :data, fontsize = 20, color=color, rotation=deg2rad(60))
        elseif ti == 2 
            text!(ax1, 0.50, 2.10, text=str, align = (:left, :center), space = :data, fontsize = 20, color=color, rotation=deg2rad(60))
        elseif ti == 3 
            text!(ax1, 0.26, 2.15, text=str, align = (:right, :center), space = :data, fontsize = 20, color=color, rotation=deg2rad(00))
        elseif ti == 4 
            text!(ax1, 0.25, 1.95, text=str, align = (:right, :center), space = :data, fontsize = 20, color=color, rotation=deg2rad(0))
        end
        local ind = findall(x-> x<date + Dates.Hour(9) + Dates.Minute(30+ti*15) + Dates.Second(1)  && x>date + Dates.Hour(9) + Dates.Minute(30+ti*15) - Dates.Second(1), TW_mag2c1hz[!, :Time]) 
        ind = ind[1]
        local str = Dates.format(date + Dates.Hour(9) + Dates.Minute(30+ti*15), "HH:MM")
        for axi in [ax1]
            scatter!(axi, TW_mag2c1hz[ind, :Probe_Position_X_MSO], sqrt.(TW_mag2c1hz[ind, :Probe_Position_Y_MSO].^2+TW_mag2c1hz[ind, :Probe_Position_Z_MSO].^2), markersize=15, label=str, color=color, transparency=true)
        end
        text!(ax1, TW_mag2c1hz[ind, :Probe_Position_X_MSO]+0.15, sqrt.(TW_mag2c1hz[ind, :Probe_Position_Y_MSO].^2+TW_mag2c1hz[ind, :Probe_Position_Z_MSO].^2)-0.08, text=str, align = (:left, :top), space = :data, fontsize = 20, color=color, rotation=deg2rad(-80))
        
    end

    text!(ax1, -2.3, 2.8, text="Tianwen-1", align = (:center, :bottom), space = :data, fontsize = 20, color=cn[1],rotation=deg2rad(10))
    text!(ax1, 1.0, 3.0, text="MAVEN", align = (:center, :center), space = :data, fontsize = 20, color=cn[3], rotation=deg2rad(-40))

    Colorbar(fig[1, 1], vertical=false, flipaxis = true,  hmVA, label="Solar Wind Bulk Velocity Elevation Angle (Deg)")

    fig
end

#Longterm solar wind
function SolarWind_OverviewBackup_plot()
    size_inches = (9, 10)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5,30, 3, 8))
    cn = cgrad(:seaborn_dark, 10, categorical = true)
    lcn = Makie.wong_colors()
    xd = date .+ Day(1).*range(-15, 3, step=3)
    #xd = date .+ Hour(1).*[-15, -13]
    #xd = DateTime(2022, 6, 11, 10, 0, 0) .+ Minute(1).*range(-30, 30, step=15)
    xtk = datetime2julian.(xd)
    xtk = Float32.(xtk .- MVN_kp[1, :JulUT])
    

    ind = findall(x->x>1.2, MVN_kp[!, 197].*1u"km"./R_Mars )
    #ind = findall(x->x<50, MVN_kp[!, 49])
    #ind = findall(x->x>680.0, sqrt.(MVN_kp[:, 43].^2 .+ MVN_kp[:, 45].^2 .+ MVN_kp[:, 47].^2))
    
    
    ax1 = Axis(fig[1, 1], backgroundcolor=:gray80, xminorticks=IntervalsBetween(3), xminorticksvisible = true, yminorticks=IntervalsBetween(5), yminorticksvisible = true, ylabel=rich("Bulk Velocity (km⋅s", superscript("-1"), ")")) #xticks=xtk, yminorticks=IntervalsBetween(5), yminorticksvisible = true, 
    #xlims!(ax1, 200, 5e3)
    ylims!(ax1, 100, 700)
    ax1.xticks = (xtk, Dates.format.(xd, "u dd"))
    xlims!(ax1, minimum(xtk), maximum(xtk))
    #scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[ind, 43])
    #scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[ind, 45])
    #scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[ind, 47])
    scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], sqrt.(MVN_kp[ind, 43].^2 .+ MVN_kp[ind, 45].^2 .+ MVN_kp[ind, 47].^2), marker=:rect, markersize=4, color=lcn[3])
    #scatter!(MVN_kp[ind, 197]./R_Mars*1e3 .+1, sqrt.(MVN_kp[ind, 43].^2 .+ MVN_kp[ind, 45].^2 .+ MVN_kp[ind, 47].^2), color=lcn[1])


    ax2 = Axis(fig[2, 1], backgroundcolor=:gray80, xminorticks=IntervalsBetween(3), xminorticksvisible = true, yminorticks=IntervalsBetween(5), yminorticksvisible = true, ylabel=rich("Proton Density (cm", superscript("-3"), ")")) #xticks=xtk, yminorticks=IntervalsBetween(5), yminorticksvisible = true, 
    #xlims!(ax1, 200, 5e3)
    #
    ax2.xticks = (xtk, Dates.format.(xd, "dd"))
    xlims!(ax2, minimum(xtk), maximum(xtk))
    ylims!(ax2, 0, 100)
    #scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[ind, 2])
    #scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[ind, 23])
    scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], abs.(MVN_kp[ind, 41]), marker=:rect, markersize=4, color=lcn[3])
    #scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], abs.(MVN_kp[ind, 54]))
    #scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], abs.(MVN_kp[ind, 56]))
    #scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], abs.(MVN_kp[ind, 58]))

    ax3 = Axis(fig[3, 1], backgroundcolor=:gray80, xminorticks=IntervalsBetween(3), xminorticksvisible = true, yminorticks=IntervalsBetween(5), yminorticksvisible = true, ylabel=rich("Proton Temperature (eV)"))#xticks=xtk, yminorticks=IntervalsBetween(5), yminorticksvisible = true, 
    #xlims!(ax1, 200, 5e3)
    #
    ax3.xticks = (xtk, Dates.format.(xd, "dd"))
    xlims!(ax3, minimum(xtk), maximum(xtk))
    ylims!(ax3, 0, 1000)
    #scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[ind, 2])
    #scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[ind, 23])
    scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[ind, 49], marker=:rect, markersize=4, color=lcn[3])


    ax4 = Axis(fig[4, 1], backgroundcolor=:gray80, xminorticks=IntervalsBetween(3), xminorticksvisible = true, yminorticks=IntervalsBetween(5), yminorticksvisible = true, xlabel="Universal Time", ylabel=rich("Magnetic Strength (nT)"))
    #xlims!(ax1, 200, 5e3)
    ylims!(ax4, 0, 40)
    ax4.xticks = (xtk, Dates.format.(xd, "u dd"))
    xlims!(ax4, minimum(xtk), maximum(xtk))
    #scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[ind, 128])
    #scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[ind, 130])
    #scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[ind, 132])
    scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], sqrt.(MVN_kp[ind, 128].^2 .+ MVN_kp[ind, 130].^2 .+ MVN_kp[ind, 132].^2), marker=:rect, markersize=4, color=lcn[3])

    ax = [ax1, ax2, ax3, ax4]
    txt = ["A", "B", "C", "D", "E", "F"]
    for (i, axi) in enumerate(ax)
        text!(axi, 0.0, 1.0, text = txt[i], font = :bold, align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 20)
    end
    
    linkxaxes!(ax1, ax2, ax3, ax4)
    hidexdecorations!(ax1, grid=false)
    hidexdecorations!(ax2, grid=false)
    hidexdecorations!(ax3, grid=false)
    fig
end

function SolarWind_Overview_plot()
    size_inches = (9, 10)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5,25, 3, 8))
    cn = cgrad(:seaborn_dark, 10, categorical = true)
    lcn = Makie.wong_colors()
    date = DateTime(2022, 06, 12, 00, 00, 00)
    xd = date .+ Day(1).*range(-22, 3, step=5)
    #xd = date .+ Hour(1).*[-15, -13]
    #xd = DateTime(2022, 6, 12, 0, 0, 0) .+ Hour(2).*range(0, 24, step=4)
    xtk = datetime2julian.(xd)
    xtk = Float32.(xtk .- MVN_kp[1, :JulUT])

    ind = findall(x->x>1.285, MVN_kp[!, 197].*1u"km"./R_Mars )

    dt = 6.0
   
    ax1 = Axis(fig[1, 1], backgroundcolor=:gray80, xminorticks=IntervalsBetween(3), xminorticksvisible = true, yminorticks=IntervalsBetween(4), yminorticksvisible = true) 
    ylims!(ax1, 140, 700)
    ax1.xticks = (xtk, Dates.format.(xd, "u dd"))
    xlims!(ax1, minimum(xtk), maximum(xtk))
    lines!(ax1, WIND_swejulian[WIND_sweind].-MVN_kp[1, :JulUT].+dt, WIND_sweVp[WIND_sweind], color=lcn[2])
    scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], sqrt.(MVN_kp[ind, 43].^2 .+ MVN_kp[ind, 45].^2 .+ MVN_kp[ind, 47].^2), marker=:rect, markersize=4, color=lcn[3])
    
    ax2 = Axis(fig[2, 1], backgroundcolor=:gray80, yscale = log10, xminorticks=IntervalsBetween(9), xminorticksvisible = true, yminorticks=IntervalsBetween(9), yminorticksvisible = true) 
    ax2.xticks = (xtk, Dates.format.(xd, "dd"))
    ax2.yticks = ([0.1, 1, 10, 100], ["0.1", "1", "10", "100"])
    xlims!(ax2, minimum(xtk), maximum(xtk))
    ylims!(ax2, 1, 80)
    lines!(ax2, WIND_swejulian[WIND_sweind].-MVN_kp[1, :JulUT].+dt, WIND_sweNp[WIND_sweind], color=lcn[2])
    scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], abs.(MVN_kp[ind, 41]), marker=:rect, markersize=4, color=lcn[3])

    ax3 = Axis(fig[3, 1], backgroundcolor=:gray80, yscale = log10, xminorticks=IntervalsBetween(3), xminorticksvisible = true, yminorticks=IntervalsBetween(9), yminorticksvisible = true)
    ax3.xticks = (xtk, Dates.format.(xd, "dd"))
    ax3.yticks = ([0.1, 1, 10, 100, 1000], ["0.1", "1", "10", "100", "1000"])
    xlims!(ax3, minimum(xtk), maximum(xtk))
    ylims!(ax3, 0.7, 900)
    lines!(ax3, WIND_swejulian[WIND_sweind].-MVN_kp[1, :JulUT].+dt, WIND_sweTp[WIND_sweind], color=lcn[2])
    scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[ind, 49], marker=:rect, markersize=4, color=lcn[3])

    ax4 = Axis(fig[4, 1], backgroundcolor=:gray80, yscale = log10, xminorticks=IntervalsBetween(8), xminorticksvisible = true, yminorticks=IntervalsBetween(9), yminorticksvisible = true, xlabel="Universal Time")
    ylims!(ax4, 1, 40)
    ax4.xticks = (xtk, Dates.format.(xd, "u dd"))
    ax4.yticks = [1, 10, 100]
    xlims!(ax4, minimum(xtk), maximum(xtk))
    lines!(ax4, WIND_mfijulian[WIND_mfiind].-MVN_kp[1, :JulUT].+dt, WIND_mfiBt[WIND_mfiind], color=lcn[2])
    scatter!(MVN_kp[ind, :JulUT].-MVN_kp[1, :JulUT], sqrt.(MVN_kp[ind, 128].^2 .+ MVN_kp[ind, 130].^2 .+ MVN_kp[ind, 132].^2), marker=:rect, markersize=4, color=lcn[3])

    ind = [1, 2, 3, 4]
    txt = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"]
    for (i, indi) in enumerate(ind)
        Label(fig[indi, -1], txt[i],font = :bold, valign=:top, tellheight=false, fontsize = 20)
    end
    Label(fig[1, 0], rich("Bulk Velocity" *"\n" * "(km⋅s", superscript("-1"), ")"), rotation = pi/2, tellheight=false)
    Label(fig[2, 0], rich("Proton Density" *"\n" * "(cm", superscript("-3"), ")"), rotation = pi/2, tellheight=false)
    Label(fig[3, 0], rich("Proton Temperature" *"\n" * "(eV)"), rotation = pi/2, tellheight=false)
    Label(fig[4, 0], rich("Magnetic Strength" *"\n" * "(nT)"), rotation = pi/2, tellheight=false)

    ax = [ax1, ax2, ax3, ax4]
    linkxaxes!(ax...)
    hidexdecorations!.(ax[1:end-1], grid=false)
    rowgap!(fig.layout, 3)
    colgap!(fig.layout, 2)
    fig
end

# orbit xyz
function orbitxyz_plot()
    size_inches = (9, 7)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5,16, 3, 8))
    cn = cgrad(:seaborn_deep, 10, categorical = true)

    ax1 = Axis(fig[1,1], aspect=DataAspect(), xlabel=rich(rich("x", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich(rich("z", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10), yticks=range(-10, 10), yautolimitmargin=(0.00f0, 0.05f0), xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticksvisible = true, backgroundcolor=:gray80)
    xlims!(ax1, -4.5, 1.5) 
    ylims!(ax1, -4, 4) 

    img = load("marsxz1.png")
    image!(ax1, (-1, 1), (-1, 1), rotr90(img))
    
    lines!(ax1, TW_mag2c1hz[!, :Probe_Position_X_MSO], TW_mag2c1hz[!, :Probe_Position_Z_MSO]; color=cn[1], overdraw = true)

    ind1 = findall(x->x>0 , MVN_mag2c1hz[!, :Probe_Position_Y_MSO])
    r = @. sqrt(MVN_mag2c1hz[!, :Probe_Position_X_MSO]^2 + MVN_mag2c1hz[!, :Probe_Position_Z_MSO]^2)
    ind2 = findall(x->x<1 , r)
    ind = intersect(ind1, ind2)
    tx = MVN_mag2c1hz[:, :Probe_Position_X_MSO]
    tz = MVN_mag2c1hz[:, :Probe_Position_Z_MSO]
    tx[ind] .= NaN64
    tz[ind] .= NaN64
    lines!(ax1, tx, tz; color=cn[3], overdraw = true)

    for ti in 0:4
        color = cn[ ti+4 == 8 ? 2 : ti+4]
        local ind = findall(x-> x<date + Dates.Hour(9) + Dates.Minute(30+ti*15) + Dates.Second(1)  && x>date + Dates.Hour(9) + Dates.Minute(30+ti*15) - Dates.Second(1), MVN_mag2c1hz[!, :Time]) 
        local str = Dates.format(date + Dates.Hour(9) + Dates.Minute(30+ti*15), "HH:MM")
        ind = ind[1]
        for axi in [ax1]
            scatter!(axi, MVN_mag2c1hz[ind, :Probe_Position_X_MSO], MVN_mag2c1hz[ind, :Probe_Position_Z_MSO], markersize=15, label=str, color=color, transparency=true)   
        end
        if ti == 0 
            text!(ax1, MVN_mag2c1hz[ind, :Probe_Position_X_MSO]-0.15, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.15, text=str, align = (:right, :top), space = :data, fontsize = 20, color=color, rotation=deg2rad(20))
        elseif ti == 1
            text!(ax1, MVN_mag2c1hz[ind, :Probe_Position_X_MSO]-0.15, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.03, text=str, align = (:right, :top), space = :data, fontsize = 20, color=color, rotation=deg2rad(20))
        elseif ti == 2
            text!(ax1, MVN_mag2c1hz[ind, :Probe_Position_X_MSO]+0.15, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.01, text=str, align = (:left, :center), space = :data, fontsize = 20, color=color, rotation=deg2rad(20))
        elseif ti == 3
            text!(ax1, MVN_mag2c1hz[ind, :Probe_Position_X_MSO]+0.15, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.09, text=str, align = (:left, :center), space = :data, fontsize = 20, color=color, rotation=deg2rad(20))
        else
            text!(ax1, MVN_mag2c1hz[ind, :Probe_Position_X_MSO]+0.15, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.05, text=str, align = (:left, :center), space = :data, fontsize = 20, color=color, rotation=deg2rad(20))
        end
        local ind = findall(x-> x<date + Dates.Hour(9) + Dates.Minute(30+ti*15) + Dates.Second(1)  && x>date + Dates.Hour(9) + Dates.Minute(30+ti*15) - Dates.Second(1), TW_mag2c1hz[!, :Time]) 
        ind = ind[1]
        local str = Dates.format(date + Dates.Hour(9) + Dates.Minute(30+ti*15), "HH:MM")
        for axi in [ax1]
            scatter!(axi, TW_mag2c1hz[ind, :Probe_Position_X_MSO], TW_mag2c1hz[ind, :Probe_Position_Z_MSO], markersize=15, label=str, color=color, transparency=true)
        end
        text!(ax1, TW_mag2c1hz[ind, :Probe_Position_X_MSO]+0.15, TW_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.08, text=str, align = (:left, :top), space = :data, fontsize = 20, color=color, rotation=deg2rad(-80))
        
    end
    text!(ax1, -2.3, 2.8, text="Tianwen-1", align = (:center, :bottom), space = :data, fontsize = 20, color=cn[1],rotation=deg2rad(10))
    text!(ax1, 0.5, -2.9, text="MAVEN", align = (:center, :bottom), space = :data, fontsize = 20, color=cn[3], rotation=deg2rad(20))



    ax2 = Axis(fig[1,2], aspect=DataAspect(), xlabel=rich(rich("y", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich(rich("z", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10), yticks=range(-10, 10), yautolimitmargin=(0.00f0, 0.05f0), xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticksvisible = true, backgroundcolor=:gray80)
    xlims!(ax2, -2, 2) 
    ylims!(ax2, -4, 4) 

    img = load("marsyz1.png")
    image!(ax2, (-1, 1), (-1, 1), rotr90(img))
    
    ind1 = findall(x->x<0 , TW_mag2c1hz[!, :Probe_Position_X_MSO])
    r = @. sqrt(TW_mag2c1hz[!, :Probe_Position_Y_MSO]^2 + TW_mag2c1hz[!, :Probe_Position_Z_MSO]^2)
    ind2 = findall(x->x<1 , r)
    ind = intersect(ind1, ind2)
    tx = TW_mag2c1hz[:, :Probe_Position_Y_MSO]
    tz = TW_mag2c1hz[:, :Probe_Position_Z_MSO]
    tx[ind] .= NaN64
    tz[ind] .= NaN64
    lines!(ax2, tx, tz; color=cn[1], overdraw = true)

 
    tx = MVN_mag2c1hz[:, :Probe_Position_Y_MSO]
    tz = MVN_mag2c1hz[:, :Probe_Position_Z_MSO]
    lines!(ax2, tx, tz; color=cn[3], overdraw = true)

    for ti in 0:4
        color = cn[ ti+4 == 8 ? 2 : ti+4]
        local ind = findall(x-> x<date + Dates.Hour(9) + Dates.Minute(30+ti*15) + Dates.Second(1)  && x>date + Dates.Hour(9) + Dates.Minute(30+ti*15) - Dates.Second(1), MVN_mag2c1hz[!, :Time]) 
        local str = Dates.format(date + Dates.Hour(9) + Dates.Minute(30+ti*15), "HH:MM")
        ind = ind[1]
        for axi in [ax2]
            scatter!(axi, MVN_mag2c1hz[ind, :Probe_Position_Y_MSO], MVN_mag2c1hz[ind, :Probe_Position_Z_MSO], markersize=15, label=str, color=color, transparency=true)   
        end
        if ti == 0 
            text!(ax2, MVN_mag2c1hz[ind, :Probe_Position_Y_MSO]-0.00, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.1, text=str, align = (:left, :center), space = :data, fontsize = 20, color=color, rotation=deg2rad(-85))
        elseif ti == 1
            text!(ax2, MVN_mag2c1hz[ind, :Probe_Position_Y_MSO]-0.00, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.1, text=str, align = (:left, :center), space = :data, fontsize = 20, color=color, rotation=deg2rad(-85))
        else
            text!(ax2, MVN_mag2c1hz[ind, :Probe_Position_Y_MSO]+0.00, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.1, text=str, align = (:left, :center), space = :data, fontsize = 20, color=color, rotation=deg2rad(-85))
        end
        local ind = findall(x-> x<date + Dates.Hour(9) + Dates.Minute(30+ti*15) + Dates.Second(1)  && x>date + Dates.Hour(9) + Dates.Minute(30+ti*15) - Dates.Second(1), TW_mag2c1hz[!, :Time]) 
        ind = ind[1]
        local str = Dates.format(date + Dates.Hour(9) + Dates.Minute(30+ti*15), "HH:MM")
        for axi in [ax2]
            scatter!(axi, TW_mag2c1hz[ind, :Probe_Position_Y_MSO], TW_mag2c1hz[ind, :Probe_Position_Z_MSO], markersize=15, label=str, color=color, transparency=true)
        end
        if ti == 0 
            text!(ax2, TW_mag2c1hz[ind, :Probe_Position_Y_MSO]+0.15, TW_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.10, text=str, align = (:left, :top), space = :data, fontsize = 20, color=color, rotation=deg2rad(-90))
        elseif ti == 1 || ti == 3
            text!(ax2, TW_mag2c1hz[ind, :Probe_Position_Y_MSO]-0.15, TW_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.1, text=str, align = (:right, :bottom), space = :data, fontsize = 20, color=color, rotation=deg2rad(-90))
        else
            text!(ax2, TW_mag2c1hz[ind, :Probe_Position_Y_MSO]+0.15, TW_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.05, text=str, align = (:left, :top), space = :data, fontsize = 20, color=color, rotation=deg2rad(-90))
        end
        
        
    end
    text!(ax2, 0.8, 2, text="Tianwen-1", align = (:center, :bottom), space = :data, fontsize = 20, color=cn[1],rotation=deg2rad(-60))
    text!(ax2, 0.5, -3.3, text="MAVEN", align = (:center, :top), space = :data, fontsize = 20, color=cn[3], rotation=deg2rad(10))


    axs = [ax1, ax2]
    txt = ["A", "B", "C", "D", "E", "F"]
    for (i, axi) in enumerate(axs)
        text!(axi, 0.0, 1.0, text = txt[i], font = :bold, align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 20)
    end

    linkyaxes!(ax1, ax2)
    hideydecorations!(ax2, grid=false)
    colsize!(fig.layout, 2, Relative(4/10))
    fig
end

function orbitxyzBackup_pre_plot()
    size_inches = (9, 7)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5,16, 3, 8))
    cn = cgrad(:seaborn_deep, 10, categorical = true)

    ax1 = Axis(fig[1,1], aspect=DataAspect(), xlabel=rich(rich("x", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich(rich("z", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10), yticks=range(-10, 10), yautolimitmargin=(0.00f0, 0.05f0), xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticksvisible = true, backgroundcolor=:gray80)
    xlims!(ax1, -4.5, 1.5) 
    ylims!(ax1, -4, 4) 

    img = load("marsxz2.png")
    image!(ax1, (-1, 1), (-1, 1), rotr90(img))
    
    lines!(ax1, TW_mag2c1hz[!, :Probe_Position_X_MSO], TW_mag2c1hz[!, :Probe_Position_Z_MSO]; color=cn[1], overdraw = true)

    ind1 = findall(x->x>0 , MVN_mag2c1hz[!, :Probe_Position_Y_MSO])
    r = @. sqrt(MVN_mag2c1hz[!, :Probe_Position_X_MSO]^2 + MVN_mag2c1hz[!, :Probe_Position_Z_MSO]^2)
    ind2 = findall(x->x<1 , r)
    ind = intersect(ind1, ind2)
    tx = MVN_mag2c1hz[:, :Probe_Position_X_MSO]
    tz = MVN_mag2c1hz[:, :Probe_Position_Z_MSO]
    tx[ind] .= NaN64
    tz[ind] .= NaN64
    lines!(ax1, tx, tz; color=cn[3], overdraw = true)

    for ti in 0:2
        local ind = findall(x-> x<date + Dates.Hour(2) + Dates.Minute(0+ti*30) + Dates.Second(1)  && x>date + Dates.Hour(2) + Dates.Minute(0+ti*30) - Dates.Second(1), MVN_mag2c1hz[!, :Time]) 
        local str = Dates.format(date + Dates.Hour(2) + Dates.Minute(0+ti*30), "HH:MM")
        ind = ind[1]
        for axi in [ax1]
            scatter!(axi, MVN_mag2c1hz[ind, :Probe_Position_X_MSO], MVN_mag2c1hz[ind, :Probe_Position_Z_MSO], markersize=15, label=str, color=cn[ti+4], transparency=true)   
        end
        if ti == 0 
            text!(ax1, MVN_mag2c1hz[ind, :Probe_Position_X_MSO]-0.17, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.20, text=str, align = (:right, :top), space = :data, fontsize = 20, color=cn[ti+4], rotation=deg2rad(20))
        elseif ti == 1
            text!(ax1, MVN_mag2c1hz[ind, :Probe_Position_X_MSO]-0.17, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.20, text=str, align = (:right, :top), space = :data, fontsize = 20, color=cn[ti+4], rotation=deg2rad(20))
        else
            text!(ax1, MVN_mag2c1hz[ind, :Probe_Position_X_MSO]+0.08, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.05, text=str, align = (:left, :center), space = :data, fontsize = 20, color=cn[ti+4], rotation=deg2rad(20))
        end
        local ind = findall(x-> x<date + Dates.Hour(2) + Dates.Minute(0+ti*30) + Dates.Second(1)  && x>date + Dates.Hour(2) + Dates.Minute(0+ti*30) - Dates.Second(1), TW_mag2c1hz[!, :Time]) 
        ind = ind[1]
        local str = Dates.format(date + Dates.Hour(2) + Dates.Minute(0+ti*30), "HH:MM")
        for axi in [ax1]
            scatter!(axi, TW_mag2c1hz[ind, :Probe_Position_X_MSO], TW_mag2c1hz[ind, :Probe_Position_Z_MSO], markersize=15, label=str, color=cn[ti+4], transparency=true)
        end
        if ti == 0 
            text!(ax1, TW_mag2c1hz[ind, :Probe_Position_X_MSO]+0.15, TW_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.08, text=str, align = (:left, :top), space = :data, fontsize = 20, color=cn[ti+4], rotation=deg2rad(-80))
        elseif ti == 1
            text!(ax1, TW_mag2c1hz[ind, :Probe_Position_X_MSO]+0.15, TW_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.08, text=str, align = (:left, :top), space = :data, fontsize = 20, color=cn[ti+4], rotation=deg2rad(-80))
        else
            text!(ax1, TW_mag2c1hz[ind, :Probe_Position_X_MSO]+0.15, TW_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.08, text=str, align = (:left, :top), space = :data, fontsize = 20, color=cn[ti+4], rotation=deg2rad(-80))
        end
        
    end
    text!(ax1, -2.3, 2.8, text="Tianwen-1", align = (:center, :bottom), space = :data, fontsize = 20, color=cn[1],rotation=deg2rad(10))
    text!(ax1, 0.5, -2.9, text="MAVEN", align = (:center, :bottom), space = :data, fontsize = 20, color=cn[3], rotation=deg2rad(20))



    ax2 = Axis(fig[1,2], aspect=DataAspect(), xlabel=rich(rich("y", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich(rich("z", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10), yticks=range(-10, 10), yautolimitmargin=(0.00f0, 0.05f0), xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticksvisible = true, backgroundcolor=:gray80)
    xlims!(ax2, -2, 2) 
    ylims!(ax2, -4, 4) 

    img = load("marsyz2.png")
    image!(ax2, (-1, 1), (-1, 1), rotr90(img))
    
    ind1 = findall(x->x<0 , TW_mag2c1hz[!, :Probe_Position_X_MSO])
    r = @. sqrt(TW_mag2c1hz[!, :Probe_Position_Y_MSO]^2 + TW_mag2c1hz[!, :Probe_Position_Z_MSO]^2)
    ind2 = findall(x->x<1 , r)
    ind = intersect(ind1, ind2)
    tx = TW_mag2c1hz[:, :Probe_Position_Y_MSO]
    tz = TW_mag2c1hz[:, :Probe_Position_Z_MSO]
    tx[ind] .= NaN64
    tz[ind] .= NaN64
    lines!(ax2, tx, tz; color=cn[1], overdraw = true)

 
    tx = MVN_mag2c1hz[:, :Probe_Position_Y_MSO]
    tz = MVN_mag2c1hz[:, :Probe_Position_Z_MSO]
    lines!(ax2, tx, tz; color=cn[3], overdraw = true)

    for ti in 0:2
        local ind = findall(x-> x<date + Dates.Hour(2) + Dates.Minute(0+ti*30) + Dates.Second(1)  && x>date + Dates.Hour(2) + Dates.Minute(0+ti*30) - Dates.Second(1), MVN_mag2c1hz[!, :Time]) 
        local str = Dates.format(date + Dates.Hour(2) + Dates.Minute(0+ti*30), "HH:MM")
        ind = ind[1]
        for axi in [ax2]
            scatter!(axi, MVN_mag2c1hz[ind, :Probe_Position_Y_MSO], MVN_mag2c1hz[ind, :Probe_Position_Z_MSO], markersize=15, label=str, color=cn[ti+4], transparency=true)   
        end
        if ti == 0 
            text!(ax2, MVN_mag2c1hz[ind, :Probe_Position_Y_MSO]-0.00, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.1, text=str, align = (:left, :center), space = :data, fontsize = 20, color=cn[ti+4], rotation=deg2rad(-85))
        elseif ti == 1
            text!(ax2, MVN_mag2c1hz[ind, :Probe_Position_Y_MSO]-0.00, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.1, text=str, align = (:left, :center), space = :data, fontsize = 20, color=cn[ti+4], rotation=deg2rad(-85))
        else
            text!(ax2, MVN_mag2c1hz[ind, :Probe_Position_Y_MSO]+0.00, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.1, text=str, align = (:left, :center), space = :data, fontsize = 20, color=cn[ti+4], rotation=deg2rad(-85))
        end
        local ind = findall(x-> x<date + Dates.Hour(2) + Dates.Minute(0+ti*30) + Dates.Second(1)  && x>date + Dates.Hour(2) + Dates.Minute(0+ti*30) - Dates.Second(1), TW_mag2c1hz[!, :Time]) 
        ind = ind[1]
        local str = Dates.format(date + Dates.Hour(2) + Dates.Minute(0+ti*30), "HH:MM")
        for axi in [ax2]
            scatter!(axi, TW_mag2c1hz[ind, :Probe_Position_Y_MSO], TW_mag2c1hz[ind, :Probe_Position_Z_MSO], markersize=15, label=str, color=cn[ti+4], transparency=true)
        end
        if ti == 0 
            text!(ax2, TW_mag2c1hz[ind, :Probe_Position_Y_MSO]-0.25, TW_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.10, text=str, align = (:right, :bottom), space = :data, fontsize = 20, color=cn[ti+4], rotation=deg2rad(-90))
        elseif ti == 1
            text!(ax2, TW_mag2c1hz[ind, :Probe_Position_Y_MSO]-0.15, TW_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.15, text=str, align = (:right, :bottom), space = :data, fontsize = 20, color=cn[ti+4], rotation=deg2rad(-90))
        else
            text!(ax2, TW_mag2c1hz[ind, :Probe_Position_Y_MSO]-0.15, TW_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.1, text=str, align = (:right, :bottom), space = :data, fontsize = 20, color=cn[ti+4], rotation=deg2rad(-90))
        end
        
        
    end
    text!(ax2, 0.8, 2, text="Tianwen-1", align = (:center, :bottom), space = :data, fontsize = 20, color=cn[1],rotation=deg2rad(-60))
    text!(ax2, 0.5, -3.3, text="MAVEN", align = (:center, :top), space = :data, fontsize = 20, color=cn[3], rotation=deg2rad(10))


    axs = [ax1, ax2]
    txt = ["A", "B", "C", "D", "E", "F"]
    for (i, axi) in enumerate(axs)
        text!(axi, 0.0, 1.0, text = txt[i], align = (:left, :top), offset = (4, -2), space = :relative, font = :bold, fontsize = 20)
    end

    linkyaxes!(ax1, ax2)
    hideydecorations!(ax2, grid=false)
    colsize!(fig.layout, 2, Relative(4/10))
    fig
end

function orbitxyz_pre_plot()
    size_inches = (9, 6.5)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5,5, 3, 8))
    cn = cgrad(:seaborn_deep, 10, categorical = true)
    lcn = Makie.wong_colors()

    ax1 = Axis(fig[1,1], aspect=DataAspect(), xlabel=rich(rich("x", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich(rich("z", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10), yticks=range(-10, 10), yautolimitmargin=(0.00f0, 0.05f0), xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticksvisible = true, backgroundcolor=:gray80)
    xlims!(ax1, -4.5, 1.5) 
    ylims!(ax1, -4, 4) 

    img = load("marsxz2.png")
    image!(ax1, (-1, 1), (-1, 1), rotr90(img))
    
    lines!(ax1, TW_mag2c1hz[!, :Probe_Position_X_MSO], TW_mag2c1hz[!, :Probe_Position_Z_MSO]; color=lcn[1], overdraw = true)

    ind1 = findall(x->x>0 , MVN_mag2c1hz[!, :Probe_Position_Y_MSO])
    r = @. sqrt(MVN_mag2c1hz[!, :Probe_Position_X_MSO]^2 + MVN_mag2c1hz[!, :Probe_Position_Z_MSO]^2)
    ind2 = findall(x->x<1 , r)
    ind = intersect(ind1, ind2)
    tx = MVN_mag2c1hz[:, :Probe_Position_X_MSO]
    tz = MVN_mag2c1hz[:, :Probe_Position_Z_MSO]
    tx[ind] .= NaN64
    tz[ind] .= NaN64
    lines!(ax1, tx, tz; color=lcn[3], overdraw = true)

    for ti in 0:2
        local ind = findall(x-> x<date + Dates.Hour(2) + Dates.Minute(0+ti*30) + Dates.Second(1)  && x>date + Dates.Hour(2) + Dates.Minute(0+ti*30) - Dates.Second(1), MVN_mag2c1hz[!, :Time]) 
        local str = Dates.format(date + Dates.Hour(2) + Dates.Minute(0+ti*30), "HH:MM")
        ind = ind[1]
        for axi in [ax1]
            scatter!(axi, MVN_mag2c1hz[ind, :Probe_Position_X_MSO], MVN_mag2c1hz[ind, :Probe_Position_Z_MSO], markersize=15, label=str, color=cn[ti+4], transparency=true)   
        end
        if ti == 0 
            text!(ax1, MVN_mag2c1hz[ind, :Probe_Position_X_MSO]-0.2, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.20, text=str, align = (:right, :top), space = :data, fontsize = 15, color=cn[ti+4], rotation=deg2rad(20))
        elseif ti == 1
            text!(ax1, MVN_mag2c1hz[ind, :Probe_Position_X_MSO]-0.2, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.20, text=str, align = (:right, :top), space = :data, fontsize = 15, color=cn[ti+4], rotation=deg2rad(20))
        else
            text!(ax1, MVN_mag2c1hz[ind, :Probe_Position_X_MSO]+0.11, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.05, text=str, align = (:left, :center), space = :data, fontsize = 15, color=cn[ti+4], rotation=deg2rad(20))
        end
        local ind = findall(x-> x<date + Dates.Hour(2) + Dates.Minute(0+ti*30) + Dates.Second(1)  && x>date + Dates.Hour(2) + Dates.Minute(0+ti*30) - Dates.Second(1), TW_mag2c1hz[!, :Time]) 
        ind = ind[1]
        local str = Dates.format(date + Dates.Hour(2) + Dates.Minute(0+ti*30), "HH:MM")
        for axi in [ax1]
            scatter!(axi, TW_mag2c1hz[ind, :Probe_Position_X_MSO], TW_mag2c1hz[ind, :Probe_Position_Z_MSO], markersize=15, label=str, color=cn[ti+4], transparency=true)
        end
        if ti == 0 
            text!(ax1, TW_mag2c1hz[ind, :Probe_Position_X_MSO]+0.15, TW_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.1, text=str, align = (:left, :top), space = :data, fontsize = 15, color=cn[ti+4], rotation=deg2rad(-80))
        elseif ti == 1
            text!(ax1, TW_mag2c1hz[ind, :Probe_Position_X_MSO]+0.15, TW_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.1, text=str, align = (:left, :top), space = :data, fontsize = 15, color=cn[ti+4], rotation=deg2rad(-80))
        else
            text!(ax1, TW_mag2c1hz[ind, :Probe_Position_X_MSO]+0.15, TW_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.1, text=str, align = (:left, :top), space = :data, fontsize = 15, color=cn[ti+4], rotation=deg2rad(-80))
        end
        
    end
    text!(ax1, -2.3, 2.8, text="Tianwen-1", align = (:center, :bottom), space = :data, fontsize = 15, color=lcn[1],rotation=deg2rad(10))
    text!(ax1, 0.5, -2.9, text="MAVEN", align = (:center, :bottom), space = :data, fontsize = 15, color=lcn[3], rotation=deg2rad(20))



    ax2 = Axis(fig[1,3], aspect=DataAspect(), xlabel=rich(rich("y", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich(rich("z", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10), yticks=range(-10, 10), yautolimitmargin=(0.00f0, 0.05f0), xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticksvisible = true, backgroundcolor=:gray80)
    xlims!(ax2, -2, 2) 
    ylims!(ax2, -4, 4) 

    img = load("marsyz2.png")
    image!(ax2, (-1, 1), (-1, 1), rotr90(img))
    
    ind1 = findall(x->x<0 , TW_mag2c1hz[!, :Probe_Position_X_MSO])
    r = @. sqrt(TW_mag2c1hz[!, :Probe_Position_Y_MSO]^2 + TW_mag2c1hz[!, :Probe_Position_Z_MSO]^2)
    ind2 = findall(x->x<1 , r)
    ind = intersect(ind1, ind2)
    tx = TW_mag2c1hz[:, :Probe_Position_Y_MSO]
    tz = TW_mag2c1hz[:, :Probe_Position_Z_MSO]
    tx[ind] .= NaN64
    tz[ind] .= NaN64
    lines!(ax2, tx, tz; color=lcn[1], overdraw = true)

 
    tx = MVN_mag2c1hz[:, :Probe_Position_Y_MSO]
    tz = MVN_mag2c1hz[:, :Probe_Position_Z_MSO]
    lines!(ax2, tx, tz; color=lcn[3], overdraw = true)

    for ti in 0:2
        local ind = findall(x-> x<date + Dates.Hour(2) + Dates.Minute(0+ti*30) + Dates.Second(1)  && x>date + Dates.Hour(2) + Dates.Minute(0+ti*30) - Dates.Second(1), MVN_mag2c1hz[!, :Time]) 
        local str = Dates.format(date + Dates.Hour(2) + Dates.Minute(0+ti*30), "HH:MM")
        ind = ind[1]
        for axi in [ax2]
            scatter!(axi, MVN_mag2c1hz[ind, :Probe_Position_Y_MSO], MVN_mag2c1hz[ind, :Probe_Position_Z_MSO], markersize=15, label=str, color=cn[ti+4], transparency=true)   
        end
        if ti == 0 
            text!(ax2, MVN_mag2c1hz[ind, :Probe_Position_Y_MSO]-0.00, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.15, text=str, align = (:left, :center), space = :data, fontsize = 15, color=cn[ti+4], rotation=deg2rad(-85))
        elseif ti == 1
            text!(ax2, MVN_mag2c1hz[ind, :Probe_Position_Y_MSO]-0.00, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.15, text=str, align = (:left, :center), space = :data, fontsize = 15, color=cn[ti+4], rotation=deg2rad(-85))
        else
            text!(ax2, MVN_mag2c1hz[ind, :Probe_Position_Y_MSO]+0.00, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.15, text=str, align = (:left, :center), space = :data, fontsize = 15, color=cn[ti+4], rotation=deg2rad(-85))
        end
        local ind = findall(x-> x<date + Dates.Hour(2) + Dates.Minute(0+ti*30) + Dates.Second(1)  && x>date + Dates.Hour(2) + Dates.Minute(0+ti*30) - Dates.Second(1), TW_mag2c1hz[!, :Time]) 
        ind = ind[1]
        local str = Dates.format(date + Dates.Hour(2) + Dates.Minute(0+ti*30), "HH:MM")
        for axi in [ax2]
            scatter!(axi, TW_mag2c1hz[ind, :Probe_Position_Y_MSO], TW_mag2c1hz[ind, :Probe_Position_Z_MSO], markersize=15, label=str, color=cn[ti+4], transparency=true)
        end
        if ti == 0 
            text!(ax2, TW_mag2c1hz[ind, :Probe_Position_Y_MSO]-0.25, TW_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.13, text=str, align = (:right, :bottom), space = :data, fontsize = 15, color=cn[ti+4], rotation=deg2rad(-90))
        elseif ti == 1
            text!(ax2, TW_mag2c1hz[ind, :Probe_Position_Y_MSO]-0.15, TW_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.18, text=str, align = (:right, :bottom), space = :data, fontsize = 15, color=cn[ti+4], rotation=deg2rad(-90))
        else
            text!(ax2, TW_mag2c1hz[ind, :Probe_Position_Y_MSO]-0.15, TW_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.13, text=str, align = (:right, :bottom), space = :data, fontsize = 15, color=cn[ti+4], rotation=deg2rad(-90))
        end
        
        
    end
    text!(ax2, 0.8, 2, text="Tianwen-1", align = (:center, :bottom), space = :data, fontsize = 15, color=lcn[1],rotation=deg2rad(-60))
    text!(ax2, 0.5, -3.3, text="MAVEN", align = (:center, :top), space = :data, fontsize = 15, color=lcn[3], rotation=deg2rad(10))

    txt = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"]
    Label(fig[1, 0], txt[1],font = :bold, valign=:top, tellheight=false, fontsize = 20)
    Label(fig[1, 2], txt[2],font = :bold, valign=:top, tellheight=false, fontsize = 20)

    linkyaxes!(ax1, ax2)
    hideydecorations!(ax2, grid=false)
    colsize!(fig.layout, 1, Relative(5.7/10))
    colsize!(fig.layout, 3, Relative(3.8/10))
    fig
end

function schematic_plot()
    size_inches = (9, 9)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5,16, 3, 8))
    cn = cgrad(:seaborn_deep, 10, categorical = true)

    ax1 = Axis(fig[1,1], aspect=DataAspect(), xlabel=rich(rich("x", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich(rich("z", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10), yticks=range(-10, 10), yautolimitmargin=(0.00f0, 0.05f0), xgridvisible = false, ygridvisible = false)
    xlims!(ax1, -5, 3) 
    ylims!(ax1, -4, 4) 

    x = range(-10, 10, length=801)
    y = range(-10, 10, length=801)
    
    z = [exp(-sqrt(tx^2+ty^2)^1.7/2.0) for tx in x, ty in y]
    ind = findall(z .> exp(-1.0/2.0))
    z[ind] .= exp(-1.0/2.0)
	z =  round.(Int, z./maximum(z)*256)
	ind = findall(x->x<1, z)
	z[ind] .= 1


    #contourf!(ax1, x, y, z, colormap=Reverse(:bilbao), levels = range(0.1, 1.0, length = 50))
    #surface!(ax1, x, y, zeros(size(z)), colormap=Reverse(:bilbao), color=z,shading = NoShading)
	#heatmap!(ax1, x, y, z, colormap=Reverse(:bilbao))
	
	imcn = reverse(cgrad(:bilbao, 256, categorical = true))
	imz = imcn[z]
    image!(ax1, (-10, 10), (-10, 10), imz)

    img = load("marsxz1.png")
    image!(ax1, (-1, 1), (-1, 1), rotr90(img))

    θ = range(0, 2π, length=3600)
    lines!(ax1, cos.(θ), sin.(θ), color=imcn[z[441,401]], linewidth=1.0, overdraw = true)

    
    lines!(ax1, xshock, ryzshock; color=cn[8], linewidth=4, overdraw = true)
    lines!(ax1, xmp, ryzmp; color=cn[8], linewidth=4,  overdraw = true)

    for ri in 2:4
        @show ri
        lx, lz, va = stream(3.0, 0.75ri)
        lines!(ax1, lx, lz; color=cn[9], overdraw = true, linewidth=2)
        lines!(ax1, lx, -lz; color=cn[9], overdraw = true, linewidth=2)
        arrows!(ax1, [lx[end]], [lz[end]], [-1.0], [0.0], arrowsize = 15, lengthscale = 0.5, color=cn[9])
        arrows!(ax1, [lx[end]], [-lz[end]], [-1.0], [0.0], arrowsize = 15, lengthscale = 0.5, color=cn[9])

    end

    for bi in -8:10
        ablines!(bi, tand(35), color=cn[7], linewidth=2, overdraw = true)
        arrows!(ax1, [3], [3*tand(35)+bi], [-cosd(35)], [-sind(35)], arrowsize = 15, lengthscale = 0.5, color=cn[7])
    end

    
    lines!(ax1, TW_mag2c1hz[!, :Probe_Position_X_MSO], TW_mag2c1hz[!, :Probe_Position_Z_MSO]; color=cn[1], overdraw = true)

    ind1 = findall(x->x>0 , MVN_mag2c1hz[!, :Probe_Position_Y_MSO])
    r = @. sqrt(MVN_mag2c1hz[!, :Probe_Position_X_MSO]^2 + MVN_mag2c1hz[!, :Probe_Position_Z_MSO]^2)
    ind2 = findall(x->x<1 , r)
    ind = intersect(ind1, ind2)
    tx = MVN_mag2c1hz[:, :Probe_Position_X_MSO]
    tz = MVN_mag2c1hz[:, :Probe_Position_Z_MSO]
    tx[ind] .= NaN64
    tz[ind] .= NaN64
    lines!(ax1, tx, tz; color=cn[3], overdraw = true)


    for ti in 0:2
        local ind = findall(x-> x<date + Dates.Hour(9) + Dates.Minute(30+ti*15) + Dates.Second(1)  && x>date + Dates.Hour(9) + Dates.Minute(30+ti*15) - Dates.Second(1), MVN_mag2c1hz[!, :Time]) 
        local str = Dates.format(date + Dates.Hour(9) + Dates.Minute(30+ti*15), "HH:MM")
        ind = ind[1]
        for axi in [ax1]
            scatter!(axi, MVN_mag2c1hz[ind, :Probe_Position_X_MSO], MVN_mag2c1hz[ind, :Probe_Position_Z_MSO], markersize=15, label=str, color=cn[ti+4], transparency=true)   
        end
        if ti == 0 
            text!(ax1, MVN_mag2c1hz[ind, :Probe_Position_X_MSO]-0.15, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.08, text=str, font = :bold, align = (:right, :top), space = :data, fontsize = 20, color=cn[ti+4], rotation=deg2rad(20))
        elseif ti == 1
            text!(ax1, MVN_mag2c1hz[ind, :Probe_Position_X_MSO]-0.15, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.03, text=str, font = :bold, align = (:right, :top), space = :data, fontsize = 20, color=cn[ti+4], rotation=deg2rad(20))
        else
            text!(ax1, MVN_mag2c1hz[ind, :Probe_Position_X_MSO]+0.15, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]+0.05, text=str, font = :bold, align = (:left, :center), space = :data, fontsize = 20, color=cn[ti+4], rotation=deg2rad(20))
        end
        local ind = findall(x-> x<date + Dates.Hour(9) + Dates.Minute(30+ti*15) + Dates.Second(1)  && x>date + Dates.Hour(9) + Dates.Minute(30+ti*15) - Dates.Second(1), TW_mag2c1hz[!, :Time]) 
        ind = ind[1]
        local str = Dates.format(date + Dates.Hour(9) + Dates.Minute(30+ti*15), "HH:MM")
        for axi in [ax1]
            scatter!(axi, TW_mag2c1hz[ind, :Probe_Position_X_MSO], TW_mag2c1hz[ind, :Probe_Position_Z_MSO], markersize=15, label=str, color=cn[ti+4], transparency=true)
        end
        text!(ax1, TW_mag2c1hz[ind, :Probe_Position_X_MSO]+0.15, TW_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.08, text=str, font = :bold, align = (:left, :top), space = :data, fontsize = 20, color=cn[ti+4], rotation=deg2rad(-80))
        
    end


    text!(ax1, 2.8, 2.95, text=rich("B", font = :italic), font = :bold, align = (:center, :top), space = :data, fontsize = 20, color=cn[7])
    text!(ax1, 2.8, 3.0, text=rich("V", font = :italic), font = :bold, align = (:center, :bottom), space = :data, fontsize = 20, color=cn[9])
    text!(ax1, -2.3, 2.8, text="Tianwen-1", font = :bold, align = (:center, :bottom), space = :data, fontsize = 20, color=cn[1],rotation=deg2rad(10))
    text!(ax1, 0.5, -2.7, text="MAVEN", font = :bold, align = (:center, :bottom), space = :data, fontsize = 20, color=cn[3], rotation=deg2rad(20))
    fig
end

function test2_schematic_plot()
    size_inches = (9, 9)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5,16, 3, 8))
    cn = cgrad(:seaborn_deep, 10, categorical = true)

    ax1 = Axis(fig[1,1], aspect=DataAspect(), xlabel=rich(rich("x", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich(rich("z", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10), yticks=range(-10, 10), yautolimitmargin=(0.00f0, 0.05f0), xgridvisible = false, ygridvisible = false)
    xlims!(ax1, -5, 3) 
    ylims!(ax1, -4, 4) 

    
    x = range(-10, 10, length=801)
    y = range(-10, 10, length=801)
    
    z = [exp(-sqrt(tx^2+ty^2)^1.7/2.0) for tx in x, ty in y]
    ind = findall(z .> exp(-1.0/2.0))
    z[ind] .= exp(-1.0/2.0)
	z =  round.(Int, z./maximum(z)*256)
	ind = findall(x->x<1, z)
	z[ind] .= 1


    #contourf!(ax1, x, y, z, colormap=Reverse(:bilbao), levels = range(0.1, 1.0, length = 50))
    #surface!(ax1, x, y, zeros(size(z)), colormap=Reverse(:bilbao), color=z,shading = NoShading)
	#heatmap!(ax1, x, y, z, colormap=Reverse(:bilbao))
	
	imcn = reverse(cgrad(:bilbao, 256, categorical = true))
	imz = imcn[z]
    image!(ax1, (-10, 10), (-10, 10), imz)

    img = load("marsxz1.png")
    image!(ax1, (-1, 1), (-1, 1), rotr90(img))

    θ = range(0, 2π, length=3600)
    lines!(ax1, cos.(θ), sin.(θ), color=imcn[z[441,401]], linewidth=1.0, overdraw = true)

    
    
    
    #=
    mars_img = load("mars-map.jpg")
    ms1 = mesh!(ax1, SphereTess(; o=Point3f(0, 0, 0), r=1.00); color=mars_img)
    Makie.rotate!(ms1, Vec3f(1, 1, 0), 0.2π)
    #Makie.rotate!(ms1, Vec3f(0, 1, 0), π/2)
    =#

	#=

    lines!(ax1, xshock, ryzshock; color=cn[8], linewidth=4, overdraw = true)
    lines!(ax1, xmp, ryzmp; color=cn[8], linewidth=4,  overdraw = true)

    
    for ri in 2:4
        @show ri
        lx, lz, va = stream(3.0, 0.75ri)
        lines!(ax1, lx, lz; color=cn[9], overdraw = true, linewidth=2)
        lines!(ax1, lx, -lz; color=cn[9], overdraw = true, linewidth=2)
        arrows!(ax1, [lx[end]], [lz[end]], [-1.0], [0.0], arrowsize = 15, lengthscale = 0.5, color=cn[9])
        arrows!(ax1, [lx[end]], [-lz[end]], [-1.0], [0.0], arrowsize = 15, lengthscale = 0.5, color=cn[9])

    end
	=#
    fig
end

function MVN_TW_Plasma_Overview_plot()
    size_inches = (9, 10)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5,20, 9, 8))
    cn = cgrad(:seaborn_dark, 10, categorical = true)
    lcn = Makie.wong_colors()
    xd = date .+ Dates.Hour(9) .+ Dates.Minute(0) .+ Dates.Minute(10) .* range(0, 9)
    
    xtk = datetime2julian.(xd)
    xtk = Float32.(xtk .- MVN_kp[1, :JulUT])

    ax11 = Axis(fig[1, 1], xticks=xtk, yminorticks=IntervalsBetween(5), yminorticksvisible = true, backgroundcolor=:gray80)
    ax11.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(ax11, minimum(xtk), maximum(xtk))
    ylims!(ax11, -350, -150)
    #title = rich("Measurements of ", rich("MAVEN ", color=lcn[1]), "and ", rich("Tianwen-1 ", color=lcn[2]), "in the MSO Coordinate system"), 
    text!(ax11, 0.0, 0.6, text = rich(rich("MAVEN", color=lcn[3]), "\n", rich("Tianwen-1", color=lcn[1])), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 20)

    lines!(ax11, MVN_kp[!, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[!, 43], label = rich(rich("V", font = :italic), subscript("x")), color=lcn[3])
    lines!(ax11, minpajulianmod1.-MVN_kp[1, :JulUT], VH[:, 1], label = rich(rich("V", font = :italic), subscript("x")), color=lcn[1])

    ax12 = Axis(fig[2, 1],  ylabel = rich("Flow Velocity (km⋅s", superscript("-1"), ")"), xticks=xtk, yminorticks=IntervalsBetween(5), yminorticksvisible = true, backgroundcolor=:gray80, ylabelpadding=10.0)
    ax12.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(ax12, minimum(xtk), maximum(xtk))
    ylims!(ax12, -100, 150)

    lines!(ax12, MVN_kp[!, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[!, 45], label = rich(rich("V", font = :italic), subscript("y")), color=lcn[3])
    #lines!(ax12, minpajulianmod1.-MVN_kp[1, :JulUT], VH[:, 2], label = rich(rich("V", font = :italic), subscript("y")), color=lcn[1])

    ax13 = Axis(fig[3, 1],  xticks=xtk, yminorticks=IntervalsBetween(5), yminorticksvisible = true, backgroundcolor=:gray80)
    ax13.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(ax13, minimum(xtk), maximum(xtk))
    ylims!(ax13, -150, 100)

    lines!(ax13, MVN_kp[!, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[!, 47], label = rich(rich("V", font = :italic), subscript("z")), color=lcn[3])
    lines!(ax13, minpajulianmod1.-MVN_kp[1, :JulUT], VH[:, 3], label = rich(rich("V", font = :italic), subscript("z")), color=lcn[1])

 
    
    ax2=Axis(fig[4,1], xlabel = "Universal Time", ylabel = rich("Density (cm", superscript("-3"), ")"), xticks=xtk, yminorticks=IntervalsBetween(5), yminorticksvisible = true, backgroundcolor=:gray80, ylabelpadding=10.0)
    ax2.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(ax2, minimum(xtk), maximum(xtk)) 
    ylims!(ax2, 0, 10)
    #ax2.yreversed = true
    lines!(ax2, MVN_kp[!, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[!, 41], label=rich(rich("N", font = :italic), subscript("p")), color=lcn[3])
    lines!(ax2, minpajulianmod1.-MVN_kp[1, :JulUT], dropdims(nH, dims=(2,)), label=rich(rich("N", font = :italic), subscript("p")), color=lcn[1])

    ax3=Axis(fig[5,1], xlabel = "Universal Time", ylabel = "Temperature (eV)", xticks=xtk, yminorticks=IntervalsBetween(5), yminorticksvisible = true, backgroundcolor=:gray80, ylabelpadding=10.0)
    ax3.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(ax3, minimum(xtk), maximum(xtk)) 
    ylims!(ax3, 50, 250)
    #ax3.yreversed = true
    lines!(ax3, MVN_kp[!, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[!, 49], label=rich(rich("T", font = :italic), subscript("p")), color=lcn[3])
    #lines!(ax3, minpajulianmod1.-MVN_kp[1, :JulUT], dropdims(TH, dims=(2,)), label=rich(rich("T", font = :italic), subscript("p")), color=lcn[1])
    

    ax = [ax11, ax12, ax13, ax2, ax3]
    txt = ["A MAVEN", "B MAVEN", "C Tianwen-1", "D Tianwen-1"]
    txt = ["A", "B", "C", "D", "E", "F"]
    for (i, axi) in enumerate(ax)
        text!(axi, 0.0, 1.0, text = txt[i], font = :bold, align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 20)
    end

    linkxaxes!(ax11, ax12, ax13, ax2, ax3)
    hidexdecorations!(ax11, grid=false)
    hidexdecorations!(ax12, grid=false)
    hidexdecorations!(ax13, grid=false)
    hidexdecorations!(ax2, grid=false)
    #rowsize!(fig.layout, 4, Relative(1.2/5))
    #rowsize!(fig.layout, 5, Relative(1.2/5))
    
    fig
end

function MVN_TW_Mag_Overview_plot()
    size_inches = (9, 11)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5,20, 9, 8))
    cn = cgrad(:seaborn_dark, 10, categorical = true)
    lcn = Makie.wong_colors()
    xd = date .+ Dates.Hour(9) .+ Dates.Minute(0) .+ Dates.Minute(10) .* range(0, 9)
    
    xtk = datetime2julian.(xd)
    xtk = Float32.(xtk .- MVN_kp[1, :JulUT])

    ax11 = Axis(fig[1, 1],  xticks=xtk, yticks=range(-16, 16, step=4), yminorticks=IntervalsBetween(4), yminorticksvisible = true, backgroundcolor=:gray80)
    #title = rich(rich("MAVEN ", color=lcn[1]), rich("Tianwen-1", color=lcn[2]), " MSO Components"),
    ax11.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(ax11, minimum(xtk), maximum(xtk))
    ylims!(ax11, -12, 8)
    text!(ax11, 0.0, 0.7, text = rich(rich("Tianwen-1", color=lcn[1]), "\n",  rich("MAVEN", color=lcn[3])), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 20)
    

    lines!(ax11, MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_mag2c1hz[!, :X_MSO], label = rich(rich("B", font = :italic), subscript("x")), color=lcn[3])
    lines!(ax11, TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_mag2c1hz[!, :X_MSO], label = rich(rich("B", font = :italic), subscript("x")), color=lcn[1])

    ax12 = Axis(fig[2, 1],  ylabel = "Magnetic Field (nT)", xticks=xtk, yticks=range(-16, 16, step=4), yminorticks=IntervalsBetween(4), yminorticksvisible = true, backgroundcolor=:gray80, ylabelpadding=10.0)
    ax12.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(ax12, minimum(xtk), maximum(xtk))
    ylims!(ax12, -8, 8)

    lines!(ax12, MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_mag2c1hz[!, :Y_MSO], label = rich(rich("B", font = :italic), subscript("y")), color=lcn[3])
    lines!(ax12, TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_mag2c1hz[!, :Y_MSO], label = rich(rich("B", font = :italic), subscript("y")), color=lcn[1])

    ax13 = Axis(fig[3, 1],  xticks=xtk, yticks=range(-16, 16, step=4), yminorticks=IntervalsBetween(4), yminorticksvisible = true, backgroundcolor=:gray80)
    ax13.xticks = (xtk, Dates.format.(xd, "HH:MM:SS"))
    xlims!(ax13, minimum(xtk), maximum(xtk))
    ylims!(ax13, -12, 4)

    lines!(ax13, MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_mag2c1hz[!, :Z_MSO], label = rich(rich("B", font = :italic), subscript("z")), color=lcn[3])
    lines!(ax13, TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_mag2c1hz[!, :Z_MSO], label = rich(rich("B", font = :italic), subscript("z")), color=lcn[1])

 
    
    ax2=Axis(fig[4,1], xlabel = "Universal Time", xticks=xtk, yscale=log2)
    ax2.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(ax2, minimum(xtk), maximum(xtk)) 
    ylims!(ax2, 2, 64)
    ax2.yreversed = true
    hp1 = heatmap!(ax2, MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_period1hz, MVN_Bpower1hz', colorscale=log10, colormap=:gist_earth, colorrange=(3e-2, 3e2))
    lines!(ax2, MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_Tcp, color = :black, linestyle=:dash)
    text!(ax2, 0.0, 0.45, text = rich(rich("MAVEN", color=:black)), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 20)

    ax3=Axis(fig[5,1], xlabel = "Universal Time", ylabel = repeat(' ', 27) * "Wave Period (s)", xticks=xtk, yscale=log2, ylabelpadding=10.0)
    ax3.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(ax3, minimum(xtk), maximum(xtk)) 
    ylims!(ax3, 2, 64)
    ax3.yreversed = true
    hp2 = heatmap!(ax3, TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_period1hz, TW_Bpower1hz', colorscale=log10, colormap=:gist_earth, colorrange=(3e-2, 3e2))
    lines!(ax3, TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_Tcp, color = :black, linestyle=:dash)
    text!(ax3, 0.0, 0.45, text = rich( rich("Tianwen-1", color=:black)), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 20)

    Colorbar(fig[6, 1], hp1, label = rich("Wavelet Power"), vertical=false, flipaxis=false)



    ax = [ax11, ax12, ax13, ax2, ax3]
    txt = ["A MAVEN", "B MAVEN", "C Tianwen-1", "D Tianwen-1"]
    txt = ["A", "B", "C", "D", "E", "F"]
    for (i, axi) in enumerate(ax)
        text!(axi, 0.0, 1.0, text = txt[i], font = :bold, align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 20)
    end

    linkxaxes!(ax11, ax12, ax13, ax2, ax3)
    hidexdecorations!(ax11, grid=false)
    hidexdecorations!(ax12, grid=false)
    hidexdecorations!(ax13, grid=false)
    hidexdecorations!(ax2, grid=false)
    #rowsize!(fig.layout, 5, Relative(1/6))
    #rowsize!(fig.layout, 6, Relative(1/6))
    #rowsize!(fig.layout, 7, Relative(1.2/6))
    
    fig
end

function MVN_TW_MagBackup_Pre_Overview_plot()
    size_inches = (9, 11)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5,30, 3, 8))
    cn = cgrad(:seaborn_dark, 10, categorical = true)
    lcn = Makie.wong_colors()
    xd = date .+ Dates.Hour(1) .+ Dates.Minute(30) .+ Dates.Minute(15) .* range(0, 8)
    
    xtk = datetime2julian.(xd)
    xtk = Float32.(xtk .- MVN_kp[1, :JulUT])

    ax11 = Axis(fig[1, 1],  xticks=xtk, yticks=range(-15, 15, step=5), yminorticks=IntervalsBetween(5), yminorticksvisible = true, backgroundcolor=:gray80)
    ax11.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(ax11, minimum(xtk), maximum(xtk))
    ylims!(ax11, -10, 10)
    text!(ax11, 0.0, 0.75, text = rich(rich("Tianwen-1", color=lcn[1]), "\n", rich("MAVEN", color=lcn[3])), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 20)

    lines!(ax11, MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_mag2c1hz[!, :X_MSO], label = rich(rich("B", font = :italic), subscript("x")), color=lcn[3])
    lines!(ax11, TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_mag2c1hz[!, :X_MSO], label = rich(rich("B", font = :italic), subscript("x")), color=lcn[1])

    ax12 = Axis(fig[2, 1],  ylabel = "Magnetic Field (nT)", xticks=xtk, yticks=range(-15, 15, step=5), yminorticks=IntervalsBetween(5), yminorticksvisible = true, backgroundcolor=:gray80, ylabelpadding=10.0)
    ax12.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(ax12, minimum(xtk), maximum(xtk))
    ylims!(ax12, -10, 10)

    lines!(ax12, MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_mag2c1hz[!, :Y_MSO], label = rich(rich("B", font = :italic), subscript("y")), color=lcn[3])
    lines!(ax12, TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_mag2c1hz[!, :Y_MSO], label = rich(rich("B", font = :italic), subscript("y")), color=lcn[1])

    ax13 = Axis(fig[3, 1],  xticks=xtk, yticks=range(-15, 15, step=5), yminorticks=IntervalsBetween(5), yminorticksvisible = true, backgroundcolor=:gray80)
    ax13.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(ax13, minimum(xtk), maximum(xtk))
    ylims!(ax13, -15, 5)

    lines!(ax13, MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_mag2c1hz[!, :Z_MSO], label = rich(rich("B", font = :italic), subscript("z")), color=lcn[3])
    lines!(ax13, TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_mag2c1hz[!, :Z_MSO], label = rich(rich("B", font = :italic), subscript("z")), color=lcn[1])

 
    
    ax2=Axis(fig[4,1], xlabel = "Universal Time", xticks=xtk, yscale=log2)
    ax2.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(ax2, minimum(xtk), maximum(xtk)) 
    ylims!(ax2, 2, 64)
    ax2.yreversed = true
    hp1 = heatmap!(ax2, MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_period1hz, MVN_Bpower1hz', colorscale=log10, colormap=:gist_earth, colorrange=(3e-2, 3e2))
    lines!(ax2, MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_Tcp, color = :black, linestyle=:dash)
    text!(ax2, 0.0, 0.45, text = rich(rich("MAVEN", color=:black)), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 20)

    ax3=Axis(fig[5,1], xlabel = "Universal Time", ylabel = repeat(' ', 27) * "Wave Period (s)", xticks=xtk, yscale=log2, ylabelpadding=10.0)
    ax3.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(ax3, minimum(xtk), maximum(xtk)) 
    ylims!(ax3, 2, 64)
    ax3.yreversed = true
    hp2 = heatmap!(ax3, TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_period1hz, TW_Bpower1hz', colorscale=log10, colormap=:gist_earth, colorrange=(3e-2, 3e2))
    lines!(ax3, TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_Tcp, color = :black, linestyle=:dash)
    text!(ax3, 0.0, 0.45, text = rich( rich("Tianwen-1", color=:black)), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 20)

    Colorbar(fig[6, 1], hp1, label = rich("Wavelet Power"), vertical=false, flipaxis=false)



    ax = [ax11, ax12, ax13, ax2, ax3]
    txt = ["A MAVEN", "B MAVEN", "C Tianwen-1", "D Tianwen-1"]
    txt = ["A", "B", "C", "D", "E", "F"]
    for (i, axi) in enumerate(ax)
        text!(axi, 0.0, 1.0, text = txt[i], font = :bold, align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 20)
    end

    linkxaxes!(ax11, ax12, ax13, ax2, ax3)
    hidexdecorations!(ax11, grid=false)
    hidexdecorations!(ax12, grid=false)
    hidexdecorations!(ax13, grid=false)
    hidexdecorations!(ax2, grid=false)
    #rowsize!(fig.layout, 5, Relative(1/6))
    #rowsize!(fig.layout, 6, Relative(1/6))
    #rowsize!(fig.layout, 7, Relative(1.2/6))
    
    fig
end

function MVN_TW_Mag_Pre_Overview_plot()
    size_inches = (9, 6.5)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5,0, 1, 2))
    cn = cgrad(:seaborn_dark, 10, categorical = true)
    lcn = Makie.wong_colors()
    xd = date .+ Dates.Hour(1) .+ Dates.Minute(30) .+ Dates.Minute(15) .* range(0, 8)
    
    xtk = datetime2julian.(xd)
    xtk = Float32.(xtk .- MVN_kp[1, :JulUT])

    ax1 = Axis(fig[1, 1],  xticks=xtk, yticks=range(-15, 15, step=5), yminorticks=IntervalsBetween(5), yminorticksvisible = true, backgroundcolor=:gray80)
    ax1.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(minimum(xtk), maximum(xtk))
    ylims!(-13, 13)

    lines!(MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_mag2c1hz[!, :X_MSO], label = rich(rich("B", font = :italic), subscript("x")), color=lcn[3])
    lines!(TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_mag2c1hz[!, :X_MSO], label = rich(rich("B", font = :italic), subscript("x")), color=lcn[1])

    ax2 = Axis(fig[2, 1],  xticks=xtk, yticks=range(-15, 15, step=5), yminorticks=IntervalsBetween(5), yminorticksvisible = true, backgroundcolor=:gray80, ylabelpadding=10.0)
    ax2.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(minimum(xtk), maximum(xtk))
    ylims!(-13, 13)

    lines!(MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_mag2c1hz[!, :Y_MSO], label = rich(rich("B", font = :italic), subscript("y")), color=lcn[3])
    lines!(TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_mag2c1hz[!, :Y_MSO], label = rich(rich("B", font = :italic), subscript("y")), color=lcn[1])

    ax3 = Axis(fig[3, 1],  xticks=xtk, yticks=range(-15, 15, step=5), yminorticks=IntervalsBetween(5), yminorticksvisible = true, backgroundcolor=:gray80)
    ax3.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(minimum(xtk), maximum(xtk))
    ylims!(-14, 7)

    lines!(MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_mag2c1hz[!, :Z_MSO], label = rich(rich("B", font = :italic), subscript("z")), color=lcn[3])
    lines!(TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_mag2c1hz[!, :Z_MSO], label = rich(rich("B", font = :italic), subscript("z")), color=lcn[1])

 
    
    ax4=Axis(fig[4,1], xlabel = "Universal Time", xticks=xtk, yscale=log2)
    ax4.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(minimum(xtk), maximum(xtk)) 
    ylims!(3, 50)
    ax4.yreversed = true
    hp1 = heatmap!(MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_period1hz, log10.(MVN_Bpower1hz'), colormap=:gist_earth, colorrange=(-1.5, 2.5))
    lines!(MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_Tcp, color = :black, linestyle=:dash)


    ax5=Axis(fig[5,1], xlabel = "Universal Time", xticks=xtk, yscale=log2, ylabelpadding=10.0)
    ax5.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(minimum(xtk), maximum(xtk)) 
    ylims!(3, 50)
    ax5.yreversed = true
    hp2 = heatmap!(TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_period1hz, log10.(TW_Bpower1hz'), colormap=:gist_earth, colorrange=(-1.5, 2.5))
    lines!(TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_Tcp, color = :black, linestyle=:dash)


    ax11 = Axis(fig[1:3, 2], xgridvisible = false, ygridvisible = false)
    ylims!(1, 10)
    xlims!(0, 2)
    bracket!(0.0, 9.5, 0.0, 1.5, text = rich(rich("Tianwen-1", color=lcn[1]),  rich("     MAVEN", color=lcn[3])), orientation = :up, style = :square, width=5, rotation=pi/2)
    hidespines!(ax11)
    hidedecorations!(ax11)

    Label(fig[4, 2], rich("MAVEN", color=lcn[3]), rotation = pi/2, tellheight=false)
    Label(fig[5, 2], rich("Tianwen-1", color=lcn[1]), rotation = pi/2, tellheight=false)
    Label(fig[1:3, 0], rich("Magnetic Field", "\n", "(nT)"), rotation = pi/2, tellheight=false)
    Label(fig[4:5, 0], rich("Wave Period", "\n", "(s)"), rotation = pi/2, tellheight=false)

    #cb = Colorbar(fig[9:10, 3], hp1, label = rich("Wavelet Power"), vertical=true, flipaxis=true, tellheight=false, height=130, labelsize=15)
    #translate!(cb.blockscene, -450, 0, 0)
    #ax9.scene.viewport
    Colorbar(fig, hp1, label = rich("Log10 Power"), 
    labelsize = 15, ticklabelsize = 15, bbox=Observable(HyperRectangle{2, Int64}([98, 46], [520, 200])), vertical=false, flipaxis=false, width=140, ticks=range(-2, 3, step=1), rightspinecolor=:white, leftspinecolor=:white, topspinecolor=:white, bottomspinecolor=:white, labelpadding=-60,
    alignmode = Outside(3), halign = :right, ticklabelcolor = :white, labelcolor = :white,tickcolor = :white)

    ind = [1, 4]
    txt = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"]
    for (i, indi) in enumerate(ind)
        Label(fig[indi, -1], txt[i],font = :bold, valign=:top, tellheight=false, fontsize = 20)
    end


    ax = [ax1, ax2, ax3, ax4, ax5]
    linkxaxes!(ax...)
    hidexdecorations!.(ax[1:end-1], grid=false)
    rowgap!(fig.layout, 3)
    colgap!(fig.layout, 0)
    colsize!(fig.layout, 2, Relative(1/20))
    
    fig
end

function MVN_STATIC_eflux_Overview_plot(im)
    date = DateTime(2022, 06, 12, 00, 00, 00)
    sp = [rich(rich("H"), superscript("+")), rich(rich("He"), superscript("++")), rich("m/q=4.573966"), rich("m/q=9.22163"), rich(rich("O"), superscript("+")), rich(rich("O"), subscript("2"), superscript("+")), rich(rich("CO"), subscript("2"), superscript("+")), rich("m/q=69.26344")]

    size_inches = (20, 30)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 20, figure_padding = (5,20, 9, 10))
    cn = cgrad(:seaborn_dark, 10, categorical = true)
    lcn = Makie.wong_colors()

    xd = date .+ Dates.Hour(9) .+ Dates.Minute(30) .+ Dates.Minute(10) .* range(0, 5)
    xtk = datetime2julian.(xd)
    xtk = Float32.(xtk .- MVN_staticjulian[1])

    ind = findall(xd[begin] .<= MVN_statictime .<= xd[end])
    @show minimum(MVN_staticswp[ind]), maximum(MVN_staticswp[ind])
    swp = MVN_staticswp[ind[1]]

    axes = [Axis(fig[j, i], yscale=log10, xticks = (xtk, Dates.format.(xd, "HH:MM")), backgroundcolor=:gray80, xticklabelrotation=-π/3) for j in 1:16, i in 1:4]

    for j in 1:16
        for i in 1:4
            ax = axes[j, i]
            heatmap!(ax, Float32.(MVN_staticjulian.-MVN_staticjulian[1]), MVN_staticenergy[im, (j-1)*4+i, :, swp, 1], transpose(MVN_staticeflux[im, (j-1)*4+i, :, :]), colormap=:gist_earth, colorscale=log10, colorrange = (1e5, 1e9))
            text!(ax, 0.0, 1.0, text = rich("θ=$i ϕ=$j"), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 20, color = :white)
            xlims!(ax, xtk[begin], xtk[end]) 
            ylims!(ax, 1e1, 1e4)
            if im==1 && i==2 && (j==13 || j==14)
                annotation!(ax, -50, -40, datetime2julian(DateTime(2022, 06, 12, 10, 20, 00))-MVN_staticjulian[1], 300, path = Ann.Paths.Corner(), text = "SW", color=cn[2], fontsize=20)
                annotation!(ax, -90, -63.5, datetime2julian(DateTime(2022, 06, 12, 10, 20, 00))-MVN_staticjulian[1], 1.3e3, path = Ann.Paths.Corner(), text = "PU", color=cn[2], fontsize=20)
            end
            if im==1 && i==2 && (j==11 || j==12)
                annotation!(ax, -130, -77.5, datetime2julian(DateTime(2022, 06, 12, 10, 20, 00))-MVN_staticjulian[1], 3e3, path = Ann.Paths.Corner(), text = "SPU", color=cn[2], fontsize=20)
                annotation!(ax, -90, -58, datetime2julian(DateTime(2022, 06, 12, 10, 20, 00))-MVN_staticjulian[1], 900, path = Ann.Paths.Corner(), text = "PU", color=cn[2], fontsize=20)
            end
        end
    end
    linkaxes!(axes...)
    hidexdecorations!.(axes[1:15,:], grid=false)
    hideydecorations!.(axes[:, 2:end], grid=false)
    colgap!(fig.layout, 30)
    Colorbar(fig[1:16, 5], colormap=:gist_earth, scale=log10, colorrange = (1e5, 1e9), label = rich("MAVEN STATIC ") * sp[im] * rich(rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")))
    Label(fig[1:16, 0], rich("Energy (eV)"), fontsize = 25, halign = :center, rotation = π/2)
    Label(fig[17, 1:5], rich("Universal Time"), fontsize = 20, halign = :center)

    fig
end


function MVN_STATIC_SWIA_efluxmean_profile_plot()    
    MVN_staticvpara, MVN_staticv3para, MVN_staticefluxpara = MVN_STATIC_eflux_slice(MVN_staticvswrough, MVN_parastatic)
    MVN_staticvperp1, MVN_staticv3perp1, MVN_staticefluxperp1 = MVN_STATIC_eflux_slice(MVN_staticvswrough, MVN_perpstatic1)
    MVN_staticvperp2, MVN_staticv3perp2, MVN_staticefluxperp2 = MVN_STATIC_eflux_slice(MVN_staticvswrough, MVN_perpstatic2)

    MVN_swiavpara, MVN_swiav3para, MVN_swiaefluxpara = MVN_SWIA_eflux_slice(MVN_swiavswrough, MVN_parastatic)
    MVN_swiavperp1, MVN_swiav3perp1, MVN_swiaefluxperp1 = MVN_SWIA_eflux_slice(MVN_swiavswrough, MVN_perpstatic1)
    MVN_swiavperp2, MVN_swiav3perp2, MVN_swiaefluxperp2 = MVN_SWIA_eflux_slice(MVN_swiavswrough, MVN_perpstatic2)

    # Define the Gaussian function：p = [A, μ, σ]
    gaussian(x, p) = @. p[1] * exp(-(x - p[2])^2 / (2 * p[3]^2))

    size_inches = (9, 10)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 10, 9, 8))
    cn = cgrad(:gist_earth, 18, categorical = true)
    cn = Makie.wong_colors()
    
    ax1 = Axis(fig[1, 1], xticks=(-800:400:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("Parallel Velocity (km s"), superscript("-1"), rich(")")), ylabel = rich(rich("MAVEN STATIC H"), superscript("+"), rich("\nDifferential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)
    
    ind = findall(0.0 .< MVN_staticvpara .< 300.0)
    x = MVN_staticvpara[ind]
    y_data = MVN_staticefluxpara[ind] ./ norm.(MVN_staticv3para[ind]).^4
    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    vpara = range(-800, 800, length=81)
    Vthpara = fit.param[3]
    Tpara = ustrip(u"eV", m_p* (Vthpara * 1u"km/s")^2)
    
    v = sqrt.(vpara.^2 .+ dot(MVN_staticvswrough, MVN_perpstatic2)^2)
    eflux = gaussian(vpara, fit.param) .* v.^4

    lines!(ax1, vpara, eflux, color = cn[2], linewidth = 2)
    scatter!(ax1, MVN_staticvpara, MVN_staticefluxpara, color = cn[1]) #, linewidth = 1, linestyle = :dash

    xlims!(ax1, -750, 750)
    ylims!(ax1, 1e3, 1e9)

    ax2 = Axis(fig[1, 2], xticks=(-800:400:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("First Perpendicular Velocity (km s"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)
    
    ind = findall(-300.0 .< MVN_staticvperp1 .< 300.0)
    x = MVN_staticvperp1[ind]
    y_data = MVN_staticefluxperp1[ind] ./ norm.(MVN_staticv3perp1[ind]).^4
    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)

    Vthperp = fit.param[3]
    Tperp = ustrip(u"eV", m_p* (Vthperp * 1u"km/s")^2)

    N = ustrip(u"cm^-3", fit.param[1]*1u"cm^-2*s^-1"/(1u"km/s")*2 * (2*π)^1.5 * Vthpara * Vthperp^2)

    fit.param[2] = 0.0
    vperp1 = range(-800, 800, length=81)

    v = norm(MVN_staticvswrough)
    eflux = gaussian(vperp1, fit.param) .* v^4

    static_paramperp = deepcopy(fit.param)

    lines!(ax2, vperp1, eflux, color = cn[2], linewidth = 2)
    scatter!(ax2, MVN_staticvperp1, MVN_staticefluxperp1, color = cn[1]) #, linewidth = 1, linestyle = :dash

    xlims!(ax2, -750, 750)
    ylims!(ax2, 1e3, 1e9)

    text!(ax1, 0, 1.0, text = rich(rich("N", font=:italic), subscript("SW"), rich("="*string(round(N, digits=1))), rich(" cm"), superscript("-3"), rich("\nT", font=:italic), subscript("SW∥"), rich("="*string(round(Tpara, digits=1))), rich(" eV"), rich("\nT", font=:italic), subscript("SW⟂"), rich("="*string(round(Tperp, digits=1))), rich(" eV")), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 15, color = :black)
    

    ax3 = Axis(fig[1, 3], xticks=(-800:400:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("Second Perpendicular Velocity (km s"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)

    
    fit.param[2] = dot(MVN_staticvswrough, MVN_perpstatic2)-20.0

    v = sqrt.(vperp2.^2 .+ dot(MVN_staticvswrough, MVN_parastatic)^2)
    eflux = gaussian(vperp2, fit.param) .* v.^4

    lines!(ax3, vperp2, eflux, color = cn[2], linewidth = 2)
    scatter!(ax3, MVN_staticvperp2, MVN_staticefluxperp2, color = cn[1])
    
    xlims!(ax3, -750, 750)
    ylims!(ax3, 1e3, 1e9)

    ax4 = Axis(fig[2, 1], xticks=(-800:400:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("Parallel\nVelocity (km s"), superscript("-1"), rich(")")), ylabel = rich(rich("MAVEN SWIA Ion"), rich("\nDifferential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)
    
    ind = findall(0.0 .< MVN_swiavpara .< 300.0)
    x = MVN_swiavpara[ind]
    y_data = MVN_swiaefluxpara[ind] ./ norm.(MVN_swiav3para[ind]).^4
    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthpara = fit.param[3]
    Tpara = ustrip(u"eV", m_p* (Vthpara * 1u"km/s")^2)
    
    N = ustrip(u"cm^-3", fit.param[1]*1u"cm^-2*s^-1"/(1u"km/s")*2 * (2*π)^1.5 * Vthpara * Vthperp^2)

    @show Vthpara, Tpara, Vthperp, Tperp, N, fit.param
    
    vpara = range(-800, 800, length=81)

    v = sqrt.(vpara.^2 .+ dot(MVN_swiavswrough, MVN_perpswia2)^2)
    eflux = gaussian(vpara, fit.param) .* v.^4
    #fit.param[1] = maximum(MVN_swiaefluxpara[ind])/maximum(eflux)*fit.param[1]
    #eflux = gaussian(vpara, fit.param) .* v.^4
    maxeflux = maximum(eflux)

    lines!(ax4, vpara, eflux, color = cn[2], linewidth = 2)
    scatter!(ax4, MVN_swiavpara, MVN_swiaefluxpara, color = cn[1])

    xlims!(ax4, -750, 750)
    ylims!(ax4, 1e3, 1e9)


    ax5 = Axis(fig[2, 2], xticks=(-800:400:800), xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("1st Perpendicular\nVelocity (km s"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)

    vperp1 = range(-800, 800, length=81)

    v = norm(MVN_swiavswrough)
    param = static_paramperp
    eflux = gaussian(vperp1, param) .* v^4
    param[1] = maxeflux/maximum(eflux)*param[1]
    param[3] = param[3]*0.85
    eflux = gaussian(vperp1, param) .* v^4
    
    lines!(ax5, vperp1, eflux, color = cn[2], linewidth = 2)
    scatter!(ax5, MVN_swiavperp1, MVN_swiaefluxperp1, color = cn[1])

    Vthperp = param[3]
    Tperp = ustrip(u"eV", m_p* (Vthperp * 1u"km/s")^2)
    N = ustrip(u"cm^-3", param[1]*1u"cm^-2*s^-1"/(1u"km/s")*2 * (2*π)^1.5 * Vthpara * Vthperp^2)

    @show Vthpara, Tpara, Vthperp, Tperp, N, param
    
    xlims!(ax5, -750, 750)
    ylims!(ax5, 1e3, 1e9)

    text!(ax4, 0, 1.0, text = rich(rich("N", font=:italic), subscript("SW"), rich("="*string(round(N, digits=1))), rich(" cm"), superscript("-3"), rich("\nT", font=:italic), subscript("SW∥"), rich("="*string(round(Tpara, digits=1))), rich(" eV"), rich("\nT", font=:italic), subscript("SW⟂"), rich("="*string(round(Tperp, digits=1))), rich(" eV")), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 15, color = :black)
    

    ax6 = Axis(fig[2, 3], xticks=(-800:400:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("2nd Perpendicular\nVelocity (km s"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)

    lines!(ax6, vperp1.-20.0, eflux, color = cn[2], linewidth = 2)
    scatter!(ax6, MVN_swiavperp2, MVN_swiaefluxperp2, color = cn[1])
    
    xlims!(ax6, -750, 750)
    ylims!(ax6, 1e3, 1e9)

    ind = [1, 2]
    txt = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"]
    for (i, indi) in enumerate(ind)
        Label(fig[indi, 0], txt[i],font = :bold, valign=:top, tellheight=false, fontsize = 20)
    end
    colsize!(fig.layout, 0, Relative(1/20))
    colgap!(fig.layout, 1, -30)

    axu = [ax1, ax2, ax3]
    linkyaxes!(axu...)
    hidexdecorations!.(axu[1:end], grid=false)
    hideydecorations!.(axu[2:end], grid=false)
    axb = [ax4, ax5, ax6]
    linkxaxes!(ax1, ax4)
    linkxaxes!(ax2, ax5)
    linkxaxes!(ax3, ax6)
    linkyaxes!(axb...)
    hideydecorations!.(axb[2:end], grid=false)
    fig

    #=
    text!(ax4, 0.0, 1.0, text = rich(rich("V", font=:italic), subscript("⟂1"), rich("="), rich("V", font=:italic), subscript("E"), rich("="*string(round(dot(MVN_staticv3para[1], MVN_perpstatic1), digits=1))), rich(" km s"), superscript("-1"), rich("\n"), rich("V", font=:italic), subscript("⟂2"), rich("="), rich("V", font=:italic), subscript("V"), rich("="*string(round(dot(MVN_staticv3para[1], MVN_perpstatic2), digits=1))), rich(" km s"), superscript("-1")), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 15, color = :black)
    text!(ax5, 0.0, 1.0, text = rich(rich("V", font=:italic), subscript("∥"), rich("="), rich("V", font=:italic), subscript("B"), rich("="*string(round(dot(MVN_staticv3perp1[1], MVN_parastatic), digits=1))), rich(" km s"), superscript("-1"), rich("\n"), rich("V", font=:italic), subscript("⟂2"), rich("="), rich("V", font=:italic), subscript("V"), rich("="*string(round(dot(MVN_staticv3perp1[1], MVN_perpstatic2), digits=1))), rich(" km s"), superscript("-1")), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 15, color = :black)
    text!(ax6, 0.0, 1.0, text = rich(rich("V", font=:italic), subscript("∥"), rich("="), rich("V", font=:italic), subscript("B"), rich("="*string(round(dot(MVN_staticv3perp2[1], MVN_parastatic), digits=1))), rich(" km s"), superscript("-1"), rich("\n"), rich("V", font=:italic), subscript("⟂1"), rich("="), rich("V", font=:italic), subscript("E"), rich("="*string(round(dot(MVN_staticv3perp2[1], MVN_perpstatic1), digits=1))), rich(" km s"), superscript("-1")), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 15, color = :black)
    =#
end

function TW_MINPA_eflux_fitting_plot()
    # Define the Gaussian function：p = [A, μ, σ]
    gaussian(x, p) = @. p[1] * exp(-(x - p[2])^2 / (2 * p[3]^2))

    size_inches = (9, 12)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 10, 9, 8))
    cn = cgrad(:gist_earth, 15, categorical = true)
    cn = Makie.wong_colors()
    ax1 = Axis(fig[1, 1], xticks=(-800:200:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("Velocity Component (km s"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)
    ylims!(ax1, 1e5, 1e9)
    xlims!(ax1, -800, 800)

    TW_minpavpara, TW_minpav3para, TW_minpaefluxpara = TW_MINPA_eflux_slice(TW_minpavswrough, TW_paraminpa)
    TW_minpavperp1, TW_minpav3perp1, TW_minpaefluxperp1 = TW_MINPA_eflux_slice(TW_minpavswrough, TW_perpminpa1)
    TW_minpavperp2, TW_minpav3perp2, TW_minpaefluxperp2 = TW_MINPA_eflux_slice(TW_minpavswrough, TW_perpminpa2)

    ind = findall(isfinite.(TW_minpaefluxpara))
    ind = findall(-300.0 .< TW_minpavpara .< 200.0)
    x = TW_minpavpara[ind]
    velc = norm.(TW_minpav3para)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* TW_minpaefluxpara[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthpara = fit.param[3]
    vdpara = fit.param[2]
    Tpara = ustrip(u"eV", m_p* (Vthpara * 1u"km/s")^2)
    Npara = fit.param[1]
    v3model = [TW_minpavswrough .+ TW_paraminpa .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, TW_paraminpa) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vparamodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #lines!(ax1, vparamodel, y_model, color = cn[1], linewidth = 2)
    scatter!(ax1, TW_minpavpara, TW_minpaefluxpara, color = cn[1], markersize = 10)

    ind = findall(isfinite.(TW_minpaefluxperp2))
    ind = findall(200.0 .< TW_minpavperp2 .< 450.0)
    x = TW_minpavperp2[ind]
    velc = norm.(TW_minpav3perp2)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* TW_minpaefluxperp2[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthperp2 = fit.param[3]
    vdperp2 = fit.param[2]
    Tperp2 = ustrip(u"eV", m_p* (Vthperp2 * 1u"km/s")^2)
    Nperp2 = fit.param[1]
    v3model = [TW_minpavswrough .+ TW_perpminpa2 .* 10.0i for i in -200:200]
    vperp2model = [dot(v3, TW_perpminpa2) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vperp2model, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #lines!(ax1, vperp2model, y_model, color = cn[3], linewidth = 2)
    scatter!(ax1, TW_minpavperp2, TW_minpaefluxperp2, color = cn[3], markersize = 10)

    ind = findall(isfinite.(TW_minpaefluxperp1))
    #ind = findall(-200.0 .< TW_minpavperp1 .< 200.0)
    x = TW_minpavperp1[ind]
    velc = norm.(TW_minpav3perp1)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* TW_minpaefluxperp1[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthperp1 = fit.param[3]
    vdperp1 = fit.param[2] 
    vdperp1 = 40.0
    fit.param[2] = vdperp1
    Tperp1 = ustrip(u"eV", m_p* (Vthperp1 * 1u"km/s")^2)
    Tperp1 = Tperp2
    Vthperp1 = Vthperp2
    fit.param[3] = Vthperp2
    Nperp1 = fit.param[1]
    v3model = [TW_minpavswrough .+ TW_perpminpa1 .* 10.0i for i in -200:200]
    vperp1model = [dot(v3, TW_perpminpa1) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vperp1model, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #lines!(ax1, vperp1model, y_model, color = cn[2], linewidth = 2)
    scatter!(ax1, TW_minpavperp1, TW_minpaefluxperp1, color = cn[2], markersize = 10)
    
    ax2 = Axis(fig[2, 1], xticks=(-800:200:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("Velocity Component (km s"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)
    ylims!(ax2, 1e5, 1e9)
    xlims!(ax2, -800, 800)

    # estimation of bulk velocity of pickup ions
    i, j, k = 1, 12, 27
    TW_minpavpurough = CartesianFromSpherical()(CoordinateTransformations.Spherical(TW_minpavelocity[1, k], deg2rad(TW_minpaphi[j]), deg2rad(TW_minpatheta[i])))
    #=
    ijk = [(1, 12, 27), (1, 13, 27), (1, 14, 27), (2, 12, 26)]
    TW_minpavpurough = zeros(3)
    for (i, j, k) in ijk 
        TW_minpavpurough .+= CartesianFromSpherical()(CoordinateTransformations.Spherical(TW_minpavelocity[1, k], deg2rad(TW_minpaphi[j]), deg2rad(TW_minpatheta[i])))
    end
    TW_minpavpurough ./= length(ijk)
    =#

    vpumfa = TW_minpa2mfa * TW_minpavpurough
    vpumfaperp = [vpumfa[1], 0.0, vpumfa[3]]
    vswmfaperp = [vdperp1, 0.0, vdperp2] 

    # perpendicular velocity respect to the solar wind core
    TW_perpminpapu = TW_minpa2mfa' * normalize(vpumfaperp .- vswmfaperp)
    TW_minpavparapu, TW_minpav3parapu, TW_minpaefluxparapu = TW_MINPA_eflux_slice(TW_minpavpurough, TW_paraminpa)
    TW_minpavperppu, TW_minpav3perppu, TW_minpaefluxperppu = TW_MINPA_eflux_slice(TW_minpavpurough, TW_perpminpapu)

    ind = findall(isfinite.(TW_minpaefluxparapu))
    ind = findall(200.0 .< TW_minpavparapu .< 500.0)
    x = TW_minpavparapu[ind]
    velc = norm.(TW_minpav3parapu)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* TW_minpaefluxparapu[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthparapu = fit.param[3]
    vdparapu = fit.param[2]
    Tparapu = ustrip(u"eV", m_p* (Vthparapu * 1u"km/s")^2)
    Nparapu = fit.param[1]
    v3model = [TW_minpavpurough .+ TW_paraminpa .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, TW_paraminpa) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vparamodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, TW_minpa2mfa, v3) for v3 in v3model]
    #lines!(ax2, vparamodel, y_model, color = cn[1], linewidth = 2)
    #lines!(ax2, vparamodel, ysw_model, color = cn[1], linewidth = 2)
    scatter!(ax2, TW_minpavparapu, TW_minpaefluxparapu, color = cn[1], markersize = 10)

    fig

    ind = findall(isfinite.(TW_minpaefluxperppu))
    ind = findall(300.0 .< TW_minpavperppu .< 600.0)
    x = TW_minpavperppu[ind]
    velc = norm.(TW_minpav3perppu)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* TW_minpaefluxperppu[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthperppu = fit.param[3]
    vdperppu = fit.param[2]
    Tperppu = ustrip(u"eV", m_p* (Vthperppu * 1u"km/s")^2)
    Nperppu = fit.param[1]
    v3model = [TW_minpavpurough .+ TW_perpminpapu .* 10.0i for i in -200:200]
    vperpmodel = [dot(v3, TW_perpminpapu) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vperpmodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, TW_minpa2mfa, v3) for v3 in v3model]
    #lines!(ax2, vperpmodel, y_model, color = cn[3], linewidth = 2)
    #lines!(ax2, vperpmodel, ysw_model, color = cn[3], linewidth = 2)
    scatter!(ax2, TW_minpavperppu, TW_minpaefluxperppu, color = cn[3], markersize = 10)

    ax3 = Axis(fig[3, 1], xticks=(-800:200:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("Velocity Component (km s"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)
    ylims!(ax3, 1e5, 1e9)
    xlims!(ax3, -800, 800)

    # estiation of bulk velocity of shock-accelerated ions
    i, j, k = 1, 12, 29
    TW_minpavscrough = CartesianFromSpherical()(CoordinateTransformations.Spherical(TW_minpavelocity[1, k], deg2rad(TW_minpaphi[j]), deg2rad(TW_minpatheta[i])))
    vscmfa = TW_minpa2mfa * TW_minpavscrough
    vscmfaperp = [vscmfa[1], 0.0, vscmfa[3]]
    vswmfaperp = [vdperp1, 0.0, vdperp2] 

    # perpendicular velocity respect to the solar wind core
    TW_perpminpasc = TW_minpa2mfa' * normalize(vscmfaperp .- vswmfaperp)
    TW_minpavparasc, TW_minpav3parasc, TW_minpaefluxparasc = TW_MINPA_eflux_slice(TW_minpavscrough, TW_paraminpa)
    TW_minpavperpsc, TW_minpav3perpsc, TW_minpaefluxperpsc = TW_MINPA_eflux_slice(TW_minpavscrough, TW_perpminpasc)

    ind = findall(isfinite.(TW_minpaefluxparasc))
    ind = findall(300.0 .< TW_minpavparasc .< 600.0)
    x = TW_minpavparasc[ind]
    velc = norm.(TW_minpav3parasc)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* TW_minpaefluxparasc[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthparasc = fit.param[3]
    vdparasc = fit.param[2]
    Tparasc = ustrip(u"eV", m_p* (Vthparasc * 1u"km/s")^2)
    Nparasc = fit.param[1]

    v3model = [TW_minpavscrough .+ TW_paraminpa .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, TW_paraminpa) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vparamodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, TW_minpa2mfa, v3) for v3 in v3model]
    #lines!(ax3, vparamodel, y_model, color = cn[1], linewidth = 2)
    #lines!(ax3, vparamodel, ysw_model, color = cn[1], linewidth = 2)
    scatter!(ax3, TW_minpavparasc, TW_minpaefluxparasc, color = cn[1], markersize = 10)

    ind = findall(isfinite.(TW_minpaefluxperpsc))
    ind = findall(350.0 .< TW_minpavperpsc .< 700.0)
    x = TW_minpavperpsc[ind]
    velc = norm.(TW_minpav3perpsc)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* TW_minpaefluxperpsc[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthperpsc = fit.param[3]
    vdperpsc = fit.param[2]
    Tperpsc = ustrip(u"eV", m_p* (Vthperpsc * 1u"km/s")^2)
    Nperpsc = fit.param[1]
    v3model = [TW_minpavscrough .+ TW_perpminpasc .* 10.0i for i in -200:200]
    vperpmodel = [dot(v3, TW_perpminpasc) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vperpmodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, TW_minpa2mfa, v3) for v3 in v3model]
    #lines!(ax3, vperpmodel, y_model, color = cn[3], linewidth = 2)
    #lines!(ax3, vperpmodel, ysw_model, color = cn[3], linewidth = 2)
    scatter!(ax3, TW_minpavperpsc, TW_minpaefluxperpsc, color = cn[3], markersize = 10)

    fig

    #solar wind core
    #vdperp1 = 0.0
    Vthperp2 = Vthperp2 + 5.0
    Tperp2 = ustrip(u"eV", m_p* (Vthperp2 * 1u"km/s")^2)

    Npara = Npara * (2*π)^1.5 * Vthpara * Vthperp2^2 / exp(-0.5 * ((dot(TW_minpavswrough, TW_perpminpa1)-vdperp1)/Vthperp2)^2) / exp(-0.5 * ((dot(TW_minpavswrough, TW_perpminpa2)-vdperp2)/Vthperp2)^2)
    
    Nperp2 = Nperp2 * (2*π)^1.5 * Vthpara * Vthperp2^2 / exp(-0.5 * ((dot(TW_minpavswrough, TW_perpminpa1)-vdperp1)/Vthperp2)^2) / exp(-0.5 * ((dot(TW_minpavswrough, TW_paraminpa)-vdpara)/Vthpara)^2)

    Nperp1 = Nperp1 * (2*π)^1.5 * Vthpara * Vthperp2^2 / exp(-0.5 * ((dot(TW_minpavswrough, TW_perpminpa2)-vdperp2)/Vthperp2)^2) / exp(-0.5 * ((dot(TW_minpavswrough, TW_paraminpa)-vdpara)/Vthpara)^2)

    Nsw = Npara * 1.1

    #Nsw = 3.0
    #Vthpara = Vthpara - 0.0
    #Tpara = ustrip(u"eV", m_p* (Vthpara * 1u"km/s")^2)
    #Vthperp2 = Vthperp2 + 20.0    
    #Tperp2 = ustrip(u"eV", m_p* (Vthperppu * 1u"km/s")^2)

    dot(TW_minpavswrough, TW_paraminpa)
    dot(TW_minpavswrough, TW_perpminpa1)
    dot(TW_minpavswrough, TW_perpminpa2)
    vdpara, vdperp1, vdperp2

    TW_minpa2mfa * TW_minpavswrough

    # solar wind pickup ions
    vdparapu = vdparapu - 40.0
    Vthparapu = Vthparapu + 35.0
    Tparapu = ustrip(u"eV", m_p* (Vthparapu * 1u"km/s")^2)
    Vthperppu = Vthperppu + 20.0    
    Tperppu = ustrip(u"eV", m_p* (Vthperppu * 1u"km/s")^2)

    vr = abs(vdperppu - dot(TW_minpa2mfa * TW_perpminpapu, vswmfaperp))
    vb = vdparapu
    vb = vb - 110.0
    vr = vr - 0.0
    a = vr/sqrt(2.0)/Vthperppu
    vperptemp = dot(TW_minpavpurough, TW_perpminpapu) - dot(TW_minpa2mfa * TW_perpminpapu, vswmfaperp)
    Nparapu = Nparapu * (2*π)^1.5 * Vthparapu * Vthperppu^2 * (exp(-a^2) + sqrt(π) * a * erfc(-a)) / exp(-0.5 * ((vperptemp-vr)/Vthperppu)^2)
    
    Nperppu = Nperppu * (2*π)^1.5 * Vthperppu * Vthparapu^2 * (exp(-a^2) + sqrt(π) * a * erfc(-a)) / exp(-0.5 * ((dot(TW_minpavpurough, TW_paraminpa)-vb)/Vthparapu)^2)

    Nswpu = Nparapu*1.2

    # shock-accelerated pickup ions
    vdparasc = vdparasc + 10.0
    Vthperpsc = Vthperpsc + 40.0
    Tperpsc = ustrip(u"eV", m_p* (Vthperpsc * 1u"km/s")^2)
    Vthparasc = Vthparasc + 10.0
    Tparasc = ustrip(u"eV", m_p* (Vthparasc * 1u"km/s")^2)
    vrsc = abs(vdperpsc - dot(TW_minpa2mfa * TW_perpminpasc, vswmfaperp)) 
    vbsc = vdparasc
    vbsc = vbsc - 10.0
    asc = vrsc/sqrt(2.0)/Vthperpsc
    vperptemp = dot(TW_minpavscrough, TW_perpminpasc) - dot(TW_minpa2mfa * TW_perpminpasc, vswmfaperp)

    Nparasc = Nparasc * (2*π)^1.5 * Vthparasc * Vthperpsc^2 * (exp(-asc^2) + sqrt(π) * asc * erfc(-asc)) / exp(-0.5 * ((vperptemp-vrsc)/Vthperpsc)^2)
    
    Nperpsc = Nperpsc * (2*π)^1.5 * Vthperpsc * Vthparasc^2 * (exp(-asc^2) + sqrt(π) * asc * erfc(-asc)) / exp(-0.5 * ((dot(TW_minpavscrough, TW_paraminpa)-vbsc)/Vthparasc)^2)

    Nswpusc = Nparasc*1.4


    # solar wind core 
    v3model = [TW_minpavswrough .+ TW_paraminpa .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, TW_paraminpa) for v3 in v3model]
    vmodel = norm.(v3model)
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, TW_minpa2mfa, v3) for v3 in v3model]
    lines!(ax1, vparamodel, ysw_model, color = cn[1], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vparamodel, yswpu_model, color = cn[1], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vparamodel, yswpusc_model, color = cn[1], linewidth = 2, linestyle = :dashdot)
    lines!(ax1, vparamodel, ysw_model .+ yswpu_model .+ yswpusc_model, color = cn[1], linewidth = 2)

    v3model = [TW_minpavswrough .+ TW_perpminpa1 .* 10.0i for i in -200:200]
    vperp1model = [dot(v3, TW_perpminpa1) for v3 in v3model]
    vmodel = norm.(v3model)
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, TW_minpa2mfa, v3) for v3 in v3model]
    lines!(ax1, vperp1model, ysw_model, color = cn[2], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vperp1model, yswpu_model, color = cn[2], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vperp1model, yswpusc_model, color = cn[2], linewidth = 2, linestyle = :dashdot)
    lines!(ax1, vperp1model, ysw_model .+ yswpu_model .+ yswpusc_model, color = cn[2], linewidth = 2)

    v3model = [TW_minpavswrough .+ TW_perpminpa2 .* 10.0i for i in -200:200]
    vperp2model = [dot(v3, TW_perpminpa2) for v3 in v3model]
    vmodel = norm.(v3model)
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, TW_minpa2mfa, v3) for v3 in v3model]
    lines!(ax1, vperp2model, ysw_model, color = cn[3], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vperp2model, yswpu_model, color = cn[3], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vperp2model, yswpusc_model, color = cn[3], linewidth = 2, linestyle = :dashdot)
    lines!(ax1, vperp2model, ysw_model .+ yswpu_model .+ yswpusc_model, color = cn[3], linewidth = 2)

    
    fig

    # solar wind pickup ions
    v3model = [TW_minpavpurough .+ TW_paraminpa .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, TW_paraminpa) for v3 in v3model]
    vmodel = norm.(v3model)
    #y_model = gaussian(vparamodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, TW_minpa2mfa, v3) for v3 in v3model]
    lines!(ax2, vparamodel, ysw_model, color = cn[1], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax2, vparamodel, yswpu_model, color = cn[1], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax2, vparamodel, yswpusc_model, color = cn[1], linewidth = 2, linestyle = :dashdot)
    lines!(ax2, vparamodel, ysw_model .+ yswpu_model .+ yswpusc_model, color = cn[1], linewidth = 2)

    v3model = [TW_minpavpurough .+ TW_perpminpapu .* 10.0i for i in -200:200]
    vperpmodel = [dot(v3, TW_perpminpapu) for v3 in v3model]
    vmodel = norm.(v3model)
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, TW_minpa2mfa, v3) for v3 in v3model]
    lines!(ax2, vperpmodel, ysw_model, color = cn[3], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax2, vperpmodel, yswpu_model, color = cn[3], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax2, vperpmodel, yswpusc_model, color = cn[3], linewidth = 2, linestyle = :dashdot)
    lines!(ax2, vperpmodel, ysw_model .+ yswpu_model .+ yswpusc_model, color = cn[3], linewidth = 2)
    fig

    # shock-accelerated pickup ions
    v3model = [TW_minpavscrough .+ TW_paraminpa .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, TW_paraminpa) for v3 in v3model]
    vmodel = norm.(v3model)
    
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, TW_minpa2mfa, v3) for v3 in v3model]
    lines!(ax3, vparamodel, ysw_model, color = cn[1], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax3, vparamodel, yswpu_model, color = cn[1], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax3, vparamodel, yswpusc_model, color = cn[1], linewidth = 2, linestyle = :dashdot)
    lines!(ax3, vparamodel, ysw_model .+ yswpu_model .+ yswpusc_model, color = cn[1], linewidth = 2)

    v3model = [TW_minpavscrough .+ TW_perpminpasc .* 10.0i for i in -200:200]
    vperpmodel = [dot(v3, TW_perpminpasc) for v3 in v3model]
    vmodel = norm.(v3model)

    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, TW_minpa2mfa, v3) for v3 in v3model]
    lines!(ax3, vperpmodel, ysw_model, color = cn[3], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax3, vperpmodel, yswpu_model, color = cn[3], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax3, vperpmodel, yswpusc_model, color = cn[3], linewidth = 2, linestyle = :dashdot)
    lines!(ax3, vperpmodel, ysw_model .+ yswpu_model .+ yswpusc_model, color = cn[3], linewidth = 2)

    #=
    # magnetosheath pickup ions
    v3model = [TW_minpavmsrough .+ TW_paraminpa .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, TW_paraminpa) for v3 in v3model]
    vmodel = norm.(v3model)
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, TW_minpa2mfa, v3) for v3 in v3model]
    lines!(ax4, vparamodel, ysw_model, color = cn[1], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax4, vparamodel, yswpu_model, color = cn[1], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax4, vparamodel, yswpusc_model, color = cn[1], linewidth = 2, linestyle = :dashdot)
    lines!(ax4, vparamodel, ysw_model .+ yswpu_model .+ yswpusc_model, color = cn[1], linewidth = 2)

    v3model = [TW_minpavmsrough .+ TW_perpminpams .* 10.0i for i in -200:200]
    vperpmodel = [dot(v3, TW_perpminpams) for v3 in v3model]
    vmodel = norm.(v3model)
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, TW_minpa2mfa, v3) for v3 in v3model]
    lines!(ax4, vperpmodel, ysw_model, color = cn[3], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax4, vperpmodel, yswpu_model, color = cn[3], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax4, vperpmodel, yswpusc_model, color = cn[3], linewidth = 2, linestyle = :dashdot)
    lines!(ax4, vperpmodel, ysw_model .+ yswpu_model .+ yswpusc_model, color = cn[3], linewidth = 2)
    fig
    =#

    hidexdecorations!.([ax1, ax2], grid=false)
    linkxaxes!(ax1, ax2, ax3)
    linkyaxes!(ax1, ax2, ax3)
    Label(fig[1:3, 0], rich(rich("Tianwen-1 MINPA H"), superscript("+"), rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), rotation = π/2)

    text!(ax1, 0.0, 1.0, text = "Crossing Peak Flux Bin of Solar Wind Protons", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.9, text = "Parallel", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[1])
    text!(ax1, 0.15, 0.9, text = "1st Perpendicular", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[2])
    text!(ax1, 0.46, 0.9, text = "2nd Perpendicular", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[3])
    vsw = TW_sc2msomean * TW_minpa2sc * TW_minpa2mfa' * [vdperp1, vdpara, vdperp2]
    text!(ax1, 0.0, 0.8, text = rich("Solar Wind Core Protons"), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.7, text = rich(rich("Density " * string(round(Nsw, digits=3)) * " cm"), superscript("-3")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.6, text = rich(rich("MSO Velocity [" * string.(round.(vsw[1], digits=1)) * ", " * string.(round.(vsw[2], digits=1)) * ", " * string.(round.(vsw[3], digits=1)) * "] km s"), superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.5, text = rich("Parallel Velocity ", string.(round.(vdpara, digits=1)), " km s", superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.4, text = rich("1st Perpendicular Velocity ", string.(round.(vdperp1, digits=1)), " km s", superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.3, text = rich("2nd Perpendicular Velocity ", string.(round.(vdperp2, digits=1)), " km s", superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.2, text = rich(rich("Parallel Temperature " * string(round(Tpara, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.1, text = rich(rich("Perpendicular Temperature " * string(round(Tperp2, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)

    text!(ax2, 0.0, 1.0, text = "Crossing Peak Flux Bin of Upsteam Solar Wind Pickup Protons", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax2, 0.0, 0.9, text = "Parallel", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[1])
    text!(ax2, 0.15, 0.9, text = "Perpendicular", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[3])
    text!(ax2, 0.0, 0.8, text = rich("Upsteam Solar Wind Pickup Protons"), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax2, 0.0, 0.7, text = rich(rich("Density " * string(round(Nswpu, digits=3)) * " cm"), superscript("-3")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax2, 0.0, 0.6, text = rich(rich("Beam Velocity " * string.(round.(vb, digits=1)) * " km s"), superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax2, 0.0, 0.5, text = rich("Ring Velocity ", string.(round.(vr, digits=1)), " km s", superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)  
    text!(ax2, 0.0, 0.4, text = rich(rich("Parallel Temperature " * string(round(Tparapu, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax2, 0.0, 0.3, text = rich(rich("Perpendicular Temperature " * string(round(Tperppu, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)

    text!(ax3, 0.0, 1.0, text = "Crossing Peak Flux Bin of Shock Accelerated Pickup Protons", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax3, 0.0, 0.9, text = "Parallel", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[1])
    text!(ax3, 0.15, 0.9, text = "Perpendicular", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[3])
    text!(ax3, 0.0, 0.8, text = rich("Shock Accelerated Pickup Protons"), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax3, 0.0, 0.7, text = rich(rich("Density " * string(round(Nswpusc, digits=3)) * " cm"), superscript("-3")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax3, 0.0, 0.6, text = rich(rich("Beam Velocity " * string.(round.(vbsc, digits=1)) * " km s"), superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax3, 0.0, 0.5, text = rich("Ring Velocity ", string.(round.(vrsc, digits=1)), " km s", superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)  
    text!(ax3, 0.0, 0.4, text = rich(rich("Parallel Temperature " * string(round(Tparasc, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax3, 0.0, 0.3, text = rich(rich("Perpendicular Temperature " * string(round(Tperpsc, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)



    fig


    swcore = (Nsw = Nsw, Tpara = Tpara, Tperp2 = Tperp2, Vthpara = Vthpara, Vthperp2 = Vthperp2, vdpara = vdpara, vdperp1 = vdperp1, vdperp2 = vdperp2)
    swpu = (Nswpu = Nswpu, Tparapu = Tparapu, Tperppu = Tperppu, Vthparapu = Vthparapu, Vthperppu = Vthperppu, vb = vb, vr = vr)
    swpusc = (Nswpusc = Nswpusc, Tparasc = Tparasc, Tperpsc = Tperpsc, Vthparasc = Vthparasc, Vthperpsc = Vthperpsc, vbsc = vbsc, vrsc = vrsc)

    return fig, swcore, swpu, swpusc
   
    #=

    ind = findall(TW_minpavparapu .< -100.0)

    TW_minpavparapu[ind[end]]

    TW_minpa2mfa * TW_minpav3parapu[ind[end]]

    vs = SphericalFromCartesian()(TW_minpav3parapu[ind[end]]) 
    vtot = vs.r
    lon = vs.θ |> rad2deg
    lat = vs.ϕ |> rad2deg
    

    ind = findall(450.0 .< TW_minpavperp2 .< 700.0)

    TW_minpavperp2[ind[1]]

    TW_minpa2mfa * TW_minpav3perp2[ind[1]]
    vdperp1, vdperp2, vr

    vs = SphericalFromCartesian()(TW_minpav3perp2[ind[1]]) 
    vtot = vs.r
    lon = vs.θ |> rad2deg
    lat = vs.ϕ |> rad2deg

    ind = findall(400.0 .< TW_minpavperppu .< 700.0)

    TW_minpavperppu[ind[1]]

    TW_minpa2mfa * TW_minpav3perppu[ind[1]]

    vs = SphericalFromCartesian()(TW_minpav3perp2[ind[1]]) 
    vtot = vs.r
    lon = vs.θ |> rad2deg
    lat = vs.ϕ |> rad2deg

    TW_minpa2mfa * TW_minpavscrough

    vs = SphericalFromCartesian()(TW_minpavscrough) 
    vtot = vs.r
    lon = vs.θ |> rad2deg
    lat = vs.ϕ |> rad2deg

    TW_minpa2mfa * TW_minpavpurough

    vs = SphericalFromCartesian()(TW_minpavpurough) 
    vtot = vs.r
    lon = vs.θ |> rad2deg
    lat = vs.ϕ |> rad2deg

    @show TW_minpaphi, TW_minpatheta
    
    axu = [ax1, ax2, ax3]
    linkxaxes!(axu...)
    #hidexdecorations!.(axu[1:end], grid=false)
    #hideydecorations!.(axu[2:end], grid=false)
    fig
    =#
end


function TW_MINPA_eflux_model_data_comparision_plot(swcore, swpu, swpusc)
    Nsw, Tpara, Tperp2, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2 = swcore
    Nswpu, Tparapu, Tperppu, Vthparapu, Vthperppu, vb, vr = swpu
    Nswpusc, Tparasc, Tperpsc, Vthparasc, Vthperpsc, vbsc, vrsc = swpusc

    size_inches = (9, 18)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 20, 9, 8))
    
    cn = cgrad(:gist_earth, 18, categorical = true)

    axs = [Axis(fig[i, 1], xminorticksvisible=true, xminorticks = IntervalsBetween(9), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich("Energy (eV)"), ylabel = rich(rich("H"), superscript("+"), rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), xscale = log10, yscale=log10, backgroundcolor=:gray80) for i in 1:4]
    for i in 1:4
        for j in 1:16
            v3model = [CartesianFromSpherical()(CoordinateTransformations.Spherical(TW_minpavelocity[1, k], deg2rad(TW_minpaphi[j]), deg2rad(TW_minpatheta[i]))) for k in 1:length(TW_minpavelocity[1, :])]
            ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, TW_minpa2mfa, v3) for v3 in v3model]
            yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
            yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, TW_minpa2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
            scatterlines!(axs[i], energymod1, ysw_model .+ yswpu_model .+ yswpusc_model, color = cn[j+1], linewidth = 2)

            #ind = findall(TW_minpaefluxmean[1, j, i, :] .> 0)
            #scatterlines!(axs[i, 1], energymod1[18:end], TW_minpaefluxmean[1, j, i, 18:end], color = cn[j+1], linewidth = 3)
        end
        text!(axs[i], 0.0, 1.0, text = "Tianwen-1 MINPA θ=$i", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
        ax1_inset = Axis(fig[i, 1], width=Relative(0.8), height=Relative(0.05), halign=0.5, valign=0.02)
        phi = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], (16, 1))
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        y = [1]
        heatmap!(ax1_inset, x, y, phi, colormap=cn, colorrange=(0,17))
        hidedecorations!(ax1_inset, grid=false)
        for ix in x
            text!(ax1_inset, ix, 1.0, text = "$ix", align = (:center, :center), offset = (0, 0), fontsize = 20, color = :black)
        end
        text!(axs[i], 0.042, 0.07, text = "φ=", align = (:left, :top), offset = (0, 0), space = :relative, fontsize = 20, color = :black)
    end
    
    xlims!.(axs, 1e0, 1e4)
    ylims!.(axs, 1e3, 1e9)
    linkaxes!(axs...)
    hidexdecorations!.(axs[1:3], grid=false)

    fig
end

function MVN_STATIC_eflux_fitting_plot()
    # Define the Gaussian function：p = [A, μ, σ]
    gaussian(x, p) = @. p[1] * exp(-(x - p[2])^2 / (2 * p[3]^2))

    size_inches = (9, 18)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 10, 9, 8))
    cn = cgrad(:gist_earth, 15, categorical = true)
    cn = Makie.wong_colors()
    ax1 = Axis(fig[1, 1], xticks=(-800:200:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("Velocity Component (km s"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)
    ylims!(ax1, 1e3, 1e9)
    xlims!(ax1, -800, 800)

    MVN_staticvpara, MVN_staticv3para, MVN_staticefluxpara = MVN_STATIC_eflux_slice(MVN_staticvswrough, MVN_parastatic)
    MVN_staticvperp1, MVN_staticv3perp1, MVN_staticefluxperp1 = MVN_STATIC_eflux_slice(MVN_staticvswrough, MVN_perpstatic1)
    MVN_staticvperp2, MVN_staticv3perp2, MVN_staticefluxperp2 = MVN_STATIC_eflux_slice(MVN_staticvswrough, MVN_perpstatic2)

    ind = findall(isfinite.(MVN_staticefluxpara))
    ind = findall(0.0 .< MVN_staticvpara .< 600.0)
    x = MVN_staticvpara[ind]
    velc = norm.(MVN_staticv3para)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* MVN_staticefluxpara[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthpara = fit.param[3]
    Vthpara = Vthpara + 0.0
    fit.param[3] = Vthpara
    vdpara = fit.param[2]
    Tpara = ustrip(u"eV", m_p* (Vthpara * 1u"km/s")^2)
    Npara = fit.param[1]
    v3model = [MVN_staticvswrough .+ MVN_parastatic .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, MVN_parastatic) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vparamodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #lines!(ax1, vparamodel, y_model, color = cn[1], linewidth = 2)
    scatter!(ax1, MVN_staticvpara, MVN_staticefluxpara, color = cn[1], markersize = 10)
    fig

    ind = findall(isfinite.(MVN_staticefluxperp2))
    ind = findall(0.0 .< MVN_staticvperp2 .< 400.0)
    x = MVN_staticvperp2[ind]
    velc = norm.(MVN_staticv3perp2)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* MVN_staticefluxperp2[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthperp2 = fit.param[3]
    Vthperp2 = Vthperp2 + 20.0
    fit.param[3] = Vthperp2
    vdperp2 = fit.param[2]
    vdperp2 = vdperp2 - 40.0
    fit.param[2] = vdperp2
    Tperp2 = ustrip(u"eV", m_p* (Vthperp2 * 1u"km/s")^2)
    Nperp2 = fit.param[1]
    v3model = [MVN_staticvswrough .+ MVN_perpstatic2 .* 10.0i for i in -200:200]
    vperp2model = [dot(v3, MVN_perpstatic2) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vperp2model, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #lines!(ax1, vperp2model, y_model, color = cn[3], linewidth = 2)
    scatter!(ax1, MVN_staticvperp2, MVN_staticefluxperp2, color = cn[3], markersize = 10)

    fig

    ind = findall(isfinite.(MVN_staticefluxperp1))
    ind = findall(-200.0 .< MVN_staticvperp1 .< 200.0)
    x = MVN_staticvperp1[ind]
    velc = norm.(MVN_staticv3perp1)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* MVN_staticefluxperp1[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthperp1 = fit.param[3]
    vdperp1 = fit.param[2] 
    #vdperp1 = 40.0
    #fit.param[2] = vdperp1
    Tperp1 = ustrip(u"eV", m_p* (Vthperp1 * 1u"km/s")^2)
    Tperp1 = Tperp2
    Vthperp1 = Vthperp2
    fit.param[3] = Vthperp2
    Nperp1 = fit.param[1]
    v3model = [MVN_staticvswrough .+ MVN_perpstatic1 .* 10.0i for i in -200:200]
    vperp1model = [dot(v3, MVN_perpstatic1) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vperp1model, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #lines!(ax1, vperp1model, y_model, color = cn[2], linewidth = 2)
    scatter!(ax1, MVN_staticvperp1, MVN_staticefluxperp1, color = cn[2], markersize = 10)

    fig
    
    ax2 = Axis(fig[2, 1], xticks=(-800:200:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("Velocity Component (km s"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)
    ylims!(ax2, 1e3, 1e9)
    xlims!(ax2, -800, 800)

    # estimation of bulk velocity of pickup ions
    i, j, k = 1, 13, 9
    #i, j, k = 2, 14, 14
    MVN_staticvpurough = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, k, MVN_staticswpind, 1], deg2rad(MVN_staticphi[1, (j-1)*4+i, k, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, (j-1)*4+i, k, MVN_staticswpind, 1])))
    MVN_staticefluxmean[1, (j-1)*4+i, k]
    MVN_staticenergy[1, (j-1)*4+i, k, MVN_staticswpind, 1]
    MVN_staticenergy[1, (j-1)*4+i, :, MVN_staticswpind, 1]
    #=
    ijk = [(1, 12, 27), (1, 13, 27), (1, 14, 27), (2, 12, 26)]
    MVN_staticvpurough = zeros(3)
    for (i, j, k) in ijk 
        MVN_staticvpurough .+= CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, k], deg2rad(MVN_staticphi[j]), deg2rad(MVN_statictheta[i])))
    end
    MVN_staticvpurough ./= length(ijk)
    =#

    vpumfa = MVN_static2mfa * MVN_staticvpurough
    vpumfaperp = [vpumfa[1], 0.0, vpumfa[3]]
    vswmfaperp = [vdperp1, 0.0, vdperp2] 

    # perpendicular velocity respect to the solar wind core
    MVN_perpstaticpu = MVN_static2mfa' * normalize(vpumfaperp .- vswmfaperp)
    MVN_staticvparapu, MVN_staticv3parapu, MVN_staticefluxparapu = MVN_STATIC_eflux_slice(MVN_staticvpurough, MVN_parastatic)
    MVN_staticvperppu, MVN_staticv3perppu, MVN_staticefluxperppu = MVN_STATIC_eflux_slice(MVN_staticvpurough, MVN_perpstaticpu)

    ind = findall(isfinite.(MVN_staticefluxparapu))
    ind = findall(200.0 .< MVN_staticvparapu .< 600.0)
    x = MVN_staticvparapu[ind]
    velc = norm.(MVN_staticv3parapu)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* MVN_staticefluxparapu[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthparapu = fit.param[3]
    vdparapu = fit.param[2]
    Tparapu = ustrip(u"eV", m_p* (Vthparapu * 1u"km/s")^2)
    Nparapu = fit.param[1]
    v3model = [MVN_staticvpurough .+ MVN_parastatic .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, MVN_parastatic) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vparamodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_static2mfa, v3) for v3 in v3model]
    #lines!(ax2, vparamodel, y_model, color = cn[1], linewidth = 2)
    #lines!(ax2, vparamodel, ysw_model, color = cn[1], linewidth = 2)
    scatter!(ax2, MVN_staticvparapu, MVN_staticefluxparapu, color = cn[1], markersize = 10)

    fig

    ind = findall(isfinite.(MVN_staticefluxperppu))
    ind = findall(00.0 .< MVN_staticvperppu .< 500.0)
    x = MVN_staticvperppu[ind]
    velc = norm.(MVN_staticv3perppu)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* MVN_staticefluxperppu[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthperppu = fit.param[3]
    vdperppu = fit.param[2]
    Tperppu = ustrip(u"eV", m_p* (Vthperppu * 1u"km/s")^2)
    Nperppu = fit.param[1]
    v3model = [MVN_staticvpurough .+ MVN_perpstaticpu .* 10.0i for i in -200:200]
    vperpmodel = [dot(v3, MVN_perpstaticpu) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vperpmodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_static2mfa, v3) for v3 in v3model]
    #lines!(ax2, vperpmodel, y_model, color = cn[3], linewidth = 2)
    #lines!(ax2, vperpmodel, ysw_model, color = cn[3], linewidth = 2)
    scatter!(ax2, MVN_staticvperppu, MVN_staticefluxperppu, color = cn[3], markersize = 10)
    fig

    ax3 = Axis(fig[3, 1], xticks=(-800:200:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("Velocity Component (km s"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)
    ylims!(ax3, 1e3, 1e9)
    xlims!(ax3, -800, 800)

    # estiation of bulk velocity of shock-accelerated ions
    i, j, k = 2, 11, 7
    MVN_staticvscrough = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, k, MVN_staticswpind, 1], deg2rad(MVN_staticphi[1, (j-1)*4+i, k, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, (j-1)*4+i, k, MVN_staticswpind, 1])))
    MVN_staticefluxmean[1, (j-1)*4+i, k]
    MVN_staticenergy[1, (j-1)*4+i, k, MVN_staticswpind, 1]
    MVN_staticenergy[1, (j-1)*4+i, :, MVN_staticswpind, 1]
    vscmfa = MVN_static2mfa * MVN_staticvscrough
    vscmfaperp = [vscmfa[1], 0.0, vscmfa[3]]
    vswmfaperp = [vdperp1, 0.0, vdperp2] 

    # perpendicular velocity respect to the solar wind core
    MVN_perpstaticsc = MVN_static2mfa' * normalize(vscmfaperp .- vswmfaperp)
    MVN_staticvparasc, MVN_staticv3parasc, MVN_staticefluxparasc = MVN_STATIC_eflux_slice(MVN_staticvscrough, MVN_parastatic)
    MVN_staticvperpsc, MVN_staticv3perpsc, MVN_staticefluxperpsc = MVN_STATIC_eflux_slice(MVN_staticvscrough, MVN_perpstaticsc)

    ind = findall(isfinite.(MVN_staticefluxparasc))
    #ind = findall(300.0 .< MVN_staticvparasc .< 600.0)
    x = MVN_staticvparasc[ind]
    velc = norm.(MVN_staticv3parasc)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* MVN_staticefluxparasc[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthparasc = fit.param[3]
    vdparasc = fit.param[2]
    Tparasc = ustrip(u"eV", m_p* (Vthparasc * 1u"km/s")^2)
    Nparasc = fit.param[1]

    v3model = [MVN_staticvscrough .+ MVN_parastatic .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, MVN_parastatic) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vparamodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_static2mfa, v3) for v3 in v3model]
    #lines!(ax3, vparamodel, y_model, color = cn[1], linewidth = 2)
    #lines!(ax3, vparamodel, ysw_model, color = cn[1], linewidth = 2)
    scatter!(ax3, MVN_staticvparasc, MVN_staticefluxparasc, color = cn[1], markersize = 10)

    fig

    ind = findall(isfinite.(MVN_staticefluxperpsc))
    ind = findall(300.0 .< MVN_staticvperpsc .< 600.0)
    x = MVN_staticvperpsc[ind]
    velc = norm.(MVN_staticv3perpsc)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* MVN_staticefluxperpsc[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthperpsc = fit.param[3]
    vdperpsc = fit.param[2]
    Tperpsc = ustrip(u"eV", m_p* (Vthperpsc * 1u"km/s")^2)
    Nperpsc = fit.param[1]
    v3model = [MVN_staticvscrough .+ MVN_perpstaticsc .* 10.0i for i in -200:200]
    vperpmodel = [dot(v3, MVN_perpstaticsc) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vperpmodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_static2mfa, v3) for v3 in v3model]
    #lines!(ax3, vperpmodel, y_model, color = cn[3], linewidth = 2)
    #lines!(ax3, vperpmodel, ysw_model, color = cn[3], linewidth = 2)
    scatter!(ax3, MVN_staticvperpsc, MVN_staticefluxperpsc, color = cn[3], markersize = 10)

    fig

    ax4 = Axis(fig[4, 1], xticks=(-800:200:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("Velocity Component (km s"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)
    ylims!(ax4, 1e3, 1e9)
    xlims!(ax4, -800, 800)

    # estiation of bulk velocity of magnetosheath-pickup ions
    i, j, k = 2, 11, 15
    MVN_staticvmsrough = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, k, MVN_staticswpind, 1], deg2rad(MVN_staticphi[1, (j-1)*4+i, k, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, (j-1)*4+i, k, MVN_staticswpind, 1])))
    MVN_staticefluxmean[1, (j-1)*4+i, k]
    MVN_staticenergy[1, (j-1)*4+i, k, MVN_staticswpind, 1]
    MVN_staticenergy[1, (j-1)*4+i, :, MVN_staticswpind, 1]
    vmsmfa = MVN_static2mfa * MVN_staticvmsrough
    vmsmfaperp = [vmsmfa[1], 0.0, vmsmfa[3]]
    vswmfaperp = [vdperp1, 0.0, vdperp2] 

    # perpendicular velocity respect to the solar wind core
    MVN_perpstaticms = MVN_static2mfa' * normalize(vmsmfaperp .- vswmfaperp)
    MVN_staticvparams, MVN_staticv3params, MVN_staticefluxparams = MVN_STATIC_eflux_slice(MVN_staticvmsrough, MVN_parastatic)
    MVN_staticvperpms, MVN_staticv3perpms, MVN_staticefluxperpms = MVN_STATIC_eflux_slice(MVN_staticvmsrough, MVN_perpstaticms)

    ind = findall(isfinite.(MVN_staticefluxparams))
    ind = findall(-200.0 .< MVN_staticvparams .< 100.0)
    x = MVN_staticvparams[ind]
    velc = norm.(MVN_staticv3params)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* MVN_staticefluxparams[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthparams = fit.param[3]
    vdparams = fit.param[2]
    Tparams = ustrip(u"eV", m_p* (Vthparams * 1u"km/s")^2)
    Nparams = fit.param[1]

    v3model = [MVN_staticvmsrough .+ MVN_parastatic .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, MVN_parastatic) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vparamodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_static2mfa, v3) for v3 in v3model]
    #lines!(ax4, vparamodel, y_model, color = cn[1], linewidth = 2)
    #lines!(ax3, vparamodel, ysw_model, color = cn[1], linewidth = 2)
    scatter!(ax4, MVN_staticvparams, MVN_staticefluxparams, color = cn[1], markersize = 10)

    fig

    ind = findall(isfinite.(MVN_staticefluxperpms))
    ind = findall(0.0 .< MVN_staticvperpms .< 200.0)
    x = MVN_staticvperpms[ind]
    velc = norm.(MVN_staticv3perpms)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* MVN_staticefluxperpms[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthperpms = fit.param[3]
    vdperpms = fit.param[2]
    Tperpms = ustrip(u"eV", m_p* (Vthperpms * 1u"km/s")^2)
    Nperpms = fit.param[1]
    v3model = [MVN_staticvmsrough .+ MVN_perpstaticms .* 10.0i for i in -200:200]
    vperpmodel = [dot(v3, MVN_perpstaticms) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vperpmodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_static2mfa, v3) for v3 in v3model]
    #lines!(ax4, vperpmodel, y_model, color = cn[3], linewidth = 2)
    #lines!(ax3, vperpmodel, ysw_model, color = cn[3], linewidth = 2)
    scatter!(ax4, MVN_staticvperpms, MVN_staticefluxperpms, color = cn[3], markersize = 10)

    fig

    #solar wind core
    #vdperp1 = 0.0
    #Vthperp2 = Vthperp2 + 5.0
    vdperp1 = vdperp1 - 30.0
    vdperp2 = vdperp2 + 0.0

    #Vthpara = Vthpara + 10.0
    Tpara = ustrip(u"eV", m_p* (Vthpara * 1u"km/s")^2)
    Vthperp2 = Vthperp2 - 0.0    
    Tperp2 = ustrip(u"eV", m_p* (Vthperp2 * 1u"km/s")^2)

    Npara = Npara * (2*π)^1.5 * Vthpara * Vthperp2^2 / exp(-0.5 * ((dot(MVN_staticvswrough, MVN_perpstatic1)-vdperp1)/Vthperp2)^2) / exp(-0.5 * ((dot(MVN_staticvswrough, MVN_perpstatic2)-vdperp2)/Vthperp2)^2)
    
    Nperp2 = Nperp2 * (2*π)^1.5 * Vthpara * Vthperp2^2 / exp(-0.5 * ((dot(MVN_staticvswrough, MVN_perpstatic1)-vdperp1)/Vthperp2)^2) / exp(-0.5 * ((dot(MVN_staticvswrough, MVN_parastatic)-vdpara)/Vthpara)^2)

    Nperp1 = Nperp1 * (2*π)^1.5 * Vthpara * Vthperp2^2 / exp(-0.5 * ((dot(MVN_staticvswrough, MVN_perpstatic2)-vdperp2)/Vthperp2)^2) / exp(-0.5 * ((dot(MVN_staticvswrough, MVN_parastatic)-vdpara)/Vthpara)^2)

    Nsw = Npara #* 0.95

    #Nsw = 3.0
    

    dot(MVN_staticvswrough, MVN_parastatic)
    dot(MVN_staticvswrough, MVN_perpstatic1)
    dot(MVN_staticvswrough, MVN_perpstatic2)
    vdpara, vdperp1, vdperp2

    MVN_static2mfa * MVN_staticvswrough

    # solar wind pickup ions

    Vthparapu = Vthparapu + 20.0
    Tparapu = ustrip(u"eV", m_p* (Vthparapu * 1u"km/s")^2)
    Vthperppu = Vthperppu + 30.0    
    Tperppu = ustrip(u"eV", m_p* (Vthperppu * 1u"km/s")^2)

    vr = abs(vdperppu - dot(MVN_static2mfa * MVN_perpstaticpu, vswmfaperp))
    vb = vdparapu
    vb = vb - 60.0
    vr = vr - 20.0
    a = vr/sqrt(2.0)/Vthperppu
    vperptemp = dot(MVN_staticvpurough, MVN_perpstaticpu) - dot(MVN_static2mfa * MVN_perpstaticpu, vswmfaperp)
    Nparapu = Nparapu * (2*π)^1.5 * Vthparapu * Vthperppu^2 * (exp(-a^2) + sqrt(π) * a * erfc(-a)) / exp(-0.5 * ((vperptemp-vr)/Vthperppu)^2)
    
    Nperppu = Nperppu * (2*π)^1.5 * Vthperppu * Vthparapu^2 * (exp(-a^2) + sqrt(π) * a * erfc(-a)) / exp(-0.5 * ((dot(MVN_staticvpurough, MVN_parastatic)-vb)/Vthparapu)^2)

    Nswpu = Nparapu * 2.5


    # shock-accelerated pickup ions
    
    Vthperpsc = Vthperpsc + 40.0
    Tparasc = ustrip(u"eV", m_p* (Vthparasc * 1u"km/s")^2)
    Vthparasc = Vthparasc + 10.0
    Tperpsc = ustrip(u"eV", m_p* (Vthperpsc * 1u"km/s")^2)
    
    vrsc = abs(vdperpsc - dot(MVN_static2mfa * MVN_perpstaticsc, vswmfaperp)) 
    vbsc = vdparasc
    #vbsc = vbsc - 10.0
    vrsc = vrsc - 80.0
    asc = vrsc/sqrt(2.0)/Vthperpsc
    vperptemp = dot(MVN_staticvscrough, MVN_perpstaticsc) - dot(MVN_static2mfa * MVN_perpstaticsc, vswmfaperp)

    Nparasc = Nparasc * (2*π)^1.5 * Vthparasc * Vthperpsc^2 * (exp(-asc^2) + sqrt(π) * asc * erfc(-asc)) / exp(-0.5 * ((vperptemp-vrsc)/Vthperpsc)^2)
    
    Nperpsc = Nperpsc * (2*π)^1.5 * Vthperpsc * Vthparasc^2 * (exp(-asc^2) + sqrt(π) * asc * erfc(-asc)) / exp(-0.5 * ((dot(MVN_staticvscrough, MVN_parastatic)-vbsc)/Vthparasc)^2)

    Nswpusc = Nperpsc *1.2

    # magnetosheath pickup ions
    
    Vthperpms = Vthperpms + 0.0
    Tperpms = ustrip(u"eV", m_p* (Vthperpms * 1u"km/s")^2)  
    Vthparams = Vthparams + 50.0
    Tparams = ustrip(u"eV", m_p* (Vthparams * 1u"km/s")^2)
    
    vrms = abs(vdperpms - dot(MVN_static2mfa * MVN_perpstaticms, vswmfaperp)) 
    vbms = vdparams
    vbms = vbms + 70.0
    vrms = vrms - 40.0
    ams = vrms/sqrt(2.0)/Vthperpms
    vperptemp = dot(MVN_staticvmsrough, MVN_perpstaticms) - dot(MVN_static2mfa * MVN_perpstaticms, vswmfaperp)

    Nparams = Nparams * (2*π)^1.5 * Vthparams * Vthperpms^2 * (exp(-ams^2) + sqrt(π) * ams * erfc(-ams)) / exp(-0.5 * ((vperptemp-vrms)/Vthperpms)^2)
    
    Nperpms = Nperpms * (2*π)^1.5 * Vthperpms * Vthparams^2 * (exp(-ams^2) + sqrt(π) * ams * erfc(-ams)) / exp(-0.5 * ((dot(MVN_staticvmsrough, MVN_parastatic)-vbms)/Vthparams)^2)

    Nmspu = Nparams  * 0.2


    # solar wind core 
    v3model = [MVN_staticvswrough .+ MVN_parastatic .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, MVN_parastatic) for v3 in v3model]
    vmodel = norm.(v3model)
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_static2mfa, v3) for v3 in v3model]
    lines!(ax1, vparamodel, ysw_model, color = cn[1], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vparamodel, yswpu_model, color = cn[1], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vparamodel, yswpusc_model, color = cn[1], linewidth = 2, linestyle = :dashdot)
    ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vparamodel, ymspu_model, color = cn[1], linewidth = 2, linestyle = :dashdotdot)
    lines!(ax1, vparamodel, ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[1], linewidth = 2)

    v3model = [MVN_staticvswrough .+ MVN_perpstatic1 .* 10.0i for i in -200:200]
    vperp1model = [dot(v3, MVN_perpstatic1) for v3 in v3model]
    vmodel = norm.(v3model)
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_static2mfa, v3) for v3 in v3model]
    lines!(ax1, vperp1model, ysw_model, color = cn[2], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vperp1model, yswpu_model, color = cn[2], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vperp1model, yswpusc_model, color = cn[2], linewidth = 2, linestyle = :dashdot)
    ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vperp1model, ymspu_model, color = cn[2], linewidth = 2, linestyle = :dashdotdot)
    lines!(ax1, vperp1model, ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[2], linewidth = 2)

    v3model = [MVN_staticvswrough .+ MVN_perpstatic2 .* 10.0i for i in -200:200]
    vperp2model = [dot(v3, MVN_perpstatic2) for v3 in v3model]
    vmodel = norm.(v3model)
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_static2mfa, v3) for v3 in v3model]
    lines!(ax1, vperp2model, ysw_model, color = cn[3], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vperp2model, yswpu_model, color = cn[3], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vperp2model, yswpusc_model, color = cn[3], linewidth = 2, linestyle = :dashdot)
    ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vperp2model, ymspu_model, color = cn[3], linewidth = 2, linestyle = :dashdotdot)
    lines!(ax1, vperp2model, ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[3], linewidth = 2)

    
    fig

    # solar wind pickup ions
    v3model = [MVN_staticvpurough .+ MVN_parastatic .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, MVN_parastatic) for v3 in v3model]
    vmodel = norm.(v3model)
    #y_model = gaussian(vparamodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_static2mfa, v3) for v3 in v3model]
    lines!(ax2, vparamodel, ysw_model, color = cn[1], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax2, vparamodel, yswpu_model, color = cn[1], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax2, vparamodel, yswpusc_model, color = cn[1], linewidth = 2, linestyle = :dashdot)
    ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax2, vparamodel, ymspu_model, color = cn[1], linewidth = 2, linestyle = :dashdotdot)
    lines!(ax2, vparamodel, ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[1], linewidth = 2)

    v3model = [MVN_staticvpurough .+ MVN_perpstaticpu .* 10.0i for i in -200:200]
    vperpmodel = [dot(v3, MVN_perpstaticpu) for v3 in v3model]
    vmodel = norm.(v3model)
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_static2mfa, v3) for v3 in v3model]
    lines!(ax2, vperpmodel, ysw_model, color = cn[3], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax2, vperpmodel, yswpu_model, color = cn[3], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax2, vperpmodel, yswpusc_model, color = cn[3], linewidth = 2, linestyle = :dashdot)
    ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax2, vperpmodel, ymspu_model, color = cn[3], linewidth = 2, linestyle = :dashdotdot)
    lines!(ax2, vperpmodel, ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[3], linewidth = 2)
    fig

    # shock-accelerated pickup ions
    v3model = [MVN_staticvscrough .+ MVN_parastatic .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, MVN_parastatic) for v3 in v3model]
    vmodel = norm.(v3model)
    
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_static2mfa, v3) for v3 in v3model]
    lines!(ax3, vparamodel, ysw_model, color = cn[1], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax3, vparamodel, yswpu_model, color = cn[1], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax3, vparamodel, yswpusc_model, color = cn[1], linewidth = 2, linestyle = :dashdot)
    ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax3, vparamodel, ymspu_model, color = cn[1], linewidth = 2, linestyle = :dashdotdot)
    lines!(ax3, vparamodel, ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[1], linewidth = 2)
    fig

    v3model = [MVN_staticvscrough .+ MVN_perpstaticsc .* 10.0i for i in -200:200]
    vperpmodel = [dot(v3, MVN_perpstaticsc) for v3 in v3model]
    vmodel = norm.(v3model)

    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_static2mfa, v3) for v3 in v3model]
    lines!(ax3, vperpmodel, ysw_model, color = cn[3], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax3, vperpmodel, yswpu_model, color = cn[3], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax3, vperpmodel, yswpusc_model, color = cn[3], linewidth = 2, linestyle = :dashdot)
    ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax3, vperpmodel, ymspu_model, color = cn[3], linewidth = 2, linestyle = :dashdotdot)
    lines!(ax3, vperpmodel, ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[3], linewidth = 2)

    fig

    
    # magnetosheath pickup ions
    v3model = [MVN_staticvmsrough .+ MVN_parastatic .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, MVN_parastatic) for v3 in v3model]
    vmodel = norm.(v3model)
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_static2mfa, v3) for v3 in v3model]
    lines!(ax4, vparamodel, ysw_model, color = cn[1], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax4, vparamodel, yswpu_model, color = cn[1], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax4, vparamodel, yswpusc_model, color = cn[1], linewidth = 2, linestyle = :dashdot)
    ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax4, vparamodel, ymspu_model, color = cn[1], linewidth = 2, linestyle = :dashdotdot)
    lines!(ax4, vparamodel, ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[1], linewidth = 2)

    v3model = [MVN_staticvmsrough .+ MVN_perpstaticms .* 10.0i for i in -200:200]
    vperpmodel = [dot(v3, MVN_perpstaticms) for v3 in v3model]
    vmodel = norm.(v3model)
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_static2mfa, v3) for v3 in v3model]
    lines!(ax4, vperpmodel, ysw_model, color = cn[3], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax4, vperpmodel, yswpu_model, color = cn[3], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax4, vperpmodel, yswpusc_model, color = cn[3], linewidth = 2, linestyle = :dashdot)
    ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax4, vperpmodel, ymspu_model, color = cn[3], linewidth = 2, linestyle = :dashdotdot)
    lines!(ax4, vperpmodel, ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[3], linewidth = 2)
    fig

    
    hidexdecorations!.([ax1, ax2, ax3], grid=false)
    linkxaxes!(ax1, ax2, ax3, ax4)
    linkyaxes!(ax1, ax2, ax3, ax4)
    Label(fig[1:4, 0], rich(rich("MAVEN STATIC H"), superscript("+"), rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), rotation = π/2)
    
    fig

    text!(ax1, 0.0, 1.0, text = "Crossing Peak Flux Bin of Solar Wind Protons", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.9, text = "Parallel", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[1])
    text!(ax1, 0.15, 0.9, text = "1st Perpendicular", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[2])
    text!(ax1, 0.46, 0.9, text = "2nd Perpendicular", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[3])
    vsw = MVN_static2msomean * MVN_static2mfa' * [vdperp1, vdpara, vdperp2]
    text!(ax1, 0.0, 0.8, text = rich("Solar Wind Core Protons"), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.7, text = rich(rich("Density " * string(round(Nsw, digits=3)) * " cm"), superscript("-3")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.6, text = rich(rich("MSO Velocity [" * string.(round.(vsw[1], digits=1)) * ", " * string.(round.(vsw[2], digits=1)) * ", " * string.(round.(vsw[3], digits=1)) * "] km s"), superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.5, text = rich("Parallel Velocity ", string.(round.(vdpara, digits=1)), " km s", superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.4, text = rich("1st Perpendicular Velocity ", string.(round.(vdperp1, digits=1)), " km s", superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.3, text = rich("2nd Perpendicular Velocity ", string.(round.(vdperp2, digits=1)), " km s", superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.2, text = rich(rich("Parallel Temperature " * string(round(Tpara, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.1, text = rich(rich("Perpendicular Temperature " * string(round(Tperp2, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)

    text!(ax2, 0.0, 1.0, text = "Crossing Peak Flux Bin of Upsteam Solar Wind Pickup Protons", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax2, 0.0, 0.9, text = "Parallel", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[1])
    text!(ax2, 0.15, 0.9, text = "Perpendicular", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[3])
    text!(ax2, 0.0, 0.8, text = rich("Upsteam Solar Wind Pickup Protons"), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax2, 0.0, 0.7, text = rich(rich("Density " * string(round(Nswpu, digits=3)) * " cm"), superscript("-3")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax2, 0.0, 0.6, text = rich(rich("Beam Velocity " * string.(round.(vb, digits=1)) * " km s"), superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax2, 0.0, 0.5, text = rich("Ring Velocity ", string.(round.(vr, digits=1)), " km s", superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)  
    text!(ax2, 0.0, 0.4, text = rich(rich("Parallel Temperature " * string(round(Tparapu, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax2, 0.0, 0.3, text = rich(rich("Perpendicular Temperature " * string(round(Tperppu, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)

    text!(ax3, 0.0, 1.0, text = "Crossing Peak Flux Bin of Shock Accelerated Pickup Protons", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax3, 0.0, 0.9, text = "Parallel", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[1])
    text!(ax3, 0.15, 0.9, text = "Perpendicular", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[3])
    text!(ax3, 0.0, 0.8, text = rich("Shock Accelerated Pickup Protons"), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax3, 0.0, 0.7, text = rich(rich("Density " * string(round(Nswpusc, digits=3)) * " cm"), superscript("-3")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax3, 0.0, 0.6, text = rich(rich("Beam Velocity " * string.(round.(vbsc, digits=1)) * " km s"), superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax3, 0.0, 0.5, text = rich("Ring Velocity ", string.(round.(vrsc, digits=1)), " km s", superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)  
    text!(ax3, 0.0, 0.4, text = rich(rich("Parallel Temperature " * string(round(Tparasc, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax3, 0.0, 0.3, text = rich(rich("Perpendicular Temperature " * string(round(Tperpsc, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)

    text!(ax4, 0.0, 1.0, text = "Crossing Peak Flux Bin of Magnetosheath Pickup Protons", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax4, 0.0, 0.9, text = "Parallel", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[1])
    text!(ax4, 0.15, 0.9, text = "Perpendicular", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[3])
    text!(ax4, 0.0, 0.8, text = rich("Magnetosheath Pickup Protons"), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax4, 0.0, 0.7, text = rich(rich("Density " * string(round(Nmspu, digits=3)) * " cm"), superscript("-3")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax4, 0.0, 0.6, text = rich(rich("Beam Velocity " * string.(round.(vbms, digits=1)) * " km s"), superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax4, 0.0, 0.5, text = rich("Ring Velocity ", string.(round.(vrms, digits=1)), " km s", superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)  
    text!(ax4, 0.0, 0.4, text = rich(rich("Parallel Temperature " * string(round(Tparams, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax4, 0.0, 0.3, text = rich(rich("Perpendicular Temperature " * string(round(Tperpms, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)



    fig


    swcore = (Nsw = Nsw, Tpara = Tpara, Tperp2 = Tperp2, Vthpara = Vthpara, Vthperp2 = Vthperp2, vdpara = vdpara, vdperp1 = vdperp1, vdperp2 = vdperp2)
    swpu = (Nswpu = Nswpu, Tparapu = Tparapu, Tperppu = Tperppu, Vthparapu = Vthparapu, Vthperppu = Vthperppu, vb = vb, vr = vr)
    swpusc = (Nswpusc = Nswpusc, Tparasc = Tparasc, Tperpsc = Tperpsc, Vthparasc = Vthparasc, Vthperpsc = Vthperpsc, vbsc = vbsc, vrsc = vrsc)
    mspu = (Nmspu = Nmspu, Tparams = Tparams, Tperpms = Tperpms, Vthparams = Vthparams, Vthperpms = Vthperpms, vbms = vbms, vrms = vrms)

    return fig, swcore, swpu, swpusc, mspu
end

function MVN_STATIC_eflux_model_data_comparision_plot(swcore, swpu, swpusc, mspu)
    Nsw, Tpara, Tperp2, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2 = swcore
    Nswpu, Tparapu, Tperppu, Vthparapu, Vthperppu, vb, vr = swpu
    Nswpusc, Tparasc, Tperpsc, Vthparasc, Vthperpsc, vbsc, vrsc = swpusc
    Nmspu, Tparams, Tperpms, Vthparams, Vthperpms, vbms, vrms = mspu

    size_inches = (9, 18)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 20, 9, 8))
    cn = cgrad(:gist_earth, 18, categorical = true)

    axs = [Axis(fig[i, 1], xminorticksvisible=true, xminorticks = IntervalsBetween(9), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich("Energy (eV)"), ylabel = rich(rich("H"), superscript("+"), rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), xscale = log10, yscale=log10, backgroundcolor=:gray80) for i in 1:4]
    for i in 1:4
        for j in 1:16
            v3model = [CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, k, MVN_staticswpind, 1], deg2rad(MVN_staticphi[1, (j-1)*4+i, k, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, (j-1)*4+i, k, MVN_staticswpind, 1]))) for k in 1:length(MVN_staticvelocity[1, (j-1)*4+i, :, MVN_staticswpind, 1])]
            ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_static2mfa, v3) for v3 in v3model]
            yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
            yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
            ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_static2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
            scatterlines!(axs[i], MVN_staticenergy[1, (j-1)*4+i, :, MVN_staticswpind, 1], ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[j+1], linewidth = 2)
        end
        text!(axs[i], 0.0, 1.0, text = "MAVEN STATIC θ=$i", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
        ax1_inset = Axis(fig[i, 1], width=Relative(0.8), height=Relative(0.05), halign=0.5, valign=0.02)
        phi = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], (16, 1))
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        y = [1]
        heatmap!(ax1_inset, x, y, phi, colormap=cn, colorrange=(0,17))
        hidedecorations!(ax1_inset, grid=false)
        for ix in x
            text!(ax1_inset, ix, 1.0, text = "$ix", align = (:center, :center), offset = (0, 0), fontsize = 20, color = :black)
        end
        text!(axs[i], 0.042, 0.07, text = "φ=", align = (:left, :top), offset = (0, 0), space = :relative, fontsize = 20, color = :black)
    end
    hidexdecorations!.(axs[begin:end-1], grid=false)
    linkxaxes!(axs)
    xlims!.(axs, 1e0, 1e4)
    ylims!.(axs, 1e3, 1e9)

    fig
end

function MVN_SWIA_eflux_fitting_plot()
    # Define the Gaussian function：p = [A, μ, σ]
    gaussian(x, p) = @. p[1] * exp(-(x - p[2])^2 / (2 * p[3]^2))

    size_inches = (9, 18)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 10, 9, 8))
    cn = cgrad(:gist_earth, 15, categorical = true)
    cn = Makie.wong_colors()
    ax1 = Axis(fig[1, 1], xticks=(-800:200:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("Velocity Component (km s"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)
    ylims!(ax1, 1e3, 1e9)
    xlims!(ax1, -800, 800)

    MVN_swiavpara, MVN_swiav3para, MVN_swiaefluxpara = MVN_SWIA_eflux_slice(MVN_swiavswrough, MVN_paraswia)
    MVN_swiavperp1, MVN_swiav3perp1, MVN_swiaefluxperp1 = MVN_SWIA_eflux_slice(MVN_swiavswrough, MVN_perpswia1)
    MVN_swiavperp2, MVN_swiav3perp2, MVN_swiaefluxperp2 = MVN_SWIA_eflux_slice(MVN_swiavswrough, MVN_perpswia2)

    ind = findall(isfinite.(MVN_swiaefluxpara))
    #ind = findall(0.0 .< MVN_swiavpara .< 600.0)
    x = MVN_swiavpara[ind]
    velc = norm.(MVN_swiav3para)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* MVN_swiaefluxpara[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthpara = fit.param[3]
    Vthpara = Vthpara + 0.0
    fit.param[3] = Vthpara
    vdpara = fit.param[2]
    Tpara = ustrip(u"eV", m_p* (Vthpara * 1u"km/s")^2)
    Npara = fit.param[1]
    v3model = [MVN_swiavswrough .+ MVN_paraswia .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, MVN_paraswia) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vparamodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #lines!(ax1, vparamodel, y_model, color = cn[1], linewidth = 2)
    scatter!(ax1, MVN_swiavpara, MVN_swiaefluxpara, color = cn[1], markersize = 10)
    fig

    ind = findall(isfinite.(MVN_swiaefluxperp2))
    #ind = findall(0.0 .< MVN_swiavperp2 .< 400.0)
    x = MVN_swiavperp2[ind]
    velc = norm.(MVN_swiav3perp2)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* MVN_swiaefluxperp2[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthperp2 = fit.param[3]
    #Vthperp2 = Vthperp2 + 20.0
    #fit.param[3] = Vthperp2
    vdperp2 = fit.param[2]
    #vdperp2 = vdperp2 - 40.0
    #fit.param[2] = vdperp2
    Tperp2 = ustrip(u"eV", m_p* (Vthperp2 * 1u"km/s")^2)
    Nperp2 = fit.param[1]
    v3model = [MVN_swiavswrough .+ MVN_perpswia2 .* 10.0i for i in -200:200]
    vperp2model = [dot(v3, MVN_perpswia2) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vperp2model, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #lines!(ax1, vperp2model, y_model, color = cn[3], linewidth = 2)
    scatter!(ax1, MVN_swiavperp2, MVN_swiaefluxperp2, color = cn[3], markersize = 10)

    fig

    ind = findall(isfinite.(MVN_swiaefluxperp1))
    ind = findall(-200.0 .< MVN_swiavperp1 .< 00.0)
    x = MVN_swiavperp1[ind]
    velc = norm.(MVN_swiav3perp1)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* MVN_swiaefluxperp1[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthperp1 = fit.param[3]
    vdperp1 = fit.param[2] 
    #vdperp1 = 40.0
    #fit.param[2] = vdperp1
    Tperp1 = ustrip(u"eV", m_p* (Vthperp1 * 1u"km/s")^2)
    #Tperp1 = Tperp2
    #Vthperp1 = Vthperp2
    #fit.param[3] = Vthperp2
    Nperp1 = fit.param[1]
    v3model = [MVN_swiavswrough .+ MVN_perpswia1 .* 10.0i for i in -200:200]
    vperp1model = [dot(v3, MVN_perpswia1) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vperp1model, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #lines!(ax1, vperp1model, y_model, color = cn[2], linewidth = 2)
    scatter!(ax1, MVN_swiavperp1, MVN_swiaefluxperp1, color = cn[2], markersize = 10)

    fig
    
    ax2 = Axis(fig[2, 1], xticks=(-800:200:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("Velocity Component (km s"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)
    ylims!(ax2, 1e3, 1e9)
    xlims!(ax2, -800, 800)

    # estimation of bulk velocity of pickup ions
    i, j, k = 3, 15, 21
    #i, j, k = 2, 14, 14
    MVN_swiavpurough = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_swiavelocity[k, 1], deg2rad(MVN_swiaphi[j, 1]), deg2rad(MVN_swiatheta[i, k, 1])))
    #=
    ijk = [(1, 12, 27), (1, 13, 27), (1, 14, 27), (2, 12, 26)]
    MVN_swiavpurough = zeros(3)
    for (i, j, k) in ijk 
        MVN_swiavpurough .+= CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_swiavelocity[1, k], deg2rad(MVN_swiaphi[j]), deg2rad(MVN_swiatheta[i])))
    end
    MVN_swiavpurough ./= length(ijk)
    =#

    vpumfa = MVN_swia2mfa * MVN_swiavpurough
    vpumfaperp = [vpumfa[1], 0.0, vpumfa[3]]
    vswmfaperp = [vdperp1, 0.0, vdperp2] 

    # perpendicular velocity respect to the solar wind core
    MVN_perpswiapu = MVN_swia2mfa' * normalize(vpumfaperp .- vswmfaperp)
    MVN_swiavparapu, MVN_swiav3parapu, MVN_swiaefluxparapu = MVN_SWIA_eflux_slice(MVN_swiavpurough, MVN_paraswia)
    MVN_swiavperppu, MVN_swiav3perppu, MVN_swiaefluxperppu = MVN_SWIA_eflux_slice(MVN_swiavpurough, MVN_perpswiapu)

    ind = findall(isfinite.(MVN_swiaefluxparapu))
    #ind = findall(200.0 .< MVN_swiavparapu .< 600.0)
    x = MVN_swiavparapu[ind]
    velc = norm.(MVN_swiav3parapu)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* MVN_swiaefluxparapu[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthparapu = fit.param[3]
    vdparapu = fit.param[2]
    Tparapu = ustrip(u"eV", m_p* (Vthparapu * 1u"km/s")^2)
    Nparapu = fit.param[1]
    v3model = [MVN_swiavpurough .+ MVN_paraswia .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, MVN_paraswia) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vparamodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #lines!(ax2, vparamodel, y_model, color = cn[1], linewidth = 2)
    scatter!(ax2, MVN_swiavparapu, MVN_swiaefluxparapu, color = cn[1], markersize = 10)

    fig

    ind = findall(isfinite.(MVN_swiaefluxperppu))
    ind = findall(200.0 .< MVN_swiavperppu .< 600.0)
    x = MVN_swiavperppu[ind]
    velc = norm.(MVN_swiav3perppu)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* MVN_swiaefluxperppu[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthperppu = fit.param[3]
    vdperppu = fit.param[2]
    Tperppu = ustrip(u"eV", m_p* (Vthperppu * 1u"km/s")^2)
    Nperppu = fit.param[1]
    v3model = [MVN_swiavpurough .+ MVN_perpswiapu .* 10.0i for i in -200:200]
    vperpmodel = [dot(v3, MVN_perpswiapu) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vperpmodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #lines!(ax2, vperpmodel, y_model, color = cn[3], linewidth = 2)
    scatter!(ax2, MVN_swiavperppu, MVN_swiaefluxperppu, color = cn[3], markersize = 10)
    fig

    ax3 = Axis(fig[3, 1], xticks=(-800:200:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("Velocity Component (km s"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)
    ylims!(ax3, 1e3, 1e9)
    xlims!(ax3, -800, 800)

    # estiation of bulk velocity of shock-accelerated ions
    i, j, k = 3, 15, 14
    MVN_swiavscrough = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_swiavelocity[k, 1], deg2rad(MVN_swiaphi[j, 1]), deg2rad(MVN_swiatheta[i, k, 1])))
    vscmfa = MVN_swia2mfa * MVN_swiavscrough
    vscmfaperp = [vscmfa[1], 0.0, vscmfa[3]]
    vswmfaperp = [vdperp1, 0.0, vdperp2] 

    # perpendicular velocity respect to the solar wind core
    MVN_perpswiasc = MVN_swia2mfa' * normalize(vscmfaperp .- vswmfaperp)
    MVN_swiavparasc, MVN_swiav3parasc, MVN_swiaefluxparasc = MVN_SWIA_eflux_slice(MVN_swiavscrough, MVN_paraswia)
    MVN_swiavperpsc, MVN_swiav3perpsc, MVN_swiaefluxperpsc = MVN_SWIA_eflux_slice(MVN_swiavscrough, MVN_perpswiasc)

    ind = findall(isfinite.(MVN_swiaefluxparasc))
    #ind = findall(300.0 .< MVN_swiavparasc .< 600.0)
    x = MVN_swiavparasc[ind]
    velc = norm.(MVN_swiav3parasc)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* MVN_swiaefluxparasc[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthparasc = fit.param[3]
    vdparasc = fit.param[2]
    Tparasc = ustrip(u"eV", m_p* (Vthparasc * 1u"km/s")^2)
    Nparasc = fit.param[1]

    v3model = [MVN_swiavscrough .+ MVN_paraswia .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, MVN_paraswia) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vparamodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #lines!(ax3, vparamodel, y_model, color = cn[1], linewidth = 2)
    scatter!(ax3, MVN_swiavparasc, MVN_swiaefluxparasc, color = cn[1], markersize = 10)

    fig

    ind = findall(isfinite.(MVN_swiaefluxperpsc))
    #ind = findall(300.0 .< MVN_swiavperpsc .< 600.0)
    x = MVN_swiavperpsc[ind]
    velc = norm.(MVN_swiav3perpsc)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* MVN_swiaefluxperpsc[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthperpsc = fit.param[3]
    vdperpsc = fit.param[2]
    Tperpsc = ustrip(u"eV", m_p* (Vthperpsc * 1u"km/s")^2)
    Nperpsc = fit.param[1]
    v3model = [MVN_swiavscrough .+ MVN_perpswiasc .* 10.0i for i in -200:200]
    vperpmodel = [dot(v3, MVN_perpswiasc) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vperpmodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #lines!(ax3, vperpmodel, y_model, color = cn[3], linewidth = 2)
    scatter!(ax3, MVN_swiavperpsc, MVN_swiaefluxperpsc, color = cn[3], markersize = 10)

    fig

    ax4 = Axis(fig[4, 1], xticks=(-800:200:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("Velocity Component (km s"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)
    ylims!(ax4, 1e3, 1e9)
    xlims!(ax4, -800, 800)

    # estiation of bulk velocity of magnetosheath-pickup ions
    i, j, k = 3, 2, 33
    MVN_swiavmsrough = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_swiavelocity[k, 1], deg2rad(MVN_swiaphi[j, 1]), deg2rad(MVN_swiatheta[i, k, 1])))
    vmsmfa = MVN_swia2mfa * MVN_swiavmsrough
    vmsmfaperp = [vmsmfa[1], 0.0, vmsmfa[3]]
    vswmfaperp = [vdperp1, 0.0, vdperp2] 

    # perpendicular velocity respect to the solar wind core
    MVN_perpswiams = MVN_swia2mfa' * normalize(vmsmfaperp .- vswmfaperp)
    MVN_swiavparams, MVN_swiav3params, MVN_swiaefluxparams = MVN_SWIA_eflux_slice(MVN_swiavmsrough, MVN_paraswia)
    MVN_swiavperpms, MVN_swiav3perpms, MVN_swiaefluxperpms = MVN_SWIA_eflux_slice(MVN_swiavmsrough, MVN_perpswiams)

    ind = findall(isfinite.(MVN_swiaefluxparams))
    ind = findall(-200.0 .< MVN_swiavparams .< 0.0)
    x = MVN_swiavparams[ind]
    velc = norm.(MVN_swiav3params)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* MVN_swiaefluxparams[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthparams = fit.param[3]
    vdparams = fit.param[2]
    Tparams = ustrip(u"eV", m_p* (Vthparams * 1u"km/s")^2)
    Nparams = fit.param[1]

    v3model = [MVN_swiavmsrough .+ MVN_paraswia .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, MVN_paraswia) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vparamodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #lines!(ax4, vparamodel, y_model, color = cn[1], linewidth = 2)
    scatter!(ax4, MVN_swiavparams, MVN_swiaefluxparams, color = cn[1], markersize = 10)

    fig

    ind = findall(isfinite.(MVN_swiaefluxperpms))
    ind = findall(0.0 .< MVN_swiavperpms .< 200.0)
    x = MVN_swiavperpms[ind]
    velc = norm.(MVN_swiav3perpms)
    y_data = 2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4) .* MVN_swiaefluxperpms[ind] ./ velc[ind].^4 

    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthperpms = fit.param[3]
    vdperpms = fit.param[2]
    Tperpms = ustrip(u"eV", m_p* (Vthperpms * 1u"km/s")^2)
    Nperpms = fit.param[1]
    v3model = [MVN_swiavmsrough .+ MVN_perpswiams .* 10.0i for i in -200:200]
    vperpmodel = [dot(v3, MVN_perpswiams) for v3 in v3model]
    vmodel = norm.(v3model)
    y_model = gaussian(vperpmodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    #lines!(ax4, vperpmodel, y_model, color = cn[3], linewidth = 2)
    scatter!(ax4, MVN_swiavperpms, MVN_swiaefluxperpms, color = cn[3], markersize = 10)

    fig

    #solar wind core
    vdperp1 = -10.0
    #vdperp2 = vdperp2 + 0.0

    #Vthpara = Vthpara + 10.0
    Tpara = ustrip(u"eV", m_p* (Vthpara * 1u"km/s")^2)
    Vthperp2 = Vthperp2 + 0.0   
    Tperp2 = ustrip(u"eV", m_p* (Vthperp2 * 1u"km/s")^2)

    Npara = Npara * (2*π)^1.5 * Vthpara * Vthperp2^2 / exp(-0.5 * ((dot(MVN_swiavswrough, MVN_perpswia1)-vdperp1)/Vthperp2)^2) / exp(-0.5 * ((dot(MVN_swiavswrough, MVN_perpswia2)-vdperp2)/Vthperp2)^2)
    
    Nperp2 = Nperp2 * (2*π)^1.5 * Vthpara * Vthperp2^2 / exp(-0.5 * ((dot(MVN_swiavswrough, MVN_perpswia1)-vdperp1)/Vthperp2)^2) / exp(-0.5 * ((dot(MVN_swiavswrough, MVN_paraswia)-vdpara)/Vthpara)^2)

    Nperp1 = Nperp1 * (2*π)^1.5 * Vthpara * Vthperp2^2 / exp(-0.5 * ((dot(MVN_swiavswrough, MVN_perpswia2)-vdperp2)/Vthperp2)^2) / exp(-0.5 * ((dot(MVN_swiavswrough, MVN_paraswia)-vdpara)/Vthpara)^2)

    Nsw = Npara #* 0.95

    #Nsw = 3.0
    

    dot(MVN_swiavswrough, MVN_paraswia)
    dot(MVN_swiavswrough, MVN_perpswia1)
    dot(MVN_swiavswrough, MVN_perpswia2)
    vdpara, vdperp1, vdperp2

    MVN_swia2mfa * MVN_swiavswrough

    # solar wind pickup ions

    #Vthparapu = Vthparapu + 20.0
    Tparapu = ustrip(u"eV", m_p* (Vthparapu * 1u"km/s")^2)
    Vthperppu = Vthperppu - 30.0    
    Tperppu = ustrip(u"eV", m_p* (Vthperppu * 1u"km/s")^2)

    vr = abs(vdperppu - dot(MVN_swia2mfa * MVN_perpswiapu, vswmfaperp))
    vb = vdparapu
    #vb = vb - 60.0
    #vr = vr - 20.0
    a = vr/sqrt(2.0)/Vthperppu
    vperptemp = dot(MVN_swiavpurough, MVN_perpswiapu) - dot(MVN_swia2mfa * MVN_perpswiapu, vswmfaperp)
    Nparapu = Nparapu * (2*π)^1.5 * Vthparapu * Vthperppu^2 * (exp(-a^2) + sqrt(π) * a * erfc(-a)) / exp(-0.5 * ((vperptemp-vr)/Vthperppu)^2)
    
    Nperppu = Nperppu * (2*π)^1.5 * Vthperppu * Vthparapu^2 * (exp(-a^2) + sqrt(π) * a * erfc(-a)) / exp(-0.5 * ((dot(MVN_swiavpurough, MVN_paraswia)-vb)/Vthparapu)^2)

    Nswpu = Nparapu * 0.2


    # shock-accelerated pickup ions
    
    #Vthperpsc = Vthperpsc + 40.0
    Tparasc = ustrip(u"eV", m_p* (Vthparasc * 1u"km/s")^2)
    #Vthparasc = Vthparasc + 10.0
    Tperpsc = ustrip(u"eV", m_p* (Vthperpsc * 1u"km/s")^2)
    
    vrsc = abs(vdperpsc - dot(MVN_swia2mfa * MVN_perpswiasc, vswmfaperp)) 
    vbsc = vdparasc
    #vbsc = vbsc - 10.0
    #vrsc = vrsc - 80.0
    asc = vrsc/sqrt(2.0)/Vthperpsc
    vperptemp = dot(MVN_swiavscrough, MVN_perpswiasc) - dot(MVN_swia2mfa * MVN_perpswiasc, vswmfaperp)

    Nparasc = Nparasc * (2*π)^1.5 * Vthparasc * Vthperpsc^2 * (exp(-asc^2) + sqrt(π) * asc * erfc(-asc)) / exp(-0.5 * ((vperptemp-vrsc)/Vthperpsc)^2)
    
    Nperpsc = Nperpsc * (2*π)^1.5 * Vthperpsc * Vthparasc^2 * (exp(-asc^2) + sqrt(π) * asc * erfc(-asc)) / exp(-0.5 * ((dot(MVN_swiavscrough, MVN_paraswia)-vbsc)/Vthparasc)^2)

    Nswpusc = Nperpsc * 0.3

    # magnetosheath pickup ions
    
    #Vthperpms = Vthperpms + 0.0
    Tperpms = ustrip(u"eV", m_p* (Vthperpms * 1u"km/s")^2)  
    #Vthparams = Vthparams + 50.0
    Tparams = ustrip(u"eV", m_p* (Vthparams * 1u"km/s")^2)
    
    vrms = abs(vdperpms - dot(MVN_swia2mfa * MVN_perpswiams, vswmfaperp)) 
    vbms = vdparams
    #vbms = vbms + 70.0
    #vrms = vrms - 40.0
    ams = vrms/sqrt(2.0)/Vthperpms
    vperptemp = dot(MVN_swiavmsrough, MVN_perpswiams) - dot(MVN_swia2mfa * MVN_perpswiams, vswmfaperp)

    Nparams = Nparams * (2*π)^1.5 * Vthparams * Vthperpms^2 * (exp(-ams^2) + sqrt(π) * ams * erfc(-ams)) / exp(-0.5 * ((vperptemp-vrms)/Vthperpms)^2)
    
    Nperpms = Nperpms * (2*π)^1.5 * Vthperpms * Vthparams^2 * (exp(-ams^2) + sqrt(π) * ams * erfc(-ams)) / exp(-0.5 * ((dot(MVN_swiavmsrough, MVN_paraswia)-vbms)/Vthparams)^2)

    Nmspu = Nparams  #* 0.2


    # solar wind core 
    v3model = [MVN_swiavswrough .+ MVN_paraswia .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, MVN_paraswia) for v3 in v3model]
    vmodel = norm.(v3model)
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_swia2mfa, v3) for v3 in v3model]
    lines!(ax1, vparamodel, ysw_model, color = cn[1], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vparamodel, yswpu_model, color = cn[1], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vparamodel, yswpusc_model, color = cn[1], linewidth = 2, linestyle = :dashdot)
    ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vparamodel, ymspu_model, color = cn[1], linewidth = 2, linestyle = :dashdotdot)
    lines!(ax1, vparamodel, ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[1], linewidth = 2)

    v3model = [MVN_swiavswrough .+ MVN_perpswia1 .* 10.0i for i in -200:200]
    vperp1model = [dot(v3, MVN_perpswia1) for v3 in v3model]
    vmodel = norm.(v3model)
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_swia2mfa, v3) for v3 in v3model]
    lines!(ax1, vperp1model, ysw_model, color = cn[2], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vperp1model, yswpu_model, color = cn[2], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vperp1model, yswpusc_model, color = cn[2], linewidth = 2, linestyle = :dashdot)
    ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vperp1model, ymspu_model, color = cn[2], linewidth = 2, linestyle = :dashdotdot)
    lines!(ax1, vperp1model, ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[2], linewidth = 2)

    v3model = [MVN_swiavswrough .+ MVN_perpswia2 .* 10.0i for i in -200:200]
    vperp2model = [dot(v3, MVN_perpswia2) for v3 in v3model]
    vmodel = norm.(v3model)
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_swia2mfa, v3) for v3 in v3model]
    lines!(ax1, vperp2model, ysw_model, color = cn[3], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vperp2model, yswpu_model, color = cn[3], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vperp2model, yswpusc_model, color = cn[3], linewidth = 2, linestyle = :dashdot)
    ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax1, vperp2model, ymspu_model, color = cn[3], linewidth = 2, linestyle = :dashdotdot)
    lines!(ax1, vperp2model, ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[3], linewidth = 2)

    
    fig

    # solar wind pickup ions
    v3model = [MVN_swiavpurough .+ MVN_paraswia .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, MVN_paraswia) for v3 in v3model]
    vmodel = norm.(v3model)
    #y_model = gaussian(vparamodel, fit.param) .* vmodel.^4 ./ (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_swia2mfa, v3) for v3 in v3model]
    lines!(ax2, vparamodel, ysw_model, color = cn[1], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax2, vparamodel, yswpu_model, color = cn[1], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax2, vparamodel, yswpusc_model, color = cn[1], linewidth = 2, linestyle = :dashdot)
    ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax2, vparamodel, ymspu_model, color = cn[1], linewidth = 2, linestyle = :dashdotdot)
    lines!(ax2, vparamodel, ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[1], linewidth = 2)

    v3model = [MVN_swiavpurough .+ MVN_perpswiapu .* 10.0i for i in -200:200]
    vperpmodel = [dot(v3, MVN_perpswiapu) for v3 in v3model]
    vmodel = norm.(v3model)
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_swia2mfa, v3) for v3 in v3model]
    lines!(ax2, vperpmodel, ysw_model, color = cn[3], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax2, vperpmodel, yswpu_model, color = cn[3], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax2, vperpmodel, yswpusc_model, color = cn[3], linewidth = 2, linestyle = :dashdot)
    ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax2, vperpmodel, ymspu_model, color = cn[3], linewidth = 2, linestyle = :dashdotdot)
    lines!(ax2, vperpmodel, ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[3], linewidth = 2)
    fig

    # shock-accelerated pickup ions
    v3model = [MVN_swiavscrough .+ MVN_paraswia .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, MVN_paraswia) for v3 in v3model]
    vmodel = norm.(v3model)
    
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_swia2mfa, v3) for v3 in v3model]
    lines!(ax3, vparamodel, ysw_model, color = cn[1], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax3, vparamodel, yswpu_model, color = cn[1], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax3, vparamodel, yswpusc_model, color = cn[1], linewidth = 2, linestyle = :dashdot)
    ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax3, vparamodel, ymspu_model, color = cn[1], linewidth = 2, linestyle = :dashdotdot)
    lines!(ax3, vparamodel, ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[1], linewidth = 2)
    fig

    v3model = [MVN_swiavscrough .+ MVN_perpswiasc .* 10.0i for i in -200:200]
    vperpmodel = [dot(v3, MVN_perpswiasc) for v3 in v3model]
    vmodel = norm.(v3model)

    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_swia2mfa, v3) for v3 in v3model]
    lines!(ax3, vperpmodel, ysw_model, color = cn[3], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax3, vperpmodel, yswpu_model, color = cn[3], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax3, vperpmodel, yswpusc_model, color = cn[3], linewidth = 2, linestyle = :dashdot)
    ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax3, vperpmodel, ymspu_model, color = cn[3], linewidth = 2, linestyle = :dashdotdot)
    lines!(ax3, vperpmodel, ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[3], linewidth = 2)

    fig

    
    # magnetosheath pickup ions
    v3model = [MVN_swiavmsrough .+ MVN_paraswia .* 10.0i for i in -200:200]
    vparamodel = [dot(v3, MVN_paraswia) for v3 in v3model]
    vmodel = norm.(v3model)
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_swia2mfa, v3) for v3 in v3model]
    lines!(ax4, vparamodel, ysw_model, color = cn[1], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax4, vparamodel, yswpu_model, color = cn[1], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax4, vparamodel, yswpusc_model, color = cn[1], linewidth = 2, linestyle = :dashdot)
    ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax4, vparamodel, ymspu_model, color = cn[1], linewidth = 2, linestyle = :dashdotdot)
    lines!(ax4, vparamodel, ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[1], linewidth = 2)

    v3model = [MVN_swiavmsrough .+ MVN_perpswiams .* 10.0i for i in -200:200]
    vperpmodel = [dot(v3, MVN_perpswiams) for v3 in v3model]
    vmodel = norm.(v3model)
    ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_swia2mfa, v3) for v3 in v3model]
    lines!(ax4, vperpmodel, ysw_model, color = cn[3], linewidth = 2, linestyle = :dot)
    yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax4, vperpmodel, yswpu_model, color = cn[3], linewidth = 2, linestyle = :dash)
    yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax4, vperpmodel, yswpusc_model, color = cn[3], linewidth = 2, linestyle = :dashdot)
    ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
    lines!(ax4, vperpmodel, ymspu_model, color = cn[3], linewidth = 2, linestyle = :dashdotdot)
    lines!(ax4, vperpmodel, ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[3], linewidth = 2)
    fig

    
    hidexdecorations!.([ax1, ax2, ax3], grid=false)
    linkxaxes!(ax1, ax2, ax3, ax4)
    linkyaxes!(ax1, ax2, ax3, ax4)
    Label(fig[1:4, 0], rich(rich("MAVEN SWIA H"), superscript("+"), rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), rotation = π/2)
    
    fig

    text!(ax1, 0.0, 1.0, text = "Crossing Peak Flux Bin of Solar Wind Protons", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.9, text = "Parallel", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[1])
    text!(ax1, 0.15, 0.9, text = "1st Perpendicular", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[2])
    text!(ax1, 0.46, 0.9, text = "2nd Perpendicular", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[3])
    vsw = MVN_swia2msomean * MVN_swia2mfa' * [vdperp1, vdpara, vdperp2]
    text!(ax1, 0.0, 0.8, text = rich("Solar Wind Core Protons"), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.7, text = rich(rich("Density " * string(round(Nsw, digits=3)) * " cm"), superscript("-3")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.6, text = rich(rich("MSO Velocity [" * string.(round.(vsw[1], digits=1)) * ", " * string.(round.(vsw[2], digits=1)) * ", " * string.(round.(vsw[3], digits=1)) * "] km s"), superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.5, text = rich("Parallel Velocity ", string.(round.(vdpara, digits=1)), " km s", superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.4, text = rich("1st Perpendicular Velocity ", string.(round.(vdperp1, digits=1)), " km s", superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.3, text = rich("2nd Perpendicular Velocity ", string.(round.(vdperp2, digits=1)), " km s", superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.2, text = rich(rich("Parallel Temperature " * string(round(Tpara, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax1, 0.0, 0.1, text = rich(rich("Perpendicular Temperature " * string(round(Tperp2, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)

    text!(ax2, 0.0, 1.0, text = "Crossing Peak Flux Bin of Upsteam Solar Wind Pickup Protons", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax2, 0.0, 0.9, text = "Parallel", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[1])
    text!(ax2, 0.15, 0.9, text = "Perpendicular", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[3])
    text!(ax2, 0.0, 0.8, text = rich("Upsteam Solar Wind Pickup Protons"), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax2, 0.0, 0.7, text = rich(rich("Density " * string(round(Nswpu, digits=3)) * " cm"), superscript("-3")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax2, 0.0, 0.6, text = rich(rich("Beam Velocity " * string.(round.(vb, digits=1)) * " km s"), superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax2, 0.0, 0.5, text = rich("Ring Velocity ", string.(round.(vr, digits=1)), " km s", superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)  
    text!(ax2, 0.0, 0.4, text = rich(rich("Parallel Temperature " * string(round(Tparapu, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax2, 0.0, 0.3, text = rich(rich("Perpendicular Temperature " * string(round(Tperppu, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)

    text!(ax3, 0.0, 1.0, text = "Crossing Peak Flux Bin of Shock Accelerated Pickup Protons", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax3, 0.0, 0.9, text = "Parallel", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[1])
    text!(ax3, 0.15, 0.9, text = "Perpendicular", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[3])
    text!(ax3, 0.0, 0.8, text = rich("Shock Accelerated Pickup Protons"), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax3, 0.0, 0.7, text = rich(rich("Density " * string(round(Nswpusc, digits=3)) * " cm"), superscript("-3")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax3, 0.0, 0.6, text = rich(rich("Beam Velocity " * string.(round.(vbsc, digits=1)) * " km s"), superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax3, 0.0, 0.5, text = rich("Ring Velocity ", string.(round.(vrsc, digits=1)), " km s", superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)  
    text!(ax3, 0.0, 0.4, text = rich(rich("Parallel Temperature " * string(round(Tparasc, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax3, 0.0, 0.3, text = rich(rich("Perpendicular Temperature " * string(round(Tperpsc, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)

    text!(ax4, 0.0, 1.0, text = "Crossing Peak Flux Bin of Magnetosheath Pickup Protons", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax4, 0.0, 0.9, text = "Parallel", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[1])
    text!(ax4, 0.15, 0.9, text = "Perpendicular", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = cn[3])
    text!(ax4, 0.0, 0.8, text = rich("Magnetosheath Pickup Protons"), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax4, 0.0, 0.7, text = rich(rich("Density " * string(round(Nmspu, digits=3)) * " cm"), superscript("-3")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax4, 0.0, 0.6, text = rich(rich("Beam Velocity " * string.(round.(vbms, digits=1)) * " km s"), superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax4, 0.0, 0.5, text = rich("Ring Velocity ", string.(round.(vrms, digits=1)), " km s", superscript("-1")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)  
    text!(ax4, 0.0, 0.4, text = rich(rich("Parallel Temperature " * string(round(Tparams, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    text!(ax4, 0.0, 0.3, text = rich(rich("Perpendicular Temperature " * string(round(Tperpms, digits=1)) * " eV")), align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)



    fig


    swcore = (Nsw = Nsw, Tpara = Tpara, Tperp2 = Tperp2, Vthpara = Vthpara, Vthperp2 = Vthperp2, vdpara = vdpara, vdperp1 = vdperp1, vdperp2 = vdperp2)
    swpu = (Nswpu = Nswpu, Tparapu = Tparapu, Tperppu = Tperppu, Vthparapu = Vthparapu, Vthperppu = Vthperppu, vb = vb, vr = vr)
    swpusc = (Nswpusc = Nswpusc, Tparasc = Tparasc, Tperpsc = Tperpsc, Vthparasc = Vthparasc, Vthperpsc = Vthperpsc, vbsc = vbsc, vrsc = vrsc)
    mspu = (Nmspu = Nmspu, Tparams = Tparams, Tperpms = Tperpms, Vthparams = Vthparams, Vthperpms = Vthperpms, vbms = vbms, vrms = vrms)

    return fig, swcore, swpu, swpusc, mspu
end

function MVN_SWIA_eflux_model_data_comparision_plot(swcore, swpu, swpusc, mspu)
    Nsw, Tpara, Tperp2, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2 = swcore
    Nswpu, Tparapu, Tperppu, Vthparapu, Vthperppu, vb, vr = swpu
    Nswpusc, Tparasc, Tperpsc, Vthparasc, Vthperpsc, vbsc, vrsc = swpusc
    Nmspu, Tparams, Tperpms, Vthparams, Vthperpms, vbms, vrms = mspu

    size_inches = (9, 18)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 20, 9, 8))
    cn = cgrad(:gist_earth, 18, categorical = true)

    axs = [Axis(fig[i, 1], xminorticksvisible=true, xminorticks = IntervalsBetween(9), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich("Energy (eV)"), ylabel = rich(rich("H"), superscript("+"), rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), xscale = log10, yscale=log10, backgroundcolor=:gray80) for i in 1:4]
    for i in 1:4
        for j in 1:16
            v3model = [CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_swiavelocity[1, (j-1)*4+i, k, MVN_swiaswpind, 1], deg2rad(MVN_swiaphi[1, (j-1)*4+i, k, MVN_swiaswpind, 1]), deg2rad(MVN_swiatheta[1, (j-1)*4+i, k, MVN_swiaswpind, 1]))) for k in 1:length(MVN_swiavelocity[1, (j-1)*4+i, :, MVN_swiaswpind, 1])]
            ysw_model = [swmaxwellian(Nsw, Vthpara, Vthperp2, vdpara, vdperp1, vdperp2, MVN_swia2mfa, v3) for v3 in v3model]
            yswpu_model = [puringbeam(Nswpu, Vthparapu, Vthperppu, vb, vr, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
            yswpusc_model = [puringbeam(Nswpusc, Vthparasc, Vthperpsc, vbsc, vrsc, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
            ymspu_model = [puringbeam(Nmspu, Vthparams, Vthperpms, vbms, vrms, MVN_swia2mfa, vdperp1, vdperp2, v3) for v3 in v3model]
            scatterlines!(axs[i], MVN_swiaenergy[1, (j-1)*4+i, :, MVN_swiaswpind, 1], ysw_model .+ yswpu_model .+ yswpusc_model .+ ymspu_model, color = cn[j+1], linewidth = 2)
        end
        text!(axs[i], 0.0, 1.0, text = "MAVEN SWIA θ=$i", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
        ax1_inset = Axis(fig[i, 1], width=Relative(0.8), height=Relative(0.05), halign=0.5, valign=0.02)
        phi = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], (16, 1))
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        y = [1]
        heatmap!(ax1_inset, x, y, phi, colormap=cn, colorrange=(0,17))
        hidedecorations!(ax1_inset, grid=false)
        for ix in x
            text!(ax1_inset, ix, 1.0, text = "$ix", align = (:center, :center), offset = (0, 0), fontsize = 20, color = :black)
        end
        text!(axs[i], 0.042, 0.07, text = "φ=", align = (:left, :top), offset = (0, 0), space = :relative, fontsize = 20, color = :black)
    end
    hidexdecorations!.(axs[begin:end-1], grid=false)
    linkxaxes!(axs)
    xlims!.(axs, 1e0, 1e4)
    ylims!.(axs, 1e3, 1e9)

    fig
end

function MVN_STATIC_efluxmean_heatmap_plot()
    size_inches = (30, 30)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 10, 9, 8))
    cn = cgrad(:gist_earth, 15, categorical = true)

    axs = [Axis(fig[j, i], xticks=(-800:100:800), xminorticksvisible=true, xminorticks = IntervalsBetween(4), yticks=(-800:100:800), yminorticksvisible=true, yminorticks = IntervalsBetween(4), backgroundcolor=:gray80, aspect=DataAspect()) for i in 1:3, j in 1:3]

    for i in 1:3
        for j in 1:3
            ipara = ((j-1)*3+i-1)*5+90
            vparai = vpara[ipara]
            heatmap!(axs[i, j], vperp1, vperp2, MVN_static_eflux_rect[:, ipara, :], colormap = :gist_earth,  colorrange = (1e3, 1e9), colorscale = log10)
            viz!(axs[i, j], Meshes.Circle(Plane((0, dot(MVN_staticvswrough, MVN_perpstatic2), 0), (0, 0, 1)), 300), transparency=true, color =:blue)
        end
    end
    fig
end

function TW_MINPA_efluxmean_slice_plot()
    using GLMakie
    GLMakie.activate!()

    size_inches = (9, 20)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 10, 9, 8))
    ls1 = LScene(fig[1, 1])

    #for ipara in 80:10:120 #length(vpara)
    #for ipara in 1:length(vpara)
    for ipara in 81:3:90
    #for ipara in 80:10:100
    #for ipara in 108:3:111
        vparai = vpara[ipara]
        heatmap!(ls1, vperp1, vperp2, TW_minpa_eflux_rect[:, ipara, :], transformation = (:xz, vparai), colormap = :gist_earth,  colorrange = (1e4, 1e8), colorscale = log10)
    end

    i, j = 1, 13
    ind = argmax(TW_minpaefluxmean[1, j, i, :])
    v1 = CartesianFromSpherical()(CoordinateTransformations.Spherical(TW_minpavelocity[1, ind], deg2rad(TW_minpaphi[j]), deg2rad(TW_minpatheta[i])))
    v1 = CartesianFromSpherical()(CoordinateTransformations.Spherical(TW_minpavelocity[1, 22], deg2rad(TW_minpaphi[j]), deg2rad(TW_minpatheta[i]+22.5)))
    ps = Point3f(0.0, 0.0, 0.0)
    pe = Point3f(dot(v1, TW_perpminpa1), dot(v1, TW_paraminpa), dot(v1, TW_perpminpa2))
    arrows3d!(ls1, ps, pe, argmode = :endpoint, color = :red, shaftradius = 0.02, tipradius = 0.05, tailradius = 0.02)

    ls2 = LScene(fig[1, 2])

    for iperp1 in 81:1:81
        vperp1i = vperp1[iperp1]
        heatmap!(ls2, vpara, vperp2, TW_minpa_eflux_rect[iperp1, :, :], transformation = (:yz, vperp1i), colormap = :gist_earth,  colorrange = (1e4, 1e8), colorscale = log10)
    end

    i, j = 1, 13
    ind = argmax(TW_minpaefluxmean[1, j, i, :])
    v1 = CartesianFromSpherical()(CoordinateTransformations.Spherical(TW_minpavelocity[1, ind], deg2rad(TW_minpaphi[j]), deg2rad(TW_minpatheta[i])))
    v1 = CartesianFromSpherical()(CoordinateTransformations.Spherical(TW_minpavelocity[1, 22], deg2rad(TW_minpaphi[j]), deg2rad(TW_minpatheta[i]+22.5)))
    ps = Point3f(0.0, 0.0, 0.0)
    pe = Point3f(dot(v1, TW_perpminpa1), dot(v1, TW_paraminpa), dot(v1, TW_perpminpa2))
    arrows3d!(ls2, ps, pe, argmode = :endpoint, color = :red, shaftradius = 0.02, tipradius = 0.05, tailradius = 0.02)
    fig

    maximum(filter(x -> !isnan(x), TW_minpa_eflux_rect))
end

function MVN_STATIC_SWIA_efluxmean_slice_plot()
    using GLMakie
    GLMakie.activate!()

    size_inches = (9, 20)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 10, 9, 8))
    ls1 = LScene(fig[1, 1])

    #for ipara in 80:10:120 #length(vpara)
    for ipara in 1:length(vpara)
    #for ipara in 90:5:110
        vparai = vpara[ipara]
        heatmap!(ls1, vperp1, vperp2, MVN_static_eflux_rect[:, ipara, :], transformation = (:xz, vparai), colormap = :gist_earth,  colorrange = (1e3, 1e9), colorscale = log10)
    end

    #=
    #probably associated with the reflected pickup protons
    i, j =2, 11
    ind = argmax(MVN_staticefluxmean[1, (j-1)*4+i, :])
    ind = 7
    v1 = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, ind, MVN_staticswpind, 1], deg2rad(MVN_staticphi[1, (j-1)*4+i, ind, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, (j-1)*4+i, ind, MVN_staticswpind, 1])))
    ps = Point3f(0.0, 0.0, 0.0)
    pe = Point3f(dot(v1, MVN_perpstatic1), dot(v1, MVN_parastatic), dot(v1, MVN_perpstatic2))
    arrows3d!(ls1, ps, pe, argmode = :endpoint, color = :red, shaftradius = 0.02, tipradius = 0.05, tailradius = 0.02)

    #probably associated with the reflected protons by the ionosphere
    i, j =1, 11
    ind = 14
    #ind = argmax(MVN_staticefluxmean[1, (j-1)*4+i, :]) 
    v3 = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, ind, MVN_staticswpind, 1], deg2rad(MVN_staticphi[1, (j-1)*4+i, ind, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, (j-1)*4+i, ind, MVN_staticswpind, 1])))
    ps = Point3f(0.0, 0.0, 0.0)
    pe = Point3f(dot(v3, MVN_perpstatic1), dot(v3, MVN_parastatic), dot(v3, MVN_perpstatic2))
    arrows3d!(ls1, ps, pe, argmode = :endpoint, color = :pink, shaftradius = 0.02, tipradius = 0.05, tailradius = 0.02)

    i, j =2, 11
    ind = 15
    #ind = argmax(MVN_staticefluxmean[1, (j-1)*4+i, :]) 
    v3 = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, ind, MVN_staticswpind, 1], deg2rad(MVN_staticphi[1, (j-1)*4+i, ind, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, (j-1)*4+i, ind, MVN_staticswpind, 1])))
    ps = Point3f(0.0, 0.0, 0.0)
    pe = Point3f(dot(v3, MVN_perpstatic1), dot(v3, MVN_parastatic), dot(v3, MVN_perpstatic2))
    arrows3d!(ls1, ps, pe, argmode = :endpoint, color = :pink, shaftradius = 0.02, tipradius = 0.05, tailradius = 0.02)

    #probably associated with the reflected protons by the ionosphere
    i, j =3, 11
    ind = 15
    #ind = argmax(MVN_staticefluxmean[1, (j-1)*4+i, :]) 
    v3 = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, ind, MVN_staticswpind, 1], deg2rad(MVN_staticphi[1, (j-1)*4+i, ind, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, (j-1)*4+i, ind, MVN_staticswpind, 1])))
    ps = Point3f(0.0, 0.0, 0.0)
    pe = Point3f(dot(v3, MVN_perpstatic1), dot(v3, MVN_parastatic), dot(v3, MVN_perpstatic2))
    arrows3d!(ls1, ps, pe, argmode = :endpoint, color = :pink, shaftradius = 0.02, tipradius = 0.05, tailradius = 0.02)

    i, j =4, 11
    ind = 15
    #ind = argmax(MVN_staticefluxmean[1, (j-1)*4+i, :]) 
    v3 = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, ind, MVN_staticswpind, 1], deg2rad(MVN_staticphi[1, (j-1)*4+i, ind, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, (j-1)*4+i, ind, MVN_staticswpind, 1])))
    ps = Point3f(0.0, 0.0, 0.0)
    pe = Point3f(dot(v3, MVN_perpstatic1), dot(v3, MVN_parastatic), dot(v3, MVN_perpstatic2))
    arrows3d!(ls1, ps, pe, argmode = :endpoint, color = :pink, shaftradius = 0.02, tipradius = 0.05, tailradius = 0.02)

    #probably associated with the pickup protons
    i, j =2, 13
    ind = 9
    v2 = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, ind, MVN_staticswpind, 1], deg2rad(MVN_staticphi[1, (j-1)*4+i, ind, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, (j-1)*4+i, ind, MVN_staticswpind, 1])))
    ps = Point3f(0.0, 0.0, 0.0)
    pe = Point3f(dot(v2, MVN_perpstatic1), dot(v2, MVN_parastatic), dot(v2, MVN_perpstatic2))
    arrows3d!(ls1, ps, pe, argmode = :endpoint, color = :blue, shaftradius = 0.02, tipradius = 0.05, tailradius = 0.02)

    

    #probably associated with the pickup protons
    i, j =2, 1
    ind = argmax(MVN_staticefluxmean[1, (j-1)*4+i, :])
    v4 = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, ind, MVN_staticswpind, 1], deg2rad(MVN_staticphi[1, (j-1)*4+i, ind, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, (j-1)*4+i, ind, MVN_staticswpind, 1])))
    ps = Point3f(0.0, 0.0, 0.0)
    pe = Point3f(dot(v4, MVN_perpstatic1), dot(v4, MVN_parastatic), dot(v4, MVN_perpstatic2))
    arrows3d!(ls1, ps, pe, argmode = :endpoint, color = :purple, shaftradius = 0.02, tipradius = 0.05, tailradius = 0.02)

    i, j =2, 10
    ind = argmax(MVN_staticefluxmean[1, (j-1)*4+i, :])
    v5 = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, ind, MVN_staticswpind, 1], deg2rad(MVN_staticphi[1, (j-1)*4+i, ind, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, (j-1)*4+i, ind, MVN_staticswpind, 1])))
    ps = Point3f(0.0, 0.0, 0.0)
    pe = Point3f(dot(v5, MVN_perpstatic1), dot(v5, MVN_parastatic), dot(v5, MVN_perpstatic2))
    arrows3d!(ls1, ps, pe, argmode = :endpoint, color = :orange, shaftradius = 0.02, tipradius = 0.05, tailradius = 0.02)

    i, j =2, 16
    ind = argmax(MVN_staticefluxmean[1, (j-1)*4+i, :])
    v6 = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, ind, MVN_staticswpind, 1], deg2rad(MVN_staticphi[1, (j-1)*4+i, ind, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, (j-1)*4+i, ind, MVN_staticswpind, 1])))
    ps = Point3f(0.0, 0.0, 0.0)
    pe = Point3f(dot(v6, MVN_perpstatic1), dot(v6, MVN_parastatic), dot(v6, MVN_perpstatic2))
    arrows3d!(ls1, ps, pe, argmode = :endpoint, color = :green, shaftradius = 0.02, tipradius = 0.05, tailradius = 0.02)

    #viz!(ls, Meshes.Torus((0, vpara[115], dot(MVN_staticvswrough, MVN_perpstatic2), (0, 1, 0)), 300), transparency=true, color =:blue)

    r = MVN_staticvelocitygridH[end-11:end-8]
    θ = deg2rad.(MVN_staticthetagridH)
    ϕ = deg2rad.(MVN_staticphigridH[vcat(11:16,1:2)])

    #r = MVN_staticvelocitygridH[begin:end-6]
    #θ = deg2rad.(MVN_staticthetagridH)
    #ϕ = deg2rad.(MVN_staticphigridH)

    # grid topology with given number of elements per dimension
    # the last coordinate is periodic so we set (false, false, true)
    topo = GridTopology((length(r)-1, length(θ)-1, length(ϕ)-1), (false, false, true))
    coordrot = RotMatrix{3}(transpose(stack([MVN_perpstatic1, MVN_parastatic, MVN_perpstatic2])))
    g = RectilinearGrid{𝔼,typeof(Meshes.Spherical(0,0,0))}((r, θ, ϕ), topo) |> Rotate(coordrot)

    #viz!(ls, g, showsegments = true, transparency=true, alpha = 0.3)

    v7 = MVN_static2msomean' * [800.0, 0.0, 0.0] #xmso
    ps = Point3f(0.0, 0.0, 0.0)
    pe = Point3f(dot(v7, MVN_perpstatic1), dot(v7, MVN_parastatic), dot(v7, MVN_perpstatic2))
    arrows3d!(ls1, ps, pe, argmode = :endpoint, color = :black, shaftradius = 0.02, tipradius = 0.05, tailradius = 0.02)


    v8 = MVN_static2msomean' * [0.0, 0.0, 800.0] #xmso
    ps = Point3f(0.0, 0.0, 0.0)
    pe = Point3f(dot(v8, MVN_perpstatic1), dot(v8, MVN_parastatic), dot(v8, MVN_perpstatic2))
    arrows3d!(ls1, ps, pe, argmode = :endpoint, color = :red, shaftradius = 0.02, tipradius = 0.05, tailradius = 0.02)

    =#
    
    ls2 = LScene(fig[1, 2])
    #for ipara in 80:10:120 #length(vpara)
    #for ipara in 1:length(vpara)
    for ipara in 90:5:110
        vparai = vpara[ipara]
        heatmap!(ls2, vperp1, vperp2, MVN_swia_eflux_rect[:, ipara, :], transformation = (:xz, vparai), colormap = :gist_earth,  colorrange = (1e3, 1e9), colorscale = log10)
    end
    

    fig

    i, j =2, 13
    ind = 9
    MVN_staticenergy[1, (j-1)*4+i, ind, MVN_staticswpind, 1]
    

    #, transparency = true, alpha = 0.5
    fig

    MVN_static_eflux_rect = Array{Float64, 3}(undef, length(vperp1), length(vpara), length(vperp2))
    cn = cgrad(:gist_earth, 15, categorical = true)
    #cn = Makie.wong_colors()
    ax1 = Axis3(fig[1, 1])
    contour!(ax1, -800..800, -800..800, -800..800, log10.(MVN_static_eflux_rect), colorrange = (3, 9), colorscale = identity, colormap = :gist_earth)
    fig

    contour!(ax1, -800..800, -800..800, -800..800, MVN_static_eflux_rect, colorrange = (1e3, 1e9), colorscale = log10, alpha = 0.05)

    vperp1 = range(-800.0, 800.0, step = 50.0)
    vperp2 = range(-800.0, 800.0, step = 50.0)
    static_eflux = Matrix{Float64}(undef, length(vperp1), length(vperp2))
    for i in 1:length(vperp1)
        for j in 1:length(vperp2)
            v3 = Array(-250 .* MVN_parastatic .+ vperp1[i] .* MVN_perpstatic1 .+ vperp2[j] .* MVN_perpstatic2)
            static_eflux[i, j] = MVN_STATIC_eflux_interp(v3)
        end
    end

    size_inches = (9, 10)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 10, 9, 8))
    cn = cgrad(:gist_earth, 15, categorical = true)
    #cn = Makie.wong_colors()
    ax1 = Axis(fig[1, 1], aspect = DataAspect())
    heatmap!(ax1, vperp1, vperp2, static_eflux, colormap = :gist_earth,  colorrange = (1e3, 1e9), colorscale = log10)
    xlims!(ax1, -800, 800)
    ylims!(ax1, -800, 800)
    contour!(ax1, vperp1, vperp2, log10.(static_eflux), levels = range(3, 9, 50), linewidth = 0.5, color = :black)

    fig

    
    ax1 = Axis(fig[1, 1], xticks=(-800:400:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("Parallel Velocity (km s"), superscript("-1"), rich(")")), ylabel = rich(rich("MAVEN STATIC H"), superscript("+"), rich("\nDifferential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)

    for iv in 0:14
        MVN_staticvpara, MVN_staticv3para, MVN_staticefluxpara = MVN_STATIC_eflux_slice(dot(MVN_staticvswrough, MVN_perpstatic2) .* MVN_perpstatic2 + (iv*50) .* MVN_perpstatic2, MVN_parastatic)

        scatterlines!(ax1, MVN_staticvpara, MVN_staticefluxpara, color = cn[iv+1], linewidth = 2)
    end

    xlims!(ax1, -300, 750)
    ylims!(ax1, 1e3, 1e9)

    fig

    

    MVN_staticvpara, MVN_staticv3para, MVN_staticefluxpara = MVN_STATIC_eflux_slice(MVN_staticvswrough, MVN_parastatic)
    MVN_staticvperp1, MVN_staticv3perp1, MVN_staticefluxperp1 = MVN_STATIC_eflux_slice(MVN_staticvswrough, MVN_perpstatic1)
    MVN_staticvperp2, MVN_staticv3perp2, MVN_staticefluxperp2 = MVN_STATIC_eflux_slice(MVN_staticvswrough, MVN_perpstatic2)

    MVN_swiavpara, MVN_swiav3para, MVN_swiaefluxpara = MVN_SWIA_eflux_slice(MVN_swiavswrough, MVN_parastatic)
    MVN_swiavperp1, MVN_swiav3perp1, MVN_swiaefluxperp1 = MVN_SWIA_eflux_slice(MVN_swiavswrough, MVN_perpstatic1)
    MVN_swiavperp2, MVN_swiav3perp2, MVN_swiaefluxperp2 = MVN_SWIA_eflux_slice(MVN_swiavswrough, MVN_perpstatic2)

    # Define the Gaussian function：p = [A, μ, σ]
    gaussian(x, p) = @. p[1] * exp(-(x - p[2])^2 / (2 * p[3]^2))

    size_inches = (9, 10)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 10, 9, 8))
    cn = cgrad(:gist_earth, 18, categorical = true)
    cn = Makie.wong_colors()
    
    ax1 = Axis(fig[1, 1], xticks=(-800:400:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("Parallel Velocity (km s"), superscript("-1"), rich(")")), ylabel = rich(rich("MAVEN STATIC H"), superscript("+"), rich("\nDifferential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)
    
    ind = findall(0.0 .< MVN_staticvpara .< 300.0)
    x = MVN_staticvpara[ind]
    y_data = MVN_staticefluxpara[ind] ./ norm.(MVN_staticv3para[ind]).^4
    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    vpara = range(-800, 800, length=81)
    Vthpara = fit.param[3]
    Tpara = ustrip(u"eV", m_p* (Vthpara * 1u"km/s")^2)
    
    v = sqrt.(vpara.^2 .+ dot(MVN_staticvswrough, MVN_perpstatic2)^2)
    eflux = gaussian(vpara, fit.param) .* v.^4

    lines!(ax1, vpara, eflux, color = cn[2], linewidth = 2)
    scatter!(ax1, MVN_staticvpara, MVN_staticefluxpara, color = cn[1]) #, linewidth = 1, linestyle = :dash

    xlims!(ax1, -750, 750)
    ylims!(ax1, 1e3, 1e9)

    ax2 = Axis(fig[1, 2], xticks=(-800:400:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("First Perpendicular Velocity (km s"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)
    
    ind = findall(-300.0 .< MVN_staticvperp1 .< 300.0)
    x = MVN_staticvperp1[ind]
    y_data = MVN_staticefluxperp1[ind] ./ norm.(MVN_staticv3perp1[ind]).^4
    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)

    Vthperp = fit.param[3]
    Tperp = ustrip(u"eV", m_p* (Vthperp * 1u"km/s")^2)

    N = ustrip(u"cm^-3", fit.param[1]*1u"cm^-2*s^-1"/(1u"km/s")*2 * (2*π)^1.5 * Vthpara * Vthperp^2)

    fit.param[2] = 0.0
    vperp1 = range(-800, 800, length=81)

    v = norm(MVN_staticvswrough)
    eflux = gaussian(vperp1, fit.param) .* v^4

    static_paramperp = deepcopy(fit.param)

    lines!(ax2, vperp1, eflux, color = cn[2], linewidth = 2)
    scatter!(ax2, MVN_staticvperp1, MVN_staticefluxperp1, color = cn[1]) #, linewidth = 1, linestyle = :dash

    xlims!(ax2, -750, 750)
    ylims!(ax2, 1e3, 1e9)

    text!(ax1, 0, 1.0, text = rich(rich("N", font=:italic), subscript("SW"), rich("="*string(round(N, digits=1))), rich(" cm"), superscript("-3"), rich("\nT", font=:italic), subscript("SW∥"), rich("="*string(round(Tpara, digits=1))), rich(" eV"), rich("\nT", font=:italic), subscript("SW⟂"), rich("="*string(round(Tperp, digits=1))), rich(" eV")), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 15, color = :black)
    

    ax3 = Axis(fig[1, 3], xticks=(-800:400:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("Second Perpendicular Velocity (km s"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)

    
    fit.param[2] = dot(MVN_staticvswrough, MVN_perpstatic2)-20.0

    v = sqrt.(vperp2.^2 .+ dot(MVN_staticvswrough, MVN_parastatic)^2)
    eflux = gaussian(vperp2, fit.param) .* v.^4

    lines!(ax3, vperp2, eflux, color = cn[2], linewidth = 2)
    scatter!(ax3, MVN_staticvperp2, MVN_staticefluxperp2, color = cn[1])
    
    xlims!(ax3, -750, 750)
    ylims!(ax3, 1e3, 1e9)

    ax4 = Axis(fig[2, 1], xticks=(-800:400:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("Parallel\nVelocity (km s"), superscript("-1"), rich(")")), ylabel = rich(rich("MAVEN SWIA Ion"), rich("\nDifferential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)
    
    ind = findall(0.0 .< MVN_swiavpara .< 300.0)
    x = MVN_swiavpara[ind]
    y_data = MVN_swiaefluxpara[ind] ./ norm.(MVN_swiav3para[ind]).^4
    # estimation
    μ_est = sum(x .* y_data) / sum(y_data)
    σ2_est = sum((x .- μ_est).^2 .* y_data) / sum(y_data)
    A_est = maximum(y_data)
    p0 = [A_est, μ_est, sqrt(σ2_est)]
    fit = curve_fit(gaussian, x, y_data, p0)
    Vthpara = fit.param[3]
    Tpara = ustrip(u"eV", m_p* (Vthpara * 1u"km/s")^2)
    
    N = ustrip(u"cm^-3", fit.param[1]*1u"cm^-2*s^-1"/(1u"km/s")*2 * (2*π)^1.5 * Vthpara * Vthperp^2)

    @show Vthpara, Tpara, Vthperp, Tperp, N, fit.param
    
    vpara = range(-800, 800, length=81)

    v = sqrt.(vpara.^2 .+ dot(MVN_swiavswrough, MVN_perpswia2)^2)
    eflux = gaussian(vpara, fit.param) .* v.^4
    #fit.param[1] = maximum(MVN_swiaefluxpara[ind])/maximum(eflux)*fit.param[1]
    #eflux = gaussian(vpara, fit.param) .* v.^4
    maxeflux = maximum(eflux)

    lines!(ax4, vpara, eflux, color = cn[2], linewidth = 2)
    scatter!(ax4, MVN_swiavpara, MVN_swiaefluxpara, color = cn[1])

    xlims!(ax4, -750, 750)
    ylims!(ax4, 1e3, 1e9)


    ax5 = Axis(fig[2, 2], xticks=(-800:400:800), xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("1st Perpendicular\nVelocity (km s"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)

    vperp1 = range(-800, 800, length=81)

    v = norm(MVN_swiavswrough)
    param = static_paramperp
    eflux = gaussian(vperp1, param) .* v^4
    param[1] = maxeflux/maximum(eflux)*param[1]
    param[3] = param[3]*0.85
    eflux = gaussian(vperp1, param) .* v^4
    
    lines!(ax5, vperp1, eflux, color = cn[2], linewidth = 2)
    scatter!(ax5, MVN_swiavperp1, MVN_swiaefluxperp1, color = cn[1])

    Vthperp = param[3]
    Tperp = ustrip(u"eV", m_p* (Vthperp * 1u"km/s")^2)
    N = ustrip(u"cm^-3", param[1]*1u"cm^-2*s^-1"/(1u"km/s")*2 * (2*π)^1.5 * Vthpara * Vthperp^2)

    @show Vthpara, Tpara, Vthperp, Tperp, N, param
    
    xlims!(ax5, -750, 750)
    ylims!(ax5, 1e3, 1e9)

    text!(ax4, 0, 1.0, text = rich(rich("N", font=:italic), subscript("SW"), rich("="*string(round(N, digits=1))), rich(" cm"), superscript("-3"), rich("\nT", font=:italic), subscript("SW∥"), rich("="*string(round(Tpara, digits=1))), rich(" eV"), rich("\nT", font=:italic), subscript("SW⟂"), rich("="*string(round(Tperp, digits=1))), rich(" eV")), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 15, color = :black)
    

    ax6 = Axis(fig[2, 3], xticks=(-800:400:800),xminorticksvisible=true, xminorticks = IntervalsBetween(4), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich(rich("2nd Perpendicular\nVelocity (km s"), superscript("-1"), rich(")")), xscale = identity, yscale=log10, backgroundcolor=:gray80)

    lines!(ax6, vperp1.-20.0, eflux, color = cn[2], linewidth = 2)
    scatter!(ax6, MVN_swiavperp2, MVN_swiaefluxperp2, color = cn[1])
    
    xlims!(ax6, -750, 750)
    ylims!(ax6, 1e3, 1e9)

    ind = [1, 2]
    txt = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"]
    for (i, indi) in enumerate(ind)
        Label(fig[indi, 0], txt[i],font = :bold, valign=:top, tellheight=false, fontsize = 20)
    end
    colsize!(fig.layout, 0, Relative(1/20))
    colgap!(fig.layout, 1, -30)

    axu = [ax1, ax2, ax3]
    linkyaxes!(axu...)
    hidexdecorations!.(axu[1:end], grid=false)
    hideydecorations!.(axu[2:end], grid=false)
    axb = [ax4, ax5, ax6]
    linkxaxes!(ax1, ax4)
    linkxaxes!(ax2, ax5)
    linkxaxes!(ax3, ax6)
    linkyaxes!(axb...)
    hideydecorations!.(axb[2:end], grid=false)
    fig
end

function MVN_STATIC_efluxBN_plot()
    size_inches = (9, 10)
    size_pt = 72 .* size_inches
    figyz = Figure(size = size_pt, fontsize = 15, figure_padding = (0, 0, 2, 0))
    
    ax1 = PolarAxis(figyz[1, 1], backgroundcolor=:gray80, rticks = (0:200:1000), gridz=1)
    vp = voronoiplot!(ax1, MVN_staticvarphi, MVN_staticvtot, MVN_staticefluxyz, show_generators = false, strokewidth = 0, colormap=:gist_earth, colorscale=log10, colorrange=(1e3,1e9))
    arrows2d!([0.0], [0.0], [angle(MVN_staticvswrough, MVN_parastatic)], [norm(MVN_staticvswrough)], argmode = :endpoint, color = :white, tiplength = 14)
    arrows2d!([0.0], [0.0], [0.0], [850], argmode = :endpoint, color = :white, tiplength = 14)
    arrows2d!([0.0], [0.0], [π/2], [850], argmode = :endpoint, color = :white, tiplength = 14)
    
    rlims!(ax1, 0.0, 850)
    hiderdecorations!(ax1, grid = false)

    for ig in 200:200:800
        text!(ax1, deg2rad(135), ig+45, text = rich(string(ig), rich(" km s"), superscript("-1")), align = (:center, :top), offset = (2, 5), space = :data, fontsize = 15, color = :black, rotation = π/4)
    end
    text!(ax1, 2.0*angle(MVN_staticvswrough, MVN_parastatic), 0.5*norm(MVN_staticvswrough), text = rich(rich("V", font=:bold), subscript("SW")), align = (:left, :top), offset = (0, 6), space = :data, fontsize = 15, color = :white, rotation = angle(MVN_staticvswrough, MVN_parastatic))
    text!(ax1, 0.0, 850.0, text = rich(rich("B", font=:bold), subscript("0")), align = (:right, :bottom), offset = (-16, 2), space = :data, fontsize = 15, color = :white)
    text!(ax1, π/2, 850.0, text = rich(rich("N", font=:bold)), align = (:left, :top), offset = (6, -12), space = :data, fontsize = 15, color = :white)

    Colorbar(figyz[2, 1], vp, vertical = false, flipaxis = false, label = rich(rich("MAVEN STATIC H"), superscript("+"), rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), tellwidth = false, width = Relative(7/9))
    return figyz
end

function MVN_STATIC_efluxEB_plot()
    size_inches = (9, 10)
    size_pt = 72 .* size_inches
    figxy = Figure(size = size_pt, fontsize = 15, figure_padding = (0, 0, 2, 0))

    ax2 = PolarAxis(figxy[1, 1], backgroundcolor=:gray80, rticks = (0:200:1000), gridz=1)
    vp = voronoiplot!(ax2, MVN_staticvarphi, MVN_staticvtot, MVN_staticefluxxy, show_generators = false, strokewidth = 0, colormap=:gist_earth, colorscale=log10, colorrange=(1e3,1e9))
    arrows2d!([0.0], [0.0], [0.0], [850], argmode = :endpoint, color = :white, tiplength = 14)
    arrows2d!([0.0], [0.0], [π/2], [850], argmode = :endpoint, color = :white, tiplength = 14)
    rlims!(ax2, 0.0, 850)
    hiderdecorations!(ax2, grid = false)

    for ig in 200:200:800
        text!(ax2, deg2rad(135), ig+45, text = rich(string(ig), rich(" km s"), superscript("-1")), align = (:center, :top), offset = (2, 5), space = :data, fontsize = 15, color = :black, rotation = π/4)
    end
    text!(ax2, 0.0, 850.0, text = rich(rich("E", font=:bold), subscript("SW")), align = (:right, :bottom), offset = (-16, 2), space = :data, fontsize = 15, color = :white)
    text!(ax2, π/2, 850.0, text = rich(rich("B", font=:bold), subscript("0")), align = (:left, :top), offset = (6, -12), space = :data, fontsize = 15, color = :white)
    Colorbar(figxy[2, 1], vp, vertical = false, flipaxis = false, label = rich(rich("MAVEN STATIC H"), superscript("+"), rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), tellwidth = false, width = Relative(7/9))
    return figxy
end

function MVN_STATIC_efluxEN_plot()   
    size_inches = (9, 10)
    size_pt = 72 .* size_inches
    figxz = Figure(size = size_pt, fontsize = 15, figure_padding = (0, 0, 2, 0))

    ax3 = PolarAxis(figxz[1, 1], backgroundcolor=:gray80, rticks = (0:200:1000), gridz=1)
    vp = voronoiplot!(ax3, MVN_staticvarphi, MVN_staticvtot, MVN_staticefluxxz, show_generators = false, strokewidth = 0, colormap=:gist_earth, colorscale=log10, colorrange=(1e3,1e9))
    arrows2d!([0.0], [0.0], [0.0], [850], argmode = :endpoint, color = :white, tiplength = 14)
    arrows2d!([0.0], [0.0], [π/2], [850], argmode = :endpoint, color = :white, tiplength = 14)
    rlims!(ax3, 0.0, 850)
    hiderdecorations!(ax3, grid = false)

    for ig in 200:200:800
        text!(ax3, deg2rad(45), ig+30, text = rich(string(ig), rich(" km s"), superscript("-1")), align = (:center, :top), offset = (8, 1), space = :data, fontsize = 15, color = :black, rotation = -π/4)
    end
    text!(ax3, 0.0, 850.0, text = rich(rich("E", font=:bold), subscript("SW")), align = (:right, :bottom), offset = (-16, 2), space = :data, fontsize = 15, color = :white)
    text!(ax3, π/2, 850.0, text = rich(rich("N", font=:bold)), align = (:left, :top), offset = (6, -12), space = :data, fontsize = 15, color = :white)
    Colorbar(figxz[2, 1], vp, vertical = false, flipaxis = false, label = rich(rich("MAVEN STATIC H"), superscript("+"), rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), tellwidth = false, width = Relative(7/9))
    return figxz  
end

function MVN_STATIC_efluxmean_plot()
    size_inches = (9, 18)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 20, 9, 8))
    cn = cgrad(:gist_earth, 18, categorical = true)

    axs = [Axis(fig[i, 1], xminorticksvisible=true, xminorticks = IntervalsBetween(9), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich("Energy (eV)"), ylabel = rich(rich("H"), superscript("+"), rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), xscale = log10, yscale=log10, backgroundcolor=:gray80) for i in 1:4]
    for i in 1:4
        for j in 1:16
            ind = findall(MVN_staticefluxmean[1, (j-1)*4+i, :] .> 0)
            scatterlines!(axs[i], MVN_staticenergy[1, (j-1)*4+i, ind, MVN_staticswpind, 1], MVN_staticefluxmean[1, (j-1)*4+i, ind], color = cn[j+1], linewidth = 3)
        end
        text!(axs[i], 0.0, 1.0, text = "MAVEN STATIC θ=$i", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
        ax1_inset = Axis(fig[i, 1], width=Relative(0.8), height=Relative(0.05), halign=0.5, valign=0.02)
        phi = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], (16, 1))
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        y = [1]
        heatmap!(ax1_inset, x, y, phi, colormap=cn, colorrange=(0,17))
        hidedecorations!(ax1_inset, grid=false)
        for ix in x
            text!(ax1_inset, ix, 1.0, text = "$ix", align = (:center, :center), offset = (0, 0), fontsize = 20, color = :black)
        end
        text!(axs[i], 0.042, 0.07, text = "φ=", align = (:left, :top), offset = (0, 0), space = :relative, fontsize = 20, color = :black)
    end
    hidexdecorations!.(axs[begin:end-1], grid=false)
    linkxaxes!(axs)
    xlims!.(axs, 1e0, 1e4)
    ylims!.(axs, 1e3, 1e9)

    fig
end



function MVN_SWIA_eflux_Overview_plot()
    date = DateTime(2022, 06, 12, 00, 00, 00)
    size_inches = (20, 30)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 20, figure_padding = (5,20, 9, 10))
    cn = cgrad(:seaborn_dark, 10, categorical = true)
    lcn = Makie.wong_colors()

    xd = date .+ Dates.Hour(9) .+ Dates.Minute(30) .+ Dates.Minute(10) .* range(0, 5)
    xtk = datetime2julian.(xd)
    xtk = Float32.(xtk .- MVN_swiajulian[1])

    axes = [Axis(fig[j, i], yscale=log10, xticks = (xtk, Dates.format.(xd, "HH:MM")), backgroundcolor=:gray80, xticklabelrotation=-π/3) for j in 1:16, i in 1:4]

    for j in 1:16
        for i in 1:4
            ax = axes[j, i]
            heatmap!(ax, Float32.(MVN_swiajulian.-MVN_swiajulian[1]), MVN_swiaenergy[:, 1], transpose(MVN_swiaeflux[j, i, :, :]), colormap=:gist_earth, colorscale=log10, colorrange = (1e5, 1e9))
            text!(ax, 0.0, 1.0, text = rich("θ=$i ϕ=$j"), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 20, color = :white)
            xlims!(ax, xtk[begin], xtk[end]) 
            ylims!(ax, 1e1, 1e4)
            
            if i==4 && (j==14 || j==15)
                annotation!(ax, -50, -40, datetime2julian(DateTime(2022, 06, 12, 10, 20, 00))-MVN_swiajulian[1], 300, path = Ann.Paths.Corner(), text = "SW", color=cn[2], fontsize=20)
            end
            if i==2 && (j==15 || j==16)
                annotation!(ax, -130, -77.5, datetime2julian(DateTime(2022, 06, 12, 10, 20, 00))-MVN_swiajulian[1], 3e3, path = Ann.Paths.Corner(), text = "SPU", color=cn[2], fontsize=20)
                annotation!(ax, -90, -58, datetime2julian(DateTime(2022, 06, 12, 10, 20, 00))-MVN_swiajulian[1], 900, path = Ann.Paths.Corner(), text = "PU", color=cn[2], fontsize=20)
            end
            
        end
    end
    linkaxes!(axes...)
    hidexdecorations!.(axes[1:15,:], grid=false)
    hideydecorations!.(axes[:, 2:end], grid=false)
    colgap!(fig.layout, 30)
    Colorbar(fig[1:16, 5], colormap=:gist_earth, scale=log10, colorrange = (1e5, 1e9), label = rich(rich("MAVEN SWIA Ion Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")))
    Label(fig[1:16, 0], rich("Energy (eV)"), fontsize = 25, halign = :center, rotation = π/2)
    Label(fig[17, 1:5], rich("Universal Time"), fontsize = 20, halign = :center)

    fig
end

function MVN_SWIA_efluxBN_plot()
    size_inches = (9, 10)
    size_pt = 72 .* size_inches
    figyz = Figure(size = size_pt, fontsize = 15, figure_padding = (0, 0, 2, 0))
    
    ax1 = PolarAxis(figyz[1, 1], backgroundcolor=:gray80, rticks = (0:200:1000), gridz=1)
    vp = voronoiplot!(ax1, MVN_swiavarphi, MVN_swiavtot, MVN_swiaefluxyz, show_generators = false, strokewidth = 0, colormap=:gist_earth, colorscale=log10, colorrange=(1e3,1e9))
    arrows2d!([0.0], [0.0], [angle(MVN_swiavswrough, MVN_paraswia)], [norm(MVN_swiavswrough)], argmode = :endpoint, color = :white, tiplength = 14)
    arrows2d!([0.0], [0.0], [0.0], [850], argmode = :endpoint, color = :white, tiplength = 14)
    arrows2d!([0.0], [0.0], [π/2], [850], argmode = :endpoint, color = :white, tiplength = 14)
    
    rlims!(ax1, 0.0, 850)
    hiderdecorations!(ax1, grid = false)

    for ig in 200:200:800
        text!(ax1, deg2rad(180), ig+25, text = rich(string(ig), rich(" km s"), superscript("-1")), align = (:center, :center), offset = (2, 5), space = :data, fontsize = 15, color = :black, rotation = π/2)
    end
    text!(ax1, 2.0*angle(MVN_swiavswrough, MVN_paraswia), 0.5*norm(MVN_swiavswrough), text = rich(rich("V", font=:bold), subscript("SW")), align = (:left, :top), offset = (0, 6), space = :data, fontsize = 15, color = :white, rotation = angle(MVN_swiavswrough, MVN_paraswia))
    text!(ax1, 0.0, 850.0, text = rich(rich("B", font=:bold), subscript("0")), align = (:right, :bottom), offset = (-16, 2), space = :data, fontsize = 15, color = :white)
    text!(ax1, π/2, 850.0, text = rich(rich("N", font=:bold)), align = (:left, :top), offset = (6, -12), space = :data, fontsize = 15, color = :white)

    Colorbar(figyz[2, 1], vp, vertical = false, flipaxis = false, label = rich(rich("MAVEN SWIA Ion"), rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), tellwidth = false, width = Relative(7/9))
    return figyz
end

function MVN_SWIA_efluxEB_plot()
    size_inches = (9, 10)
    size_pt = 72 .* size_inches
    figxy = Figure(size = size_pt, fontsize = 15, figure_padding = (0, 0, 2, 0))

    ax2 = PolarAxis(figxy[1, 1], backgroundcolor=:gray80, rticks = (0:200:1000), gridz=1)
    vp = voronoiplot!(ax2, MVN_swiavarphi, MVN_swiavtot, MVN_swiaefluxxy, show_generators = false, strokewidth = 0, colormap=:gist_earth, colorscale=log10, colorrange=(1e3,1e9))
    arrows2d!([0.0], [0.0], [0.0], [850], argmode = :endpoint, color = :white, tiplength = 14)
    arrows2d!([0.0], [0.0], [π/2], [850], argmode = :endpoint, color = :white, tiplength = 14)
    rlims!(ax2, 0.0, 850)
    hiderdecorations!(ax2, grid = false)

    for ig in 200:200:800
        text!(ax2, deg2rad(45), ig, text = rich(string(ig), rich(" km s"), superscript("-1")), align = (:center, :top), offset = (16, 8), space = :data, fontsize = 15, color = :black, rotation = -π/4)
    end
    text!(ax2, 0.0, 850.0, text = rich(rich("E", font=:bold), subscript("SW")), align = (:right, :bottom), offset = (-16, 2), space = :data, fontsize = 15, color = :white)
    text!(ax2, π/2, 850.0, text = rich(rich("B", font=:bold), subscript("0")), align = (:left, :top), offset = (6, -12), space = :data, fontsize = 15, color = :white)
    Colorbar(figxy[2, 1], vp, vertical = false, flipaxis = false, label = rich(rich("MAVEN SWIA Ion"), rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), tellwidth = false, width = Relative(7/9))
    return figxy
end

function MVN_SWIA_efluxEN_plot()   
    size_inches = (9, 10)
    size_pt = 72 .* size_inches
    figxz = Figure(size = size_pt, fontsize = 15, figure_padding = (0, 0, 2, 0))

    ax3 = PolarAxis(figxz[1, 1], backgroundcolor=:gray80, rticks = (0:200:1000), gridz=1)
    vp = voronoiplot!(ax3, MVN_swiavarphi, MVN_swiavtot, MVN_swiaefluxxz, show_generators = false, strokewidth = 0, colormap=:gist_earth, colorscale=log10, colorrange=(1e3,1e9))
    arrows2d!([0.0], [0.0], [0.0], [850], argmode = :endpoint, color = :white, tiplength = 14)
    arrows2d!([0.0], [0.0], [π/2], [850], argmode = :endpoint, color = :white, tiplength = 14)
    rlims!(ax3, 0.0, 850)
    hiderdecorations!(ax3, grid = false)

    for ig in 200:200:800
        text!(ax3, deg2rad(-45), ig, text = rich(string(ig), rich(" km s"), superscript("-1")), align = (:center, :top), offset = (2, 5), space = :data, fontsize = 15, color = :black, rotation = π/4)
    end
    text!(ax3, 0.0, 850.0, text = rich(rich("E", font=:bold), subscript("SW")), align = (:right, :bottom), offset = (-16, 2), space = :data, fontsize = 15, color = :white)
    text!(ax3, π/2, 850.0, text = rich(rich("N", font=:bold)), align = (:left, :top), offset = (6, -12), space = :data, fontsize = 15, color = :white)
    Colorbar(figxz[2, 1], vp, vertical = false, flipaxis = false, label = rich(rich("MAVEN SWIA Ion"), rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), tellwidth = false, width = Relative(7/9))

    return figxz  
end

function MVN_SWIA_efluxmean_plot()
    size_inches = (9, 18)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 20, 9, 8))
    
    cn = cgrad(:gist_earth, 18, categorical = true)

    axs = [Axis(fig[i, 1], xminorticksvisible=true, xminorticks = IntervalsBetween(9), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich("Energy (eV)"), ylabel = rich(rich("Ion Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), xscale = log10, yscale=log10, backgroundcolor=:gray80) for i in 1:4]
    for i in 1:4
        for j in 1:16
            ind = findall(MVN_swiaefluxmean[j, i, :] .> 0)
            scatterlines!(axs[i], MVN_swiaenergy[ind, 1], MVN_swiaefluxmean[j, i, ind], color = cn[j+1], linewidth = 3)
            #scatterlines!(axs[i], MVN_swiaenergy[ind, 1], MVN_swiaefluxmid[j, i, ind], color = cn[j+1], linewidth = 3)
        end
        text!(axs[i], 0.0, 1.0, text = "MAVEN SWIA θ=$i", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
        ax1_inset = Axis(fig[i, 1], width=Relative(0.8), height=Relative(0.05), halign=0.5, valign=0.02)
        phi = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], (16, 1))
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        y = [1]
        heatmap!(ax1_inset, x, y, phi, colormap=cn, colorrange=(0,17))
        hidedecorations!(ax1_inset, grid=false)
        for ix in x
            text!(ax1_inset, ix, 1.0, text = "$ix", align = (:center, :center), offset = (0, 0), fontsize = 20, color = :black)
        end
        text!(axs[i], 0.042, 0.07, text = "φ=", align = (:left, :top), offset = (0, 0), space = :relative, fontsize = 20, color = :black)
    end
    
    xlims!.(axs, 1e0, 1e4)
    ylims!.(axs, 1e3, 1e9)
    linkaxes!(axs...)
    hidexdecorations!.(axs[begin:end-1], grid=false)

    fig
end

function MVN_STATIC_SWIA_efluxmean_plot()
    size_inches = (9, 12)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 20, 9, 8))
    cn = cgrad(:gist_earth, 18, categorical = true)

    ax1 = Axis(fig[1, 1], xminorticksvisible=true, xminorticks = IntervalsBetween(9), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich("Energy (eV)"), ylabel = rich(rich("H"), superscript("+"), rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), xscale = log10, yscale=log10, backgroundcolor=:gray80)
    for j in 1:16
        ind = findall(MVN_staticefluxmean[1, (j-1)*4+2, :] .> 0)
        scatterlines!(ax1, MVN_staticenergy[1, (j-1)*4+2, ind, MVN_staticswpind, 1], MVN_staticefluxmean[1, (j-1)*4+2, ind], color = cn[j+1], linewidth = 3)
    end
    text!(ax1, 0.0, 1.0, text = "MAVEN STATIC θ=2", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    xlims!(ax1, 1e0, 1e4)
    ylims!(ax1, 1e3, 1e9)

    ax1_inset = Axis(fig[1, 1], width=Relative(0.8), height=Relative(0.05), halign=0.5, valign=0.02)
    phi = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], (16, 1))
    x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    y = [1]
    heatmap!(ax1_inset, x, y, phi, colormap=cn, colorrange=(0,17))
    hidedecorations!(ax1_inset, grid=false)
    for ix in x
        text!(ax1_inset, ix, 1.0, text = "$ix", align = (:center, :center), offset = (0, 0), fontsize = 20, color = :black)
    end
    text!(ax1, 0.042, 0.07, text = "φ=", align = (:left, :top), offset = (0, 0), space = :relative, fontsize = 20, color = :black)

    cn = cgrad(:gist_earth, 6, categorical = true)
    ax2 = Axis(fig[2, 1], xminorticksvisible=true, xminorticks = IntervalsBetween(9), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich("Energy (eV)"), ylabel = rich(rich("Ion Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), xscale = log10, yscale=log10, backgroundcolor=:gray80)
    for i in 1:4
        ind = findall(MVN_swiaefluxmean[15, i, :] .> 0)
        scatterlines!(ax2, MVN_swiaenergy[:, 1], MVN_swiaefluxmean[15, i, ind], color = cn[i+1], linewidth = 3)
    end
    text!(ax2, 0.0, 1.0, text = "MAVEN SWIA φ=15", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
    xlims!(ax2, 1e0, 1e4)
    ylims!(ax2, 1e3, 1e9)

    ax2_inset = Axis(fig[2, 1], width=Relative(0.8), height=Relative(0.05), halign=0.5, valign=0.02)
    θ = reshape([1, 2, 3, 4], (4, 1))
    x = [1, 2, 3, 4]
    y = [1]
    heatmap!(ax2_inset, x, y, θ, colormap=cn, colorrange=(0,5))
    hidedecorations!(ax2_inset, grid=false)
    for ix in x
        text!(ax2_inset, ix, 1.0, text = "$ix", align = (:center, :center), offset = (0, 0), fontsize = 20, color = :black)
    end
    text!(ax2, 0.05, 0.07, text = "θ=", align = (:left, :top), offset = (0, 0), space = :relative, fontsize = 20, color = :black)

    linkaxes!(ax1, ax2)
    hidexdecorations!(ax1, grid=false)

    ind = [1, 2]
    txt = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"]
    for (i, indi) in enumerate(ind)
        Label(fig[indi, 0], txt[i],font = :bold, valign=:top, tellheight=false, tellwidth=false, fontsize = 20)
    end
    colsize!(fig.layout, 0, Relative(1/40))

    fig
end


function TW_MINPA_eflux_Overview_plot(im)
    date = DateTime(2022, 06, 12, 00, 00, 00)
    size_inches = (20, 30)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 20, figure_padding = (5,20, 9, 10))
    cn = cgrad(:seaborn_dark, 10, categorical = true)
    lcn = Makie.wong_colors()

    xd = date .+ Dates.Hour(9) .+ Dates.Minute(30) .+ Dates.Minute(10) .* range(0, 5)
    xtk = datetime2julian.(xd)
    xtk = Float32.(xtk .- minpajulianmod1[1])

    axes = [Axis(fig[j, i], yscale=log10, xticks = (xtk, Dates.format.(xd, "HH:MM")), backgroundcolor=:gray80, xticklabelrotation=-π/3) for j in 1:16, i in 1:4]

    for j in 1:16
        for i in 1:4
            ax = axes[j, i]
            heatmap!(ax, Float32.(minpajulianmod1.-minpajulianmod1[1]), energymod1, TW_minpaefluxmod1[:,im,j,i,:], colormap=:gist_earth, colorscale=log10, colorrange = (1e5, 1e9))
            text!(ax, 0.0, 1.0, text = rich("θ=$i ϕ=$j"), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 20, color = :white)
            xlims!(ax, xtk[begin], xtk[end]) 
            ylims!(ax, 1e1, 1e4)
            
            if im == 1 && i==1 && j==12
                annotation!(ax, -80, -40, datetime2julian(DateTime(2022, 06, 12, 10, 00, 00))-minpajulianmod1[1], 600, path = Ann.Paths.Corner(), text = "SW", color=cn[2], fontsize=20)
                annotation!(ax, -120, -60, datetime2julian(DateTime(2022, 06, 12, 10, 00, 00))-minpajulianmod1[1], 2000, path = Ann.Paths.Corner(), text = "PU", color=cn[2], fontsize=20)
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

    fig
end

function TW_MINPA_efluxmean_plot()
    size_inches = (9, 18)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 20, 9, 8))
    
    cn = cgrad(:gist_earth, 18, categorical = true)

    axs = [Axis(fig[i, 1], xminorticksvisible=true, xminorticks = IntervalsBetween(9), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich("Energy (eV)"), ylabel = rich(rich("H"), superscript("+"), rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), xscale = log10, yscale=log10, backgroundcolor=:gray80) for i in 1:4]
    for i in 1:4
        for j in 1:16
            ind = findall(TW_minpaefluxmean[1, j, i, :] .> 0)
            scatterlines!(axs[i], energymod1[ind], TW_minpaefluxmean[1, j, i, ind], color = cn[j+1], linewidth = 3)
            #scatterlines!(axs[i], MVN_swiaenergy[ind, 1], MVN_swiaefluxmid[j, i, ind], color = cn[j+1], linewidth = 3)
        end
        text!(axs[i], 0.0, 1.0, text = "Tianwen-1 MINPA θ=$i", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
        ax1_inset = Axis(fig[i, 1], width=Relative(0.8), height=Relative(0.05), halign=0.5, valign=0.02)
        phi = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], (16, 1))
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        y = [1]
        heatmap!(ax1_inset, x, y, phi, colormap=cn, colorrange=(0,17))
        hidedecorations!(ax1_inset, grid=false)
        for ix in x
            text!(ax1_inset, ix, 1.0, text = "$ix", align = (:center, :center), offset = (0, 0), fontsize = 20, color = :black)
        end
        text!(axs[i], 0.042, 0.07, text = "φ=", align = (:left, :top), offset = (0, 0), space = :relative, fontsize = 20, color = :black)
    end
    
    xlims!.(axs, 1e0, 1e4)
    ylims!.(axs, 1e3, 1e9)
    linkaxes!(axs...)
    hidexdecorations!.(axs[begin:end-1], grid=false)

    fig
end

function NGITemperature_plot()
    size_inches = (9, 10)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5,25, 3, 8))
    cn = cgrad(:seaborn_dark, 10, categorical = true)
    lcn = Makie.wong_colors()

    date=DateTime(2022, 05, 25, 00, 00, 00)
    xd = date .+ Day(1).*range(0, 30, step=5)
    #xd = date .+ Hour(1).*[-15, -13]
    #xd = DateTime(2022, 6, 11, 10, 0, 0) .+ Minute(1).*range(-30, 30, step=15)
    xtk = datetime2julian.(xd)
    xtk = Float32.(xtk .- ngi_julian[1])

    ax1 = Axis(fig[1, 1], backgroundcolor=:gray80, xminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticks=IntervalsBetween(5), yminorticksvisible = true) #xticks=xtk, yminorticks=IntervalsBetween(5), yminorticksvisible = true, 
    #xlims!(ax1, 200, 5e3)
    ylims!(ax1, 1, 390)
    ax1.xticks = (xtk, Dates.format.(xd, "u dd"))
    xlims!(ax1, minimum(xtk), maximum(xtk))

    barplot!(ngi_julian[begin:2:end-1].- ngi_julian[1], ngi[begin:2:end-1, :Temperature], color = lcn[2])

    ax2 = Axis(fig[2, 1], backgroundcolor=:gray80, xminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticks=IntervalsBetween(5), yminorticksvisible = true) #xticks=xtk, yminorticks=IntervalsBetween(5), yminorticksvisible = true, 
    #xlims!(ax1, 200, 5e3)
    ylims!(ax2, 1, 390)
    ax2.xticks = (xtk, Dates.format.(xd, "u dd"))
    xlims!(ax2, minimum(xtk), maximum(xtk))
    barplot!(ngi_julian[begin:2:end-1].- ngi_julian[1], ngi[begin:2:end-1, :Temperature_error], color = lcn[2])

    ax3 = Axis(fig[3, 1], backgroundcolor=:gray80, xminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticks=IntervalsBetween(5), yminorticksvisible = true) #xticks=xtk, yminorticks=IntervalsBetween(5), yminorticksvisible = true, 
    #xlims!(ax1, 200, 5e3)
    ylims!(ax3, 1, 390)
    ax3.xticks = (xtk, Dates.format.(xd, "u dd"))
    xlims!(ax3, minimum(xtk), maximum(xtk))

    barplot!(ngi_julian[begin+1:2:end].- ngi_julian[1], ngi[begin+1:2:end, :Temperature], color = lcn[4])

    ax4 = Axis(fig[4, 1], xlabel = "Universal Time",  backgroundcolor=:gray80, xminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticks=IntervalsBetween(5), yminorticksvisible = true) #xticks=xtk, yminorticks=IntervalsBetween(5), yminorticksvisible = true, 
    #xlims!(ax1, 200, 5e3)
    ylims!(ax4, 1, 190)
    ax4.xticks = (xtk, Dates.format.(xd, "u dd"))
    xlims!(ax4, minimum(xtk), maximum(xtk))
    barplot!(ngi_julian[begin+1:2:end].- ngi_julian[1], ngi[begin+1:2:end, :Temperature_error], color = lcn[4])

    ind = [1, 2, 3, 4]
    txt = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"]
    for (i, indi) in enumerate(ind)
        Label(fig[indi, -1], txt[i],font = :bold, valign=:top, tellheight=false, fontsize = 20)
    end
    Label(fig[1, 0], rich("Ar Temperature" *"\n" * "(K)"), rotation = pi/2, tellheight=false)
    Label(fig[2, 0], rich("Ar Temperature Error" *"\n" * "(K)"), rotation = pi/2, tellheight=false)
    Label(fig[3, 0], rich("CO", subscript("2"), " Temperature" *"\n" * "(K)"), rotation = pi/2, tellheight=false)
    Label(fig[4, 0], rich("CO", subscript("2"), " Temperature Error" *"\n" * "(K)"), rotation = pi/2, tellheight=false)


    ax = [ax1, ax2, ax3, ax4]
    linkxaxes!(ax...)
    hidexdecorations!.(ax[1:end-1], grid=false)
    rowgap!(fig.layout, 3)
    colgap!(fig.layout, 2)
    fig
end

function MVN_TW_Overview_plot()
    size_inches = (9, 12)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5,0, 1, 2))
    cn = cgrad(:seaborn_dark, 10, categorical = true)
    lcn = Makie.wong_colors()
    xd = date .+ Dates.Hour(9) .+ Dates.Minute(0) .+ Dates.Minute(10) .* range(0, 9)
    
    xtk = datetime2julian.(xd)
    xtk = Float32.(xtk .- MVN_kp[1, :JulUT])

    ax1 = Axis(fig[1, 1], xticks=xtk, yminorticks=IntervalsBetween(5), yminorticksvisible = true, backgroundcolor=:gray80)
    ax1.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(minimum(xtk), maximum(xtk))
    ylims!(-330, -160)
    

    lines!(MVN_kp[!, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[!, 43], label = rich(rich("V", font = :italic), subscript("x")), color=lcn[3])
    lines!(minpajulianmod1.-MVN_kp[1, :JulUT], VH[:, 1], label = rich(rich("V", font = :italic), subscript("x")), color=lcn[1])

    ax2 = Axis(fig[2, 1],   xticks=xtk, yticks=range(-500, 500, step=50), yminorticks=IntervalsBetween(5), yminorticksvisible = true, backgroundcolor=:gray80, ylabelpadding=10.0)
    ax2.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(minimum(xtk), maximum(xtk))
    ylims!(-95, 45)

    lines!(MVN_kp[!, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[!, 45], label = rich(rich("V", font = :italic), subscript("y")), color=lcn[3])
    #lines!(ax12, minpajulianmod1.-MVN_kp[1, :JulUT], VH[:, 2], label = rich(rich("V", font = :italic), subscript("y")), color=lcn[1])

    ax3 = Axis(fig[3, 1],  xticks=xtk, yminorticks=IntervalsBetween(5), yminorticksvisible = true, backgroundcolor=:gray80)
    ax3.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(minimum(xtk), maximum(xtk))
    ylims!(-130, 70)

    lines!(MVN_kp[!, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[!, 47], label = rich(rich("V", font = :italic), subscript("z")), color=lcn[3])
    lines!(minpajulianmod1.-MVN_kp[1, :JulUT], VH[:, 3], label = rich(rich("V", font = :italic), subscript("z")), color=lcn[1])

 
    
    ax4=Axis(fig[4,1], xlabel = "Universal Time", xticks=xtk, yticks=range(0, 20, step=5), yminorticks=IntervalsBetween(5), yminorticksvisible = true, backgroundcolor=:gray80, ylabelpadding=10.0)
    ax4.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(minimum(xtk), maximum(xtk)) 
    ylims!(0.1, 13)
    #ax2.yreversed = true
    #swia
    lines!(MVN_kp[!, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[!, 41], label=rich(rich("N", font = :italic), subscript("p")), color=lcn[3])
    #static
    lines!(MVN_kp[!, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[!, 54], label=rich(rich("N", font = :italic), subscript("p")), color=lcn[2])
    lines!(minpajulianmod1.-MVN_kp[1, :JulUT], dropdims(nH, dims=(2,)), label=rich(rich("N", font = :italic), subscript("p")), color=lcn[1])

    ax5=Axis(fig[5,1], xlabel = "Universal Time", xticks=xtk, yticks=range(-500, 500, step=50), yminorticks=IntervalsBetween(5), yminorticksvisible = true, backgroundcolor=:gray80, ylabelpadding=10.0)
    ax5.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(minimum(xtk), maximum(xtk)) 
    ylims!(40, 220)
    #ax3.yreversed = true
    lines!(MVN_kp[!, :JulUT].-MVN_kp[1, :JulUT], MVN_kp[!, 49], label=rich(rich("T", font = :italic), subscript("p")), color=lcn[3])
    #lines!(ax3, minpajulianmod1.-MVN_kp[1, :JulUT], dropdims(TH, dims=(2,)), label=rich(rich("T", font = :italic), subscript("p")), color=lcn[1])
    
    

    xd = date .+ Dates.Hour(9) .+ Dates.Minute(0) .+ Dates.Minute(10) .* range(0, 9)
    
    xtk = datetime2julian.(xd)
    xtk = Float32.(xtk .- MVN_kp[1, :JulUT])

    ax6 = Axis(fig[6, 1],  xticks=xtk, yticks=range(-50, 50, step=5), yminorticks=IntervalsBetween(5), yminorticksvisible = true, backgroundcolor=:gray80)
    #title = rich(rich("MAVEN ", color=lcn[1]), rich("Tianwen-1", color=lcn[2]), " MSO Components"),
    ax6.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(minimum(xtk), maximum(xtk))
    ylims!(-13, 9)
    
    

    lines!(MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_mag2c1hz[!, :X_MSO], label = rich(rich("B", font = :italic), subscript("x")), color=lcn[3])
    lines!(TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_mag2c1hz[!, :X_MSO], label = rich(rich("B", font = :italic), subscript("x")), color=lcn[1])

    ax7 = Axis(fig[7, 1],  xticks=xtk, yticks=range(-50, 50, step=5), yminorticks=IntervalsBetween(5),  yminorticksvisible = true, backgroundcolor=:gray80, ylabelpadding=10.0)
    ax7.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(minimum(xtk), maximum(xtk))
    ylims!(-8.5, 8.5)

    lines!(MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_mag2c1hz[!, :Y_MSO], label = rich(rich("B", font = :italic), subscript("y")), color=lcn[3])
    lines!(TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_mag2c1hz[!, :Y_MSO], label = rich(rich("B", font = :italic), subscript("y")), color=lcn[1])

    #=
    vals = [xtk[5]]
    cen = [-1.0]
    lows = cen .- 5.0
    highs = cen .+ 5.0
    rangebars!(ax7, vals, lows, highs, color = lcn[2], whiskerwidth = 10, direction = :y)
    text!(ax7, 0.45, 0.45, text = "±5", font = :bold, align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 20)
    =#

    ax8 = Axis(fig[8, 1],  xticks=xtk, yticks=range(-50, 50, step=5), yminorticks=IntervalsBetween(5),  yminorticksvisible = true, backgroundcolor=:gray80)
    ax8.xticks = (xtk, Dates.format.(xd, "HH:MM:SS"))
    xlims!(minimum(xtk), maximum(xtk))
    ylims!(-13, 6)

    lines!(MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_mag2c1hz[!, :Z_MSO], label = rich(rich("B", font = :italic), subscript("z")), color=lcn[3])
    lines!(TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_mag2c1hz[!, :Z_MSO], label = rich(rich("B", font = :italic), subscript("z")), color=lcn[1])

 
    
    ax9=Axis(fig[9,1], xlabel = "Universal Time", xticks=xtk, yscale=log2)
    ax9.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(minimum(xtk), maximum(xtk)) 
    ylims!(3, 50)
    ax9.yreversed = true
    hp1 = heatmap!(MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_period1hz, log10.(MVN_Bpower1hz'), colormap=:gist_earth, colorrange=(-1.5, 2.5))
    lines!(MVN_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], MVN_Tcp, color = :black, linestyle=:dash)
   

    ax10=Axis(fig[10,1], xlabel = "Universal Time",  xticks=xtk, yscale=log2, ylabelpadding=10.0)
    ax10.xticks = (xtk, Dates.format.(xd, "HH:MM"))
    xlims!(minimum(xtk), maximum(xtk)) 
    ylims!(3, 50)
    ax10.yreversed = true
    hp2 = heatmap!(TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_period1hz, log10.(TW_Bpower1hz'), colormap=:gist_earth, colorrange=(-1.5, 2.5))
    lines!(TW_mag2c1hz[!, :JulUT].-MVN_kp[1, :JulUT], TW_Tcp, color = :black, linestyle=:dash)
    

    

    ax11 = Axis(fig[1:8, 2], xgridvisible = false, ygridvisible = false)
    ylims!(1, 10)
    xlims!(0, 2)
    #bracket!(0.0, 9.5, 0.0, 1.5, text = rich(rich("Tianwen-1", color=lcn[1]),  rich("     MAVEN", color=lcn[3])), orientation = :up, style = :square, width=5, rotation=pi/2)
    bracket!(0.0, 9.5, 0.0, 1.5, style = :square, width=5)
    Label(fig[4, 2], rich("MAVEN", color=lcn[3]), rotation = pi/2, tellheight=false)
    Label(fig[5, 2], rich("Tianwen-1", color=lcn[1]), rotation = pi/2, tellheight=false)
    hidespines!(ax11)
    hidedecorations!(ax11)

    Label(fig[9, 2], rich("MAVEN", color=lcn[3]), rotation = pi/2, tellheight=false)
    Label(fig[10, 2], rich("Tianwen-1", color=lcn[1]), rotation = pi/2, tellheight=false)
    Label(fig[1:3, 0], rich("Flow Velocity", "\n", "(km⋅s", superscript("-1"), ")"), rotation = pi/2, tellheight=false)
    Label(fig[4, 0], rich("Density", "\n", "(cm", superscript("-3"), ")"), rotation = pi/2, tellheight=false)
    Label(fig[5, 0], rich("Temperature", "\n", "(eV)"), rotation = pi/2, tellheight=false)
    Label(fig[6:8, 0], rich("Magnetic Field", "\n", "(nT)"), rotation = pi/2, tellheight=false)
    Label(fig[9:10, 0], rich("Wave Period", "\n", "(s)"), rotation = pi/2, tellheight=false)

    #cb = Colorbar(fig[9:10, 3], hp1, label = rich("Wavelet Power"), vertical=true, flipaxis=true, tellheight=false, height=130, labelsize=15)
    #translate!(cb.blockscene, -450, 0, 0)
    #ax9.scene.viewport
    Colorbar(fig, hp1, label = rich("Log10 Power"), 
    labelsize = 15, ticklabelsize = 15, bbox=Observable(HyperRectangle{2, Int64}([98, 46], [520, 200])), vertical=false, flipaxis=false, width=140, ticks=range(-2, 3, step=1), rightspinecolor=:white, leftspinecolor=:white, topspinecolor=:white, bottomspinecolor=:white, labelpadding=-60,
    alignmode = Outside(3), halign = :right, ticklabelcolor = :white, labelcolor = :white,tickcolor = :white)
    
    ind = [1, 4, 5, 6, 9]
    txt = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"]
    for (i, indi) in enumerate(ind)
        Label(fig[indi, -1], txt[i],font = :bold, valign=:top, tellheight=false, fontsize = 20)
    end
    ax = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10]
    linkxaxes!(ax...)
    hidexdecorations!.(ax[1:end-1], grid=false)
    rowgap!(fig.layout, 3)
    colgap!(fig.layout, 0)
    colsize!(fig.layout, 2, Relative(1/20))
    #colsize!(fig.layout, 3, Relative(1/20))
    fig
end

function orbit_plot()
    size_inches = (9, 12)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5,0, 0, 0))
    cn = cgrad(:seaborn_deep, 10, categorical = true)
    lcn = Makie.wong_colors()

    ax1 = Axis(fig[1,1:5], aspect=DataAspect(), xlabel=rich(rich("x", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich(rich("z", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10), yticks=range(-10, 10), yautolimitmargin=(0.00f0, 0.05f0), xgridvisible = false, ygridvisible = false)
    xlims!(-5, 3) 
    ylims!(-4, 4) 

    
    x = range(-10, 10, length=801)
    y = range(-10, 10, length=801)
    
    z = [exp(-sqrt(tx^2+ty^2)^1.7/2.0) for tx in x, ty in y]
    ind = findall(z .> exp(-1.0/2.0))
    z[ind] .= exp(-1.0/2.0)
	z =  round.(Int, z./maximum(z)*256)
	ind = findall(x->x<1, z)
	z[ind] .= 1

	
	imcn = reverse(cgrad(:bilbao, 256, categorical = true))
	imz = imcn[z]
    image!((-10, 10), (-10, 10), imz)

    img = load("marsxz1.png")
    image!((-1, 1), (-1, 1), rotr90(img))

    θ = range(0, 2π, length=3600)
    lines!(cos.(θ), sin.(θ), color=imcn[z[441,401]], linewidth=1.0, overdraw = true)



    ax2 = Axis(fig[2,1], aspect=DataAspect(), xlabel=rich(rich("x", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich(rich("z", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10, step=2), yticks=range(-10, 10, step=2), yautolimitmargin=(0.00f0, 0.05f0), xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticksvisible = true, backgroundcolor=:gray80)
    xlims!(-4.5, 1.5) 
    ylims!(-4, 4) 

    img = load("marsxz1.png")
    image!((-1, 1), (-1, 1), rotr90(img))
    
    lines!(TW_mag2c1hz[!, :Probe_Position_X_MSO], TW_mag2c1hz[!, :Probe_Position_Z_MSO]; color=lcn[1], overdraw = true)

    ind1 = findall(x->x>0 , MVN_mag2c1hz[!, :Probe_Position_Y_MSO])
    r = @. sqrt(MVN_mag2c1hz[!, :Probe_Position_X_MSO]^2 + MVN_mag2c1hz[!, :Probe_Position_Z_MSO]^2)
    ind2 = findall(x->x<1 , r)
    ind = intersect(ind1, ind2)
    tx = MVN_mag2c1hz[:, :Probe_Position_X_MSO]
    tz = MVN_mag2c1hz[:, :Probe_Position_Z_MSO]
    tx[ind] .= NaN64
    tz[ind] .= NaN64
    lines!(tx, tz; color=lcn[3], overdraw = true)

    for ti in 0:4
        color = cn[ ti+4 == 8 ? 2 : ti+4]
        local ind = findall(x-> x<date + Dates.Hour(9) + Dates.Minute(30+ti*15) + Dates.Second(1)  && x>date + Dates.Hour(9) + Dates.Minute(30+ti*15) - Dates.Second(1), MVN_mag2c1hz[!, :Time]) 
        local str = Dates.format(date + Dates.Hour(9) + Dates.Minute(30+ti*15), "HH:MM")
        ind = ind[1]
        scatter!(MVN_mag2c1hz[ind, :Probe_Position_X_MSO], MVN_mag2c1hz[ind, :Probe_Position_Z_MSO], markersize=10, label=str, color=color, transparency=true)   
       
        local ind = findall(x-> x<date + Dates.Hour(9) + Dates.Minute(30+ti*15) + Dates.Second(1)  && x>date + Dates.Hour(9) + Dates.Minute(30+ti*15) - Dates.Second(1), TW_mag2c1hz[!, :Time]) 
        ind = ind[1]
        local str = Dates.format(date + Dates.Hour(9) + Dates.Minute(30+ti*15), "HH:MM")
        scatter!(TW_mag2c1hz[ind, :Probe_Position_X_MSO], TW_mag2c1hz[ind, :Probe_Position_Z_MSO], markersize=10, label=str, color=color, transparency=true)
        text!(TW_mag2c1hz[ind, :Probe_Position_X_MSO]+0.35, TW_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.20, text=str, align = (:left, :top), space = :data, fontsize = 10, color=color, rotation=deg2rad(-60-ti*7))
        
        
    end
    



    ax3 = Axis(fig[2,2], aspect=DataAspect(), xlabel=rich(rich("y", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich(rich("z", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10, step=2), yticks=range(-10, 10, step=2), yautolimitmargin=(0.00f0, 0.05f0), xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticksvisible = true, backgroundcolor=:gray80)
    xlims!(-2, 2) 
    ylims!(-4, 4) 

    img = load("marsyz1.png")
    image!((-1, 1), (-1, 1), rotr90(img))
    
    ind1 = findall(x->x<0 , TW_mag2c1hz[!, :Probe_Position_X_MSO])
    r = @. sqrt(TW_mag2c1hz[!, :Probe_Position_Y_MSO]^2 + TW_mag2c1hz[!, :Probe_Position_Z_MSO]^2)
    ind2 = findall(x->x<1 , r)
    ind = intersect(ind1, ind2)
    tx = TW_mag2c1hz[:, :Probe_Position_Y_MSO]
    tz = TW_mag2c1hz[:, :Probe_Position_Z_MSO]
    tx[ind] .= NaN64
    tz[ind] .= NaN64
    lines!(tx, tz; color=lcn[1], overdraw = true)

 
    tx = MVN_mag2c1hz[:, :Probe_Position_Y_MSO]
    tz = MVN_mag2c1hz[:, :Probe_Position_Z_MSO]
    lines!(tx, tz; color=lcn[3], overdraw = true)

    for ti in 0:4
        color = cn[ ti+4 == 8 ? 2 : ti+4]
        local ind = findall(x-> x<date + Dates.Hour(9) + Dates.Minute(30+ti*15) + Dates.Second(1)  && x>date + Dates.Hour(9) + Dates.Minute(30+ti*15) - Dates.Second(1), MVN_mag2c1hz[!, :Time]) 
        local str = Dates.format(date + Dates.Hour(9) + Dates.Minute(30+ti*15), "HH:MM")
        ind = ind[1]
        
        scatter!(MVN_mag2c1hz[ind, :Probe_Position_Y_MSO], MVN_mag2c1hz[ind, :Probe_Position_Z_MSO], markersize=10, label=str, color=color, transparency=true) 
          
        
        local ind = findall(x-> x<date + Dates.Hour(9) + Dates.Minute(30+ti*15) + Dates.Second(1)  && x>date + Dates.Hour(9) + Dates.Minute(30+ti*15) - Dates.Second(1), TW_mag2c1hz[!, :Time]) 
        ind = ind[1]
        local str = Dates.format(date + Dates.Hour(9) + Dates.Minute(30+ti*15), "HH:MM")
        scatter!(TW_mag2c1hz[ind, :Probe_Position_Y_MSO], TW_mag2c1hz[ind, :Probe_Position_Z_MSO], markersize=10, label=str, color=color, transparency=true)

        text!(MVN_mag2c1hz[ind, :Probe_Position_Y_MSO]-0.00, MVN_mag2c1hz[ind, :Probe_Position_Z_MSO]-0.30, text=str, align = (:left, :center), space = :data, fontsize = 10, color=color, rotation=deg2rad(-85))  
    end
    

    

    ax4 = Axis(fig[2, 4], aspect=DataAspect(), xlabel=rich(rich("x", font = :italic), " (", rich("R", font = :italic), subscript("M"), ")"), ylabel=rich("(", rich("y", font = :italic), superscript("2"), "+", rich("z", font = :italic), superscript("2"), ")", superscript("1/2"), " (", rich("R", font = :italic), subscript("M"), ")"), xticks=range(-10, 10), yticks=range(-10, 10), yautolimitmargin=(0.00f0, 0.05f0), xminorticks=IntervalsBetween(5), yminorticks=IntervalsBetween(5), xminorticksvisible = true, yminorticksvisible = true)
    xlims!(-5, 3) 
    ylims!(0, 4) 

    

    hmVA = heatmap!(X, Y, MVN_VA, colormap=:gist_earth, nan_color=:grey80, colorrange=(140, 180))

    img = load("marsxz1.png")
    image!((-1, 1), (-1, 1), rotr90(img))
    lines!(xshock, ryzshock; color=cn[8], linewidth=4, overdraw = true)
    lines!(xmp, ryzmp; color=cn[8], linewidth=4,  overdraw = true)
    lines!(TW_mag2c1hz[!, :Probe_Position_X_MSO], sqrt.(TW_mag2c1hz[!, :Probe_Position_Y_MSO].^2+TW_mag2c1hz[!, :Probe_Position_Z_MSO].^2); color=lcn[1], overdraw = true)
    lines!(MVN_mag2c1hz[!, :Probe_Position_X_MSO], sqrt.(MVN_mag2c1hz[!, :Probe_Position_Y_MSO].^2+MVN_mag2c1hz[!, :Probe_Position_Z_MSO].^2); color=lcn[3], overdraw = true)

    for ri in 2:4
        @show ri
        lx, lz, va = stream(3.0, 0.75ri)
        lines!(lx, lz; color=cn[9], overdraw = true, linewidth=2)
        lines!(lx, -lz; color=cn[9], overdraw = true, linewidth=2)
        arrows!([3.0], [lz[end]], [-1.0], [0.0], arrowsize = 15, lengthscale = 0.5, color=cn[9])
    end

    for ti in 0:4
        local ind = findall(x-> x<date + Dates.Hour(9) + Dates.Minute(30+ti*15) + Dates.Second(1)  && x>date + Dates.Hour(9) + Dates.Minute(30+ti*15) - Dates.Second(1), MVN_mag2c1hz[!, :Time]) 
        local str = Dates.format(date + Dates.Hour(9) + Dates.Minute(30+ti*15), "HH:MM")
        ind = ind[1]
        color = cn[ ti+4 == 8 ? 2 : ti+4]
        scatter!(MVN_mag2c1hz[ind, :Probe_Position_X_MSO], sqrt.(MVN_mag2c1hz[ind, :Probe_Position_Y_MSO].^2+MVN_mag2c1hz[ind, :Probe_Position_Z_MSO].^2), markersize=10, label=str, color=color, transparency=true)   
        if ti == 0 
            text!(0.13, 2.3, text=str, align = (:right, :center), space = :data, fontsize = 10, color=color, rotation=deg2rad(-40))
        elseif ti == 1
            text!(0.35, 2.38, text=str, align = (:left, :center), space = :data, fontsize = 10, color=color, rotation=deg2rad(70))
        elseif ti == 2 
            text!(0.46, 2.22, text=str, align = (:left, :center), space = :data, fontsize = 10, color=color, rotation=deg2rad(40))
        elseif ti == 3 
            text!(0.46, 2.12, text=str, align = (:left, :center), space = :data, fontsize = 10, color=color, rotation=deg2rad(-10))
        elseif ti == 4 
            text!(0.44, 1.90, text=str, align = (:left, :center), space = :data, fontsize = 10, color=color, rotation=deg2rad(-10))
        end

        local ind = findall(x-> x<date + Dates.Hour(9) + Dates.Minute(30+ti*15) + Dates.Second(1)  && x>date + Dates.Hour(9) + Dates.Minute(30+ti*15) - Dates.Second(1), TW_mag2c1hz[!, :Time]) 
        ind = ind[1]
        local str = Dates.format(date + Dates.Hour(9) + Dates.Minute(30+ti*15), "HH:MM")
        scatter!(TW_mag2c1hz[ind, :Probe_Position_X_MSO], sqrt.(TW_mag2c1hz[ind, :Probe_Position_Y_MSO].^2+TW_mag2c1hz[ind, :Probe_Position_Z_MSO].^2), markersize=10, label=str, color=color, transparency=true)
        text!(TW_mag2c1hz[ind, :Probe_Position_X_MSO]-0.18, sqrt.(TW_mag2c1hz[ind, :Probe_Position_Y_MSO].^2+TW_mag2c1hz[ind, :Probe_Position_Z_MSO].^2)+0.22, text=str, align = (:left, :top), space = :data, fontsize = 10, color=color, rotation=deg2rad(70))
        
    end

    cb = Colorbar(fig[2, 5], vertical=false, flipaxis = true,  hmVA, label="Flow                               " * "\n" * "Angle                              " * "\n" * "(Deg)                              ", tellwidth=false, width=100, tellheight=false, labelsize=10, ticklabelsize=10)
    translate!(cb.blockscene, -243, -60, 0)

    str=[Dates.format(date + Dates.Hour(9) + Dates.Minute(30+ti*15), "HH:MM ") for ti in 0:4]
    color = [cn[ ti+4 == 8 ? 2 : ti+4] for ti in 0:4]
    st = rich(rich(str[1], color=color[1]), rich(str[2], color=color[2]), rich(str[3], color=color[3]), rich(str[4], color=color[4]), rich(str[5], color=color[5]))

    text!(ax2, 0.0, 1.0, text = rich(rich("Tianwen-1", color=lcn[1])), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 10)
    text!(ax2, 0.2, -2.5, text = rich(rich("MAVEN", color=lcn[3])), align = (:right, :top), fontsize = 10)

    text!(ax4, -3.5, 1.7, text = rich(rich("Tianwen-1", color=lcn[1])), align = (:left, :bottom), fontsize = 10)
    text!(ax4, -0.3, 1.0, text = rich(rich("MAVEN", color=lcn[3])), align = (:right, :bottom), fontsize = 10)

    
    #Label(fig[2, 1:3], st, tellheight=false, fontsize=15)

    linkyaxes!(ax2, ax3)
    hideydecorations!(ax3, grid=false)
    rowsize!(fig.layout, 1, Relative(3.17/4))
    #rowsize!(fig.layout, 2, Relative(0.01/4))
    colsize!(fig.layout, 1, Relative(4.37/2/10))
    colsize!(fig.layout, 2, Relative(4.37/3/10))
    colsize!(fig.layout, 5, Relative(0.001/10))
    rowgap!(fig.layout, 5)
    colgap!(fig.layout, 1, 3)
    colgap!(fig.layout, 2, 5)
    
    txt = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m"]
    ax = [ax1, ax2, ax3, ax4]
    

    Label(fig[1, 0], txt[1],font = :bold, valign=:top, tellheight=false, fontsize = 20)

    Label(fig[2, 0], txt[2],font = :bold, valign=:top, tellheight=false, fontsize = 20)

    Label(fig[2, 3], txt[3], font = :bold, valign=:top, tellheight=false, fontsize = 20)

    fig 
end



#=
Ephem = Ephem_plot()
save("Ephem"* datestr * ".pdf", Ephem, pt_per_unit = 1)

SWSta = SWSta_plot()
save("SWSta"* ".pdf", SWSta, pt_per_unit = 1)

Mach = Mach_plot()
save("Mach"* ".pdf", Mach, pt_per_unit = 1)

SolarCycle = SolarCycle_plot()
save("SolarCycle"* ".pdf", SolarCycle, pt_per_unit = 1)

streamline = streamline_plot()
save("stream" * ".pdf", streamline, pt_per_unit = 1)

SolarWind_Overview = SolarWind_Overview_plot()
save("SolarWind_Overview.pdf", SolarWind_Overview, pt_per_unit = 1)

orbit = orbit_plot()
save("TW_MVN_orbit"* datestr * ".pdf", orbit, pt_per_unit = 1)

orbit = orbitxyz_plot()
save("TW_MVN_orbitxyz"* datestr * ".pdf", orbit, pt_per_unit = 1)

orbit = orbitxyz_pre_plot()
save("TW_MVN_orbitxyz_pre"* datestr * ".pdf", orbit, pt_per_unit = 1)

schematic = test2_schematic_plot()
save("test2_schematic"* datestr * ".pdf", schematic, pt_per_unit = 1)

schematic = schematic_plot()
save("schematic"* datestr * ".pdf", schematic, pt_per_unit = 1)

MVN_TW_Plasma_Overview = MVN_TW_Plasma_Overview_plot()
save("MVN_TW_Plasma_Overview"* datestr * ".pdf", MVN_TW_Plasma_Overview, pt_per_unit = 1)

MVN_TW_Mag_Overview = MVN_TW_Mag_Overview_plot()
save("MVN_TW_Mag_Overview"* datestr * ".pdf", MVN_TW_Mag_Overview, pt_per_unit = 1)

MVN_TW_Overview = MVN_TW_Overview_plot()
save("MVN_TW_Overview"* datestr * ".pdf", MVN_TW_Overview, pt_per_unit = 1)

MVN_TW_Mag_Pre_Overview = MVN_TW_Mag_Pre_Overview_plot()
save("MVN_TW_Mag_Pre_Overview"* datestr * ".pdf", MVN_TW_Mag_Pre_Overview, pt_per_unit = 1)

im=1
MVN_STATIC_eflux_Overview = MVN_STATIC_eflux_Overview_plot(im)
save("MVN_STATIC_H+eflux_Overview"* datestr * ".pdf", MVN_STATIC_eflux_Overview, pt_per_unit = 1)
im=2
MVN_STATIC_eflux_Overview = MVN_STATIC_eflux_Overview_plot(im)
save("MVN_STATIC_He++eflux_Overview"* datestr * ".pdf", MVN_STATIC_eflux_Overview, pt_per_unit = 1)
im=5
MVN_STATIC_eflux_Overview = MVN_STATIC_eflux_Overview_plot(im)
save("MVN_STATIC_O+eflux_Overview"* datestr * ".pdf", MVN_STATIC_eflux_Overview, pt_per_unit = 1)
im=6
MVN_STATIC_eflux_Overview = MVN_STATIC_eflux_Overview_plot(im)
save("MVN_STATIC_O2+eflux_Overview"* datestr * ".pdf", MVN_STATIC_eflux_Overview, pt_per_unit = 1)
im=7
MVN_STATIC_eflux_Overview = MVN_STATIC_eflux_Overview_plot(im)
save("MVN_STATIC_CO2+eflux_Overview"* datestr * ".pdf", MVN_STATIC_eflux_Overview, pt_per_unit = 1)


MVN_STATIC_efluxBN = MVN_STATIC_efluxBN_plot()
save("MVN_STATIC_efluxBN"* datestr * ".pdf", MVN_STATIC_efluxBN, pt_per_unit = 1)
MVN_STATIC_efluxEB = MVN_STATIC_efluxEB_plot()
save("MVN_STATIC_efluxEB"* datestr * ".pdf", MVN_STATIC_efluxEB, pt_per_unit = 1)
MVN_STATIC_efluxEN = MVN_STATIC_efluxEN_plot()
save("MVN_STATIC_efluxEN"* datestr * ".pdf", MVN_STATIC_efluxEN, pt_per_unit = 1)


MVN_SWIA_eflux_Overview = MVN_SWIA_eflux_Overview_plot()
save("MVN_SWIA_eflux_Overview"* datestr * ".pdf", MVN_SWIA_eflux_Overview, pt_per_unit = 1)

MVN_SWIA_efluxBN = MVN_SWIA_efluxBN_plot()
save("MVN_SWIA_efluxBN"* datestr * ".pdf", MVN_SWIA_efluxBN, pt_per_unit = 1)
MVN_SWIA_efluxEB = MVN_SWIA_efluxEB_plot()
save("MVN_SWIA_efluxEB"* datestr * ".pdf", MVN_SWIA_efluxEB, pt_per_unit = 1)
MVN_SWIA_efluxEN = MVN_SWIA_efluxEN_plot()
save("MVN_SWIA_efluxEN"* datestr * ".pdf", MVN_SWIA_efluxEN, pt_per_unit = 1)

MVN_STATIC_SWIA_efluxmean = MVN_STATIC_SWIA_efluxmean_plot()
save("MVN_STATIC_SWIA_efluxmean"* datestr * ".pdf", MVN_STATIC_SWIA_efluxmean, pt_per_unit = 1)

MVN_STATIC_SWIA_efluxmean_profile = MVN_STATIC_SWIA_efluxmean_profile_plot()
save("MVN_STATIC_SWIA_efluxmean_profile"* datestr * ".pdf", MVN_STATIC_SWIA_efluxmean_profile, pt_per_unit = 1)

MVN_STATIC_efluxmean = MVN_STATIC_efluxmean_plot()
save("MVN_STATIC_efluxmean"* datestr * ".pdf", MVN_STATIC_efluxmean, pt_per_unit = 1)
MVN_SWIA_efluxmean = MVN_SWIA_efluxmean_plot()
save("MVN_SWIA_efluxmean"* datestr * ".pdf", MVN_SWIA_efluxmean, pt_per_unit = 1)

MVN_STATIC_efluxmean_heatmap = MVN_STATIC_efluxmean_heatmap_plot()
save("MVN_STATIC_efluxmean_heatmap"* datestr * ".pdf", MVN_STATIC_efluxmean_heatmap, pt_per_unit = 1)

im = 1
TW_MINPA_eflux_Overview = TW_MINPA_eflux_Overview_plot(im)
save("TW_MINPA_H+eflux_Overview"* datestr * ".pdf", TW_MINPA_eflux_Overview, pt_per_unit = 1)
im = 2
TW_MINPA_eflux_Overview = TW_MINPA_eflux_Overview_plot(im)
save("TW_MINPA_He++eflux_Overview"* datestr * ".pdf", TW_MINPA_eflux_Overview, pt_per_unit = 1)
im = 4
TW_MINPA_eflux_Overview = TW_MINPA_eflux_Overview_plot(im)
save("TW_MINPA_O+eflux_Overview"* datestr * ".pdf", TW_MINPA_eflux_Overview, pt_per_unit = 1)
im = 6
TW_MINPA_eflux_Overview = TW_MINPA_eflux_Overview_plot(im)
save("TW_MINPA_O2+eflux_Overview"* datestr * ".pdf", TW_MINPA_eflux_Overview, pt_per_unit = 1)
im = 7
TW_MINPA_eflux_Overview = TW_MINPA_eflux_Overview_plot(im)
save("TW_MINPA_CO2+eflux_Overview"* datestr * ".pdf", TW_MINPA_eflux_Overview, pt_per_unit = 1)

TW_MINPA_efluxmean = TW_MINPA_efluxmean_plot()
save("TW_MINPA_efluxmean"* datestr * ".pdf", TW_MINPA_efluxmean, pt_per_unit = 1)

fig, swcore, swpu, swpusc = TW_MINPA_eflux_fitting_plot()
save("TW_MINPA_eflux_fitting"* datestr * ".pdf", fig, pt_per_unit = 1)
fig = TW_MINPA_eflux_model_data_comparision_plot(swcore, swpu, swpusc)
save("TW_MINPA_eflux_model_data_comparision"* datestr * ".pdf", fig, pt_per_unit = 1)


fig, swcore, swpu, swpusc = MVN_STATIC_eflux_fitting_plot()
save("MVN_STATIC_eflux_fitting"* datestr * ".pdf", fig, pt_per_unit = 1)
fig = MVN_STATIC_eflux_model_data_comparision_plot(swcore, swpu, swpusc, mspu)
save("MVN_STATIC_eflux_model_data_comparision"* datestr * ".pdf", fig, pt_per_unit = 1)

NGITemperature = NGITemperature_plot()
save("NGITemperature.pdf", NGITemperature, pt_per_unit = 1)



=#