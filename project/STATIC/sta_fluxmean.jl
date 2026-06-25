using Glob
using CommonDataFormat
using Dates
using Statistics
using CairoMakie
include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_plot.jl")
using .MAVEN_load
using .MAVEN_plot

date = Date(2022, 8, 24)
datestr = Dates.format.(date, "yyyymmdd")
year_num = Dates.format(date, "yyyy")
month_num = Dates.format(date, "mm")

data_dir = "E:\\Data"
output_dir = "E:\\Output"
MAVEN_path = joinpath(data_dir, "MAVEN")
STATIC_path = joinpath(MAVEN_path, "STATIC", "l2", year_num, month_num)
MAVEN_output_path = joinpath(output_dir, "MAVEN")
STATIC_output_path = joinpath(MAVEN_output_path, "STATIC")
fluxmean_path = joinpath(STATIC_output_path, "fluxmean", year_num, month_num)
mkpath(fluxmean_path)

STA_file = glob("mvn_sta_l2*32e4d16a8m*$datestr*.cdf", STATIC_path)
println("STATIC file: ", STA_file)
println("Output path: ", fluxmean_path)

sta_data = MAVEN_load.load_STATIC(STA_file[1])

sp_title = ["H+", "He++", "m/q=4.573966", "m/q=9.22163", 
    "O+", "O2+", "CO2+", "m/q=69.26344"]
sp = [rich(rich("H"), superscript("+")), rich(rich("He"), superscript("++")), 
    rich("m/q=4.573966"), rich("m/q=9.22163"), rich(rich("O"), superscript("+")), 
    rich(rich("O"), subscript("2"), superscript("+")), rich(rich("CO"), 
    subscript("2"), superscript("+")), rich("m/q=69.26344")]

img_list = [1, 2, 5, 6, 7]
time_range = DateTime(date, Time(7, 30, 00)).+Dates.Minute(10) .* range(0, 1)
idx = findall(time_range[1] .<= sta_data[:epoch] .<= time_range[end])
swp = sta_data[:swp_ind][idx[1]]

using CairoMakie
CairoMakie.activate!()
set_theme!(;
    Axis=(xscale = log10, yscale = log10,
        backgroundcolor = :oldlace,
        bordercolor=:black, bordersize=15, spinewidth = 2, 
        xlabelsize=18, ylabelsize=18, xticklabelsize=16, yticklabelsize=16, 
        xminorticksvisible = true, yminorticksvisible = true,
        xminorticks = IntervalsBetween(10), yminorticks = IntervalsBetween(10))
)

import Logging
Logging.with_logger(Logging.NullLogger()) do
for img in img_list
    energy = sta_data[:energy][img,:,:,swp][1,:]
    phi = sta_data[:phi][img,:,:,swp]
    theta = sta_data[:theta][img,:,:,swp]
    flux = sta_data[:eflux][idx, img, :,:]
    fluxmean = dropdims(mean(flux, dims=(1)), dims=(1))
    fluxmean_phi = dropdims(mean(reshape(fluxmean, 4, 16, 32), dims=1), dims=1)
    fluxmean_theta = [mean(fluxmean[i:4:end, j]) for i in 1:4, j in 1:32]

    title = "MAVEN STATIC " * sp[img] * "Energy Flux Mean " * 
        Dates.format(time_range[1], "yyyymmddTHH:MM")
    name = "MVN_STATIC_" * sp_title[img] * "eflux_mean_" * 
        Dates.format(time_range[1], "yyyymmddTHHMM") * ".pdf"

    fig = Figure(resolution = (800, 1200))
    ax1 = Axis(fig[1, 1];xlabel="Energy (eV)" ,ylabel="Energy Flux")
    ax2 = Axis(fig[2, 1];xlabel="Energy (eV)" ,ylabel="Energy Flux")
    colors = cgrad(:gist_earth, 16, categorical = true)
    ylims!(ax1, 1e3, 1e9)
    for i in 1:16
    scatterlines!(ax1, energy, fluxmean_phi[i, :], color = colors[i], label = "$i", 
        linewidth = 3, marker = :circle, markersize = 10)
    end
    # 添加图例
    Legend(fig[1, 2], ax1, "Channel", nbanks = 2)
    # --- 绘制下图 (示例：只选前 4 个通道) ---
    ylims!(ax2, 1e2, 1e9)
    for i in 1:4
    scatterlines!(ax2, energy, fluxmean_theta[i, :], color = colors[i], label = "$i", 
        linewidth = 3, marker = :circle, markersize = 10)
    end
    Legend(fig[2, 2], ax2, "Channel")
    Label(fig[0, :], title, fontsize=20, tellwidth=false)

    save(joinpath(fluxmean_path, name), fig, pt_per_unit = 1)
    println("Saved figure: ", joinpath(fluxmean_path, name))
end
end