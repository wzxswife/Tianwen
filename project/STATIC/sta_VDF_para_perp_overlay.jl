using Glob
using CommonDataFormat
using Dates
using Statistics, LinearAlgebra
using Rotations
using CairoMakie
include("../../src/scripts/MAVEN_STATIC.jl")
include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_plot.jl")
using .MAVEN_STATIC
using .MAVEN_load
using .MAVEN_plot

date = Date(2022, 8, 24)
datestr = Dates.format.(date, "yyyymmdd")
year_num = Dates.format(date, "yyyy")
month_num = Dates.format(date, "mm")

ts = DateTime.(date + Dates.Hour(7) + Dates.Minute(20))
time_stemp = ts .+ Dates.Minute(1) .* range(0, 10)

data_dir = "E:\\Data"
output_dir = "E:\\Output"
MAVEN_path = joinpath(data_dir, "MAVEN")
MAG_path = joinpath(MAVEN_path, "MAG", "l2", year_num, month_num)
STATIC_path = joinpath(MAVEN_path, "STATIC", "l2", year_num, month_num)
MAVEN_output_path = joinpath(output_dir, "MAVEN")
STATIC_output_path = joinpath(MAVEN_output_path, "STATIC")
VDF_path = joinpath(STATIC_output_path, "VDF", year_num, month_num)
mkpath(VDF_path)
MAG_file = glob("mvn_mag_l2*$datestr*.sts", MAG_path)
STA_file = glob("mvn_sta_l2*d1*32e4d16a8m*$datestr*.cdf", STATIC_path)
println("MAG file: ", MAG_file)
println("STATIC file: ", STA_file)
println("Output path: ", VDF_path)

raw_sta = load_cdf(STA_file[1])
mag_data = MAVEN_load.load_mag_l2(MAG_file[1])

plot_range = 1000  # km/s
colorrange = (1e-10, 1e9)

fig_name = "MAVEN_STATIC_H+_VDF_overlay_" * Dates.format(time_stemp[1], "yyyymmddTHHMM")*
    "-"*Dates.format(time_stemp[end], "HHMM")*".pdf"
title = "MAVEN STATIC H+ VDF overlay " * Dates.format(time_stemp[1], "yyyymmddTHH:MM")*
    "-"*Dates.format(time_stemp[end], "HH:MM")*".pdf"


fig = Figure(size = (1000, 600))
ax = Axis(fig[1, 1], aspect=DataAspect(),title=title, 
    xlabel=L"$V_{\parallel}$ (km/s)", ylabel=L"$V_{\perp}$ (km/s)",
    limits=(-plot_range, plot_range, 0, plot_range),
    xticks=[-1000, -500, 0, 500, 1000], yticks=[0, 500, 1000])

v_par_all = Float64[]
v_perp_all = Float64[]
df_all = Float64[]
vel_par_all = Float64[]
vel_perp_all = Float64[]

for (col_idx, II) in enumerate(eachindex(time_stemp))
    println("==== $(time_stemp[II]) ====")

    _, time_ind = findmin(d -> abs(d - time_stemp[II]), raw_sta[:epoch])
    swi_time = raw_sta[:epoch][time_ind]

    _, time_b_ind = findmin(d -> abs(d - swi_time), mag_data[:epoch])
    b0 = mag_data[:B][time_b_ind, :]
    vsc = (mag_data[:position][time_b_ind+1, :] .- mag_data[:position][time_b_ind-1, :]) ./
          (datetime2unix(mag_data[:epoch][time_b_ind+1]) - datetime2unix(mag_data[:epoch][time_b_ind-1]))

    slip = MAVEN_STATIC.static_slip(raw_sta, time_ind)
    slip[:valid] == 0 && continue
    MAVEN_STATIC.STA_count2df(slip)
    MAVEN_STATIC.static_rotation(slip)

    b_hat = slip[:magf] ./ norm(slip[:magf])
    mass_range = [0, 1.55]
    m_int = 1

    mass_mask = (slip[:mass_arr] .>= mass_range[1]) .& (slip[:mass_arr] .<= mass_range[2])

    energy_grid = slip[:energy][1, :, :] .+ slip[:sc_pot]
    energy_grid[energy_grid .<= 0] .= 0.001
    theta_grid = slip[:theta][1, :, :]
    phi_grid   = slip[:phi][1, :, :]

    v_mag = MAVEN_STATIC.ion_energy2v.(energy_grid, m_int) ./ 1000

    v_x = v_mag .* cosd.(theta_grid) .* cosd.(phi_grid) .+ vsc[1]
    v_y = v_mag .* cosd.(theta_grid) .* sind.(phi_grid) .+ vsc[2]
    v_z = v_mag .* sind.(theta_grid) .+ vsc[3]

    df_species = copy(slip[:df])
    df_species[.!mass_mask] .= 0.0
    df_species = sum(df_species; dims=1)[1, :, :]

    v_flat = hcat(vec(v_x), vec(v_y), vec(v_z))
    df_flat = vec(df_species)
    v_par_flat = v_flat * b_hat
    v_perp_flat = vec(sqrt.(sum((v_flat .- v_par_flat * b_hat').^2, dims=2)))

    valid = (df_flat .> 1e-30) .& (.!isnan.(df_flat)) .& (v_perp_flat .< plot_range * 1.1)
    append!(v_par_all,   v_par_flat[valid])
    append!(v_perp_all,  v_perp_flat[valid])
    append!(df_all,      df_flat[valid])

    vel_sp, _, _ = MAVEN_STATIC.sta_v_4d(slip; mass_range=mass_range, m_int=m_int, unit_cover=true)
    vel_total = vel_sp .+ vsc
    push!(vel_par_all,  dot(vel_total, b_hat))
    push!(vel_perp_all, norm(vel_total - dot(vel_total, b_hat) * b_hat))
end

# ---- 统一绘图 ----
if length(v_par_all) > 0
    cr_in = log10.(colorrange)
    v_in  = log10.(df_all)
    colors_idx = clamp.(round.(Int, ((v_in .- cr_in[1]) ./ (cr_in[2] - cr_in[1])) .* 255 .+ 1), 1, 256)

    scatter!(ax, v_par_all, v_perp_all;
        color=colors_idx, colormap=:viridis,
        markersize=10, strokewidth=0)
end

# 漂移速度
# scatter!(ax, vel_par_all, vel_perp_all;
#     color=:cyan, markersize=12, marker=:cross, strokewidth=3)

# 参考线
# vlines!(ax, 0; color=:white, linestyle=:dash, linewidth=0.8)
# hlines!(ax, 0; color=:white, linestyle=:dash, linewidth=0.8)

# 色标
Colorbar(fig[1, 2]; limits=colorrange,
    label=L"f(\mathbf{v})\ (\mathrm{s^3/km^6})",
    colormap=:viridis, scale=log10)

save(joinpath(VDF_path, fig_name), fig)
println("Saved: ", joinpath(VDF_path, fig_name))
