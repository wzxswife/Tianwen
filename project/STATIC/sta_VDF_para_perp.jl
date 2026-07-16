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

ts = DateTime.(date + Dates.Hour(7) + Dates.Minute(30))
time_stemp = ts .+ Dates.Minute(3) .* range(0, 3)

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

raw_sta = load_cdf(STA_file[1])             # 原始 CDF，保留 :mass_arr
mag_data = MAVEN_load.load_mag_l2(MAG_file[1]) # MAG L2

species = [
    Dict(:name => "H⁺",  :mass_range => [0, 1.55], :m_int => 1,  :color => :red),
    # Dict(:name => "He²⁺",:mass_range => [1.55, 2.7],:m_int => 2,  :color => :blue),
    # Dict(:name => "O⁺",  :mass_range => [14, 20],   :m_int => 16, :color => :green),
    # Dict(:name => "O₂⁺", :mass_range => [24, 40],   :m_int => 32, :color => :orange),
]

n_species = length(species)
n_times = length(time_stemp)

plot_range = 1000  # km/s
colorrange = (1e-13, 1e6)

fig_VDF = Figure(size = (250 * n_times, 250 * n_species))
fig_name = "MAVEN_STATIC_H+_VDF_" * Dates.format(time_stemp[1], "yyyymmddTHHMM")*
    "-"*Dates.format(time_stemp[end], "HHMM")*".pdf"

function map_color(vals, cr; scaler="log10")
    cr_in = scaler == "log10" ? log10.(cr) : cr
    v_in  = scaler == "log10" ? log10.(vals) : vals
    mapped = round.(Int, ((v_in .- cr_in[1]) ./ (cr_in[2] - cr_in[1])) .* 255 .+ 1)
    clamp!(mapped, 1, 256)
    return mapped
end

for (col_idx, II) in enumerate(eachindex(time_stemp))
    println("==== $(time_stemp[II]) ====")

    _, time_ind = findmin(d -> abs(d - time_stemp[II]), raw_sta[:epoch])
    swi_time = raw_sta[:epoch][time_ind]

    _, time_b_ind = findmin(d -> abs(d - swi_time), mag_data[:epoch])
    b0 = mag_data[:B][time_b_ind, :]
    vsc = (mag_data[:position][time_b_ind+1, :] .- mag_data[:position][time_b_ind-1, :]) ./
          (datetime2unix(mag_data[:epoch][time_b_ind+1]) - datetime2unix(mag_data[:epoch][time_b_ind-1]))

    # STATIC 时间切片
    slip = MAVEN_STATIC.static_slip(raw_sta, time_ind)
    slip[:valid] == 0 && continue
    MAVEN_STATIC.STA_count2df(slip)       # 添加 :df（相空间密度）
    MAVEN_STATIC.static_rotation(slip)                 # theta/phi/magf → MSO

    b_hat = slip[:magf] ./ norm(slip[:magf])

    for (sp_idx, sp) in enumerate(species)
        mass_range = sp[:mass_range]
        m_int = sp[:m_int]

        mass_mask = (slip[:mass_arr] .>= mass_range[1]) .& (slip[:mass_arr] .<= mass_range[2])
        nmass, nbins, nenergy = size(slip[:mass_arr])

        # 能量/角度各质量通道相同，取第一质量通道（去除质量维）
        energy_grid = slip[:energy][1, :, :] .+ slip[:sc_pot]
        energy_grid[energy_grid .<= 0] .= 0.001
        theta_grid = slip[:theta][1, :, :]
        phi_grid   = slip[:phi][1, :, :]

        # 速度大小 (km/s)
        v_mag = MAVEN_STATIC.ion_energy2v.(energy_grid, m_int) ./ 1000

        # 笛卡尔速度 → Vsc 修正（nbins × nenergy）
        v_x = v_mag .* cosd.(theta_grid) .* cosd.(phi_grid) .+ vsc[1]
        v_y = v_mag .* cosd.(theta_grid) .* sind.(phi_grid) .+ vsc[2]
        v_z = v_mag .* sind.(theta_grid) .+ vsc[3]

        # 物种相空间密度（按质量求和）
        df_species = copy(slip[:df])
        df_species[.!mass_mask] .= 0.0
        df_species = sum(df_species; dims=1)[1, :, :]

        # 展平
        v_flat = hcat(vec(v_x), vec(v_y), vec(v_z))
        df_flat = vec(df_species)
        v_par_flat = v_flat * b_hat
        v_perp_flat = vec(sqrt.(sum((v_flat .- v_par_flat * b_hat').^2, dims=2)))

        # 剔除无效/超范围点
        valid = (df_flat .> 1e-30) .& (.!isnan.(df_flat)) .& (v_perp_flat .< plot_range * 1.1)
        v_par   = v_par_flat[valid]
        v_perp  = v_perp_flat[valid]
        df_vals = df_flat[valid]

        # 该物种整体漂移速度 (sta_v_4d 返回速度在仪器坐标系, 加 vsc 转 MSO)
        vel_sp, _, _ = MAVEN_STATIC.sta_v_4d(slip; mass_range=mass_range, m_int=m_int, unit_cover=true)
        vel_total = vel_sp .+ vsc
        vel_par  = dot(vel_total, b_hat)
        vel_perp = norm(vel_total - vel_par * b_hat)

        # ---- 绘图 ----
        sp_title = col_idx == 1 ? sp[:name] : ""
        ax = Axis(fig_VDF[sp_idx, col_idx], title=sp_title, aspect=DataAspect(),
            xlabel=L"$V_{\parallel}$ (km/s)", ylabel=L"$V_{\perp}$ (km/s)",
            limits=(-plot_range, plot_range, 0, plot_range),
            xticks=[-1000, -500, 0, 500, 1000], yticks=[0, 500, 1000])

        if length(v_par) > 0
            colors_idx = map_color(df_vals, colorrange; scaler="log10")
            scatter!(ax, v_par, v_perp; color=colors_idx, colormap=:viridis,
                markersize=4, strokewidth=0)
        end

        # scatter!(ax, [vel_par], [vel_perp]; color=sp[:color], markersize=15, marker='x', strokewidth=3)

        # vlines!(ax, 0; color=:white, linestyle=:dash, linewidth=0.8)
        # hlines!(ax, 0; color=:white, linestyle=:dash, linewidth=0.8)

        if sp_idx == 1
            ax.title = Dates.format(swi_time, "HH:MM")
        end

        if col_idx > 1
            hideydecorations!(ax, grid=false, ticks=false)
        end
    end
end

# 全局色标
Colorbar(fig_VDF[1:n_species, n_times+1]; limits=colorrange,
    label=L"f(\mathbf{v})\ (\mathrm{s^3/km^6})",
    colormap=:viridis, scale=log10)

save(joinpath(VDF_path, fig_name), fig_VDF)
println("Saved figure: ", joinpath(VDF_path, fig_name))
