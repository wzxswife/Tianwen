using CommonDataFormat
using Dates, CSV, DataFrames
using Statistics, LinearAlgebra
using CairoMakie, DelaunayTriangulation
include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_plot.jl")
include("../../src/scripts/MAVEN_SWIA.jl")
include("../../src/scripts/MAVEN_STATIC.jl")
using .MAVEN_load
using .MAVEN_plot
using .MAVEN_SWIA
using .MAVEN_STATIC

date = DateTime(2022, 8, 24)
# date = DateTime(2015, 10, 29)
ts = DateTime.(date + Dates.Hour(7) + Dates.Minute(20))
# ts = DateTime.(date + Dates.Hour(11) + Dates.Minute(30) + Dates.Second(50))
time_stemp = ts .+ Dates.Minute(3) .* range(0, 3)

date_str = Dates.format.(time_stemp[1], "yyyymmdd")
date_str_dye = "$(year(time_stemp[1]))$(lpad(dayofyear(time_stemp[1]), 3, '0'))"

save_name = "MAVEN_SWIA_Ion_VDF_ParPerp_$(replace(string(time_stemp[1]), ":" => "-")).pdf"
save_path = joinpath(pwd(), "results", "IonBeam", save_name)
data_path = joinpath(pwd(), "Data", "MAVEN", "mvn_swi_l2_coarsesvy3d_$(date_str)_v02_r01.cdf")
quat_path = joinpath(pwd(), "Data", "MAVEN", "mvn_spice_swia_qu_$(date_str).csv")
mag_path = joinpath(pwd(), "Data", "MAVEN", "mvn_mag_l3_$(date_str_dye)ss1s_$(date_str)_v01_r01.f77_unformatted")

data = load_cdf(data_path)
quat_data = load_quat(quat_path)
quat_data[:data_load_flag] = true
swi_data = get_3dc!(data, quat_data = quat_data)
mag_data = load_mag_l3(mag_path)
time_b = mag_data[:epoch]
B_ss = mag_data[:B]
position_mso = mag_data[:position]

fig_VDF = Figure(size = (1200, 500))

spice = Dict(:mass_range => [0, 2], :m_int => 1)

colorrange = (1e0, 1e8)
spice_name = "H+"
V_30ev  = MAVEN_STATIC.ion_energy2v(30, spice[:m_int]) / 1000
V_494ev = MAVEN_STATIC.ion_energy2v(494, spice[:m_int]) / 1000
plot_range = 800

# 将 flux 值映射到颜色索引 1-256
function map_color(vals, cr; scaler="log10")
    cr_in = scaler == "log10" ? log10.(cr) : cr
    v_in  = scaler == "log10" ? log10.(vals) : vals
    mapped = round.(Int, ((v_in .- cr_in[1]) ./ (cr_in[2] - cr_in[1])) .* 255 .+ 1)
    clamp!(mapped, 1, 256)
    return mapped
end

for (col_idx, II) in enumerate(eachindex(time_stemp))
    println("Doing $(time_stemp[II]) $spice_name")
    _, time_ind = findmin(d -> abs(d - time_stemp[II]), swi_data[:epoch])
    swi_time = swi_data[:epoch][time_ind]

    _, time_b_ind = findmin(d -> abs(d - swi_time), time_b)
    local vsc = (position_mso[time_b_ind+1, :] .- position_mso[time_b_ind-1, :]) ./
                (datetime2unix(time_b[time_b_ind+1]) - datetime2unix(time_b[time_b_ind-1]))
    b0 = B_ss[time_b_ind, :]

    vel, _, _ = MAVEN_SWIA.v_3d(swi_data, time_ind; energy_range=[494, 1e8])
    vel = vel[1, :] .+ vsc

    # 速度网格（SWIA 16 阳极 × 4 偏转 × 48 能档）
    energy = reshape(swi_data[:energy], 1, 1, 48)
    theta  = reshape(swi_data[:theta][time_ind, :, :], 1, 4, 48)
    phi    = reshape(swi_data[:phi], 16, 1, 1)
    v0 = MAVEN_SWIA.ion_energy2v.(energy, 1) ./ 1000
    vv = vec(MAVEN_SWIA.sphere2xyz_for_SWIA.(v0, theta, phi))
    v = hcat([x[1] for x in vv], [x[2] for x in vv], [x[3] for x in vv])

    dF = vec(swi_data[:diff_en_fluxes][time_ind, :, :, :]) .* 1e-3

    # ========== 沿磁场分解速度 ==========
    b_hat = b0 ./ norm(b0)
    v_par  = v * b_hat                                           # V∥
    v_perp = vec(sqrt.(sum((v .- v_par * b_hat').^2, dims=2)))    # V⟂

    # 整体漂移速度在磁场坐标系下
    vel_par  = dot(vel, b_hat)
    vel_perp = norm(vel - vel_par * b_hat)

    # 剔除无效点
    valid = dF .> 1e-20
    v_par  = v_par[valid]
    v_perp = v_perp[valid]
    dF     = dF[valid]

    title = Dates.format(swi_time, "HH:MM:SS")
    ax = Axis(fig_VDF[1, col_idx], title=title, aspect=DataAspect(),
        xlabel=L"$V_{\parallel}$ (km/s)", ylabel=L"$V_{\perp}$ (km/s)",
        limits=(-plot_range, plot_range, 0, plot_range),
        xticks=[-600, -300, 0, 300, 600], yticks=[0, 300, 600])

    # # Voronoi 着色（与 VDF_2d_slip 相同方法）
    # pts = [Point2f(v_par[i], v_perp[i]) for i in eachindex(v_par)]
    # if length(pts) > 3
    #     tri = voronoi(triangulate(pts))
    #     colors_idx = map_color(dF, colorrange; scaler="log10")
    #     voronoiplot!(ax, tri, color=colors_idx, colormap=:viridis,
    #         strokewidth=0, markersize=0)
    # end
    
    # 散点图，颜色表示通量大小
    colors_idx = map_color(dF, colorrange; scaler="log10")
    scatter!(ax, v_par, v_perp; color=colors_idx, colormap=:viridis,
        markersize=6, strokewidth=0)

    # 能壳（等能量圆）
    lines!(ax, Circle(Point2f(0, 0), V_30ev),  color=:white, linewidth=1.5)
    lines!(ax, Circle(Point2f(0, 0), V_494ev), color=:white, linewidth=1.5, linestyle=:dash)

    # 漂移速度（整体运动）
    scatter!(ax, [vel_par], [vel_perp], color=:red, markersize=18, marker='x', strokewidth=3)

    # 参考线
    vlines!(ax, 0, color=:white, linestyle=:dash, linewidth=0.8)
    hlines!(ax, 0, color=:white, linestyle=:dash, linewidth=0.8)

    if col_idx > 1
        hideydecorations!(ax, grid=false, ticks=false)
    end
end

Colorbar(fig_VDF[1, length(time_stemp)+1], limits=colorrange,
    label=rich("i", superscript("+"), " Flux (cm", superscript("-2"),
        " s", superscript("-1"), " sr", superscript("-1)")),
    colormap=:viridis, scale=log10)

save(save_path, fig_VDF)
println("Saved: ", save_path)
