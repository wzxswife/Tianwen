"""
STATIC H+ VDF 散点图 (坐标系: x̂=B̂, ŷ=V̂sw, ẑ=B̂×V̂sw) — 逐帧

逐帧版: 遍历时间窗口内每个时刻, 将所有帧的测量点画在同一张散点图上。
坐标系 (B̂, V̂sw) 用窗口平均值作为固定参考帧, 不随帧变化。

数据处理:
    1. 平均磁场经四元数旋转到 STATIC 仪器系 → B̂
    2. 太阳风速度: 平均 eflux/E² 峰值能道
    3. 构建基矢量: x̂=B̂, ŷ=V̂sw, ẑ=B̂×V̂sw (STATIC 系)
    4. 逐帧: 每个时刻的 64 方向 × 32 能道 bin → 速度矢量 → 投影 → 散点

注意:
    - 速度均为 STATIC 仪器系 (未加飞船速度), 与 vdf_3D.jl 一致
    - 每帧的 swp_ind 可能不同, 方向/速度大小随 swp 变化
"""

using Glob
using Dates
using LinearAlgebra
using Statistics
using Unitful
using ThreadsX
using Rotations
using CoordinateTransformations
using CairoMakie
using LaTeXStrings
using DelaunayTriangulation
DelaunayTriangulation.toggle_warn_on_dupes!()
import PhysicalConstants.CODATA2018: m_p, c_0

include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_STATIC.jl")
include("../../src/scripts/SphericalSplines.jl")
using .MAVEN_STATIC
using .MAVEN_load
using .SphericalSplines

# ============================================================
# 参数设置
# ============================================================

# --- 时间范围 ---
date = Date(2022, 8, 24)
ts = DateTime(date) + Dates.Hour(7) + Dates.Minute(27)
te = ts + Dates.Minute(3)
# --- 指定绘图时间点 (HH:MM:SS 格式), 留空则绘制全部帧 ---
plot_times = [
    "07:27:20",
    "07:27:40",
    "07:28:00",
    "07:28:20",
    "07:28:40",
    "07:29:00",
]

# --- 太阳风速度估算参数 (eflux/E² 峰值搜寻 bin) ---
sw_i = 2
sw_j = 14

# --- 过滤参数 ---
flux_threshold = 10^5   # eflux 阈值 (eV/cm²/s/sr/eV)
plot_range = (-1000, 1000)  # 绘图速度范围 (km/s)

# --- 能量排除范围 (eV): 设置为 nothing 则不排除 ---
energy_exclude_min = nothing   # 排除能量下限 (eV), nothing 表示不限
energy_exclude_max = nothing # 排除能量上限 (eV), nothing 表示不限

# --- 热力图参数 ---
heatmap_nbin = 100  # 热力图网格分辨率


# ============================================================
# 路径构造与数据加载
# ============================================================
datestr = Dates.format(date, "yyyymmdd")
year_num = Dates.format(date, "yyyy")
month_num = Dates.format(date, "mm")

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

MVN_mag2c1hz = MAVEN_load.load_mag_l2(MAG_file[1])
MVN_static = MAVEN_load.load_cdf(STA_file[1])
println("Data loaded.")

# ============================================================
# MVN MAG: 时间窗口内平均磁场 (MSO)
# ============================================================
MVN_BMSO1hz = MVN_mag2c1hz[:B]
ind = ThreadsX.findall(ts .<= MVN_mag2c1hz[:epoch] .<= te)
MVN_Bmeanmso = dropdims(mean(Matrix(MVN_BMSO1hz[ind, :]), dims=(1,)), dims=(1,))

# ============================================================
# MVN STATIC: 数据提取
# ============================================================
MVN_staticswp = MVN_static[:swp_ind] |> Array
MVN_staticswp .+= 1
MVN_statictime = MVN_static[:epoch] |> Array
MVN_staticeflux = MVN_static[:eflux] |> Array
MVN_staticmass = MVN_static[:mass_arr] |> Array
MVN_staticenergy = MVN_static[:energy] |> Array
MVN_statictheta = MVN_static[:theta] |> Array # -45 to 45 deg, 纬度
MVN_staticphi = MVN_static[:phi] |> Array     # -180 to 180 deg, 经度
MVN_staticvelocity = sqrt.(MVN_staticenergy .* (2.0 / ustrip(u"eV", m_p * c_0^2)) ./ MVN_staticmass) .* ustrip(u"km/s", c_0)

MVN_staticquat_mso = MVN_static[:quat_mso] |> Array # STATIC→MSO 四元数

MVN_staticind = ThreadsX.findall(ts .<= MVN_statictime .<= te)
MVN_staticswpind = MVN_staticswp[MVN_staticind[1]] # 窗口首帧扫描档

# 平均 eflux (仅用于估算太阳风速度, 不用于绘图)
MVN_staticefluxmean = dropdims(mean(MVN_staticeflux[MVN_staticind, :, :, :], dims=(1,)), dims=(1,))

# 平均时刻旋转: 磁场从 MSO 转回 STATIC 仪器系
MVN_staticindmean = round(Int64, mean(MVN_staticind))
MVN_static2msomean = QuatRotation(MVN_staticquat_mso[MVN_staticindmean, :])
MVN_Bmeanstatic = MVN_static2msomean' * MVN_Bmeanmso
MVN_parastatic = MVN_Bmeanstatic ./ norm(MVN_Bmeanstatic) # B̂

# ============================================================
# 太阳风速度 (STATIC 系): 平均 eflux/E² 峰值能道
# ============================================================
ind = argmax(MVN_staticefluxmean[1, (sw_j - 1) * 4 + sw_i, :] ./ MVN_staticenergy[1, 1, :, MVN_staticswpind, 1].^2)
MVN_staticvswrough = CartesianFromSpherical()(
    CoordinateTransformations.Spherical(
        MVN_staticvelocity[1, (sw_j - 1) * 4 + sw_i, ind, MVN_staticswpind, 1],
        deg2rad(MVN_staticphi[1, (sw_j - 1) * 4 + sw_i, ind, MVN_staticswpind, 1]),
        deg2rad(MVN_statictheta[1, (sw_j - 1) * 4 + sw_i, ind, MVN_staticswpind, 1])))
# println("太阳风速度 (STATIC系): ", round.(MVN_staticvswrough, digits = 1), " km/s, |Vsw| = ", round(norm(MVN_staticvswrough), digits = 1), " km/s")

# ============================================================
# 构建 MVA 坐标系: x̂=B̂, ŷ=B̂×V̂sw, ẑ=B̂×(B̂×V̂sw) (STATIC 系)
# ============================================================
e_x = MVN_parastatic                                                    # B̂
e_y = cross(e_x, MVN_staticvswrough ./ norm(MVN_staticvswrough))       # B̂ × V̂sw
e_y = e_y ./ norm(e_y)
e_z = cross(e_x, e_y)                                                   # B̂ × (B̂ × V̂sw)
e_z = e_z ./ norm(e_z)

# 基矢量转换到 MSO 坐标系 (用于热力图坐标转换)
e_x_mso = MVN_static2msomean * e_x  # B̂ in MSO
e_y_mso = MVN_static2msomean * e_y  # B̂×V̂sw in MSO
e_z_mso = MVN_static2msomean * e_z  # B̂×(B̂×V̂sw) in MSO

# ============================================================
# 弓激波法向计算 (Edberg et al., 2008 模型)
# ============================================================
const Rm = 3390.0  # 火星半径 (km)
const BS_XF = 0.55; const BS_EPSILON = 1.05; const BS_L = 2.10

function bowshock_r(x)
    temp = (BS_EPSILON^2 - 1.0) * (x - BS_XF)^2 -
           2 * BS_EPSILON * BS_L * (x - BS_XF) + BS_L^2
    return temp >= 0 ? sqrt(temp) : Inf64
end

function bowshock_dr_dx(x)
    r_bs = bowshock_r(x)
    (!isfinite(r_bs) || r_bs < 1e-10) && return Inf64
    numerator = (BS_EPSILON^2 - 1.0) * (x - BS_XF) - BS_EPSILON * BS_L
    return numerator / r_bs
end

function bowshock_normal(x, y, z)
    r_perp = sqrt(y^2 + z^2)
    r_bs = bowshock_r(x)
    (!isfinite(r_bs)) && return [0.0, 0.0, 0.0]
    drdx = bowshock_dr_dx(x)
    (!isfinite(drdx)) && return [0.0, 0.0, 0.0]
    if r_perp < 1e-10
        nx = drdx > 0 ? 1.0 : -1.0
        return [nx, 0.0, 0.0]
    end
    nx = -drdx; ny = y / r_perp; nz = z / r_perp
    n = [nx, ny, nz] / norm([nx, ny, nz])
    n[1] > 0 && (n = -n)
    return n
end

# 计算时间窗口中心时刻的位置和激波法向
MVN_position = MVN_mag2c1hz[:position]
mid_idx = ind[div(length(ind), 2) + 1]
pos_mean_mso = Float64.(MVN_position[mid_idx, :])
pos_mean_Rm = pos_mean_mso ./ Rm

shock_normal_mso = bowshock_normal(pos_mean_Rm[1], pos_mean_Rm[2], pos_mean_Rm[3])
shock_surface_mso = cross(shock_normal_mso, [0.0, 0.0, 1.0])  # 激波面切线 (近似)
if norm(shock_surface_mso) < 1e-10
    shock_surface_mso = cross(shock_normal_mso, [0.0, 1.0, 0.0])
end
shock_surface_mso = shock_surface_mso / norm(shock_surface_mso)

# 转换到 (B̂, V̂sw, B̂×V̂sw) 坐标系
# 由于基矢量不正交, 使用 MSO→(B̂, V̂sw, B̂×V̂sw) 的转换矩阵
# [e_x_mso, e_y_mso, e_z_mso] 构成变换矩阵
Basis_mso = hcat(e_x_mso, e_y_mso, e_z_mso)
shock_normal_bvsw = Basis_mso' * shock_normal_mso
shock_surface_bvsw = Basis_mso' * shock_surface_mso

# 在 xz 和 yz 平面的投影 (仅保留平面内的分量)
shock_xz_dir = [shock_surface_bvsw[1], shock_surface_bvsw[3]]  # (B̂, B̂×V̂sw) 平面
shock_yz_dir = [shock_surface_bvsw[2], shock_surface_bvsw[3]]  # (V̂sw, B̂×V̂sw) 平面
if norm(shock_xz_dir) > 1e-10
    shock_xz_dir = shock_xz_dir / norm(shock_xz_dir)
end
if norm(shock_yz_dir) > 1e-10
    shock_yz_dir = shock_yz_dir / norm(shock_yz_dir)
end

# ============================================================
# 逐帧遍历: 每个时刻 × 每个 bin → 速度矢量 + eflux
# ============================================================
n_frames = length(MVN_staticind)
n_bins_per_frame = 64 * 32

# 按帧存储数据
V_frames = Vector{Matrix{Float64}}(undef, n_frames)
eflux_frames = Vector{Vector{Float64}}(undef, n_frames)
energy_frames = Vector{Vector{Float64}}(undef, n_frames)
time_labels = Vector{String}(undef, n_frames)

for (fi, t_idx) in enumerate(MVN_staticind)
    swp = MVN_staticswp[t_idx]
    V_frame = Array{Float64}(undef, n_bins_per_frame, 3)
    eflux_frame = Array{Float64}(undef, n_bins_per_frame)
    energy_frame = Array{Float64}(undef, n_bins_per_frame)
    for (bin_idx, (ij, k)) in enumerate(Iterators.product(1:64, 1:32))
        dir = CartesianFromSpherical()(
            CoordinateTransformations.Spherical(1.0,
                deg2rad(MVN_staticphi[1, ij, k, swp, 1]),
                deg2rad(MVN_statictheta[1, ij, k, swp, 1])))
        V_frame[bin_idx, :] = MVN_staticvelocity[1, ij, k, swp, 1] .* Array(dir)
        eflux_frame[bin_idx] = MVN_staticeflux[t_idx, 1, ij, k]
        energy_frame[bin_idx] = MVN_staticenergy[1, ij, k, swp, 1]
    end
    V_frames[fi] = V_frame
    eflux_frames[fi] = eflux_frame
    energy_frames[fi] = energy_frame
    time_labels[fi] = Dates.format(MVN_statictime[t_idx], "HH:MM:SS")
end

# ============================================================
# 选择指定时间点最近的帧
# ============================================================
if isempty(plot_times)
    selected_fi = collect(1:n_frames)
else
    selected_fi = Int[]
    t_start = MVN_statictime[MVN_staticind[1]]
    t_end = MVN_statictime[MVN_staticind[end]]
    for pt in plot_times
        parts = split(pt, ":")
        pt_dt = DateTime(date) + Dates.Hour(parse(Int, parts[1])) +
                Dates.Minute(parse(Int, parts[2])) + Dates.Second(parse(Int, parts[3]))
        diffs = abs.(MVN_statictime[MVN_staticind] .- pt_dt)
        _, closest = findmin(diffs)
        if pt_dt < t_start || pt_dt > t_end
            @warn "时刻 $pt 超出数据范围 ($(Dates.format(t_start,"HH:MM:SS")) ~ $(Dates.format(t_end,"HH:MM:SS"))), 使用最近帧 $(time_labels[closest])"
        end
        push!(selected_fi, closest)
    end
    selected_fi = unique(sort(selected_fi))
end
n_selected = length(selected_fi)
println("n_frames=$n_frames, selected_fi=$selected_fi, n_selected=$n_selected")
println("Time range: ", MVN_statictime[MVN_staticind[1]], " ~ ", MVN_statictime[MVN_staticind[end]])

# ============================================================
# 投影到 MVA 坐标系: (v·B̂, v·(B̂×V̂sw), v·(B̂×(B̂×V̂sw)))
# ============================================================
Basis = hcat(e_x, e_y, e_z)
proj_frames = [V_frames[fi] * Basis for fi in 1:n_frames]

# ============================================================
# 绘图 (CairoMakie 散点): 2行n列子图, 每列一个时间点
# ============================================================
CairoMakie.activate!()

fig_name = "MAVEN_STATIC_H+_VDF_MVA_xz_yz_frame_" *
    Dates.format(ts, "yyyymmddTHHMM") * "-" *
    Dates.format(te, "HHMM") * ".pdf"

crange = (flux_threshold, 1e9)
ticks = [-1000, -500, 0, 500, 1000]

fig = Figure(resolution = (400 * n_selected, 1200), fontsize = 12)

for (col, fi) in enumerate(selected_fi)
    proj_fi = proj_frames[fi]
    eflux_fi = eflux_frames[fi]
    energy_fi = energy_frames[fi]

    # 逐帧过滤
    e_keep = (energy_exclude_min === nothing || energy_exclude_max === nothing) ? trues(n_bins_per_frame) :
        .!(energy_exclude_min .<= energy_fi .<= energy_exclude_max)
    f_keep = (eflux_fi .> flux_threshold) .& e_keep .&
        (proj_fi[:, 1] .> plot_range[1]) .& (proj_fi[:, 1] .< plot_range[2]) .&
        (proj_fi[:, 3] .> plot_range[1]) .& (proj_fi[:, 3] .< plot_range[2])

    # ---- xz 平面 (第1行) ----
    ax_xz = Axis(fig[1, col],
        xlabel = col == 1 ? L"$v \cdot \hat{B}$ (km/s)" : "",
        ylabel = L"$v \cdot (\hat{B} \times (\hat{B} \times \hat{V}_{\mathrm{sw}}))$ (km/s)",
        title = time_labels[fi],
        xticks = ticks, yticks = ticks)
    xlims!(ax_xz, plot_range[1], plot_range[2])
    ylims!(ax_xz, plot_range[1], plot_range[2])
    hlines!(ax_xz, 0, color = :white, linestyle = :dash)
    vlines!(ax_xz, 0, color = :white, linestyle = :dash)
    if count(f_keep) > 0
        scatter!(ax_xz, proj_fi[f_keep, 1], proj_fi[f_keep, 3];
            color = eflux_fi[f_keep], colormap = :gist_earth, colorscale = log10,
            colorrange = crange, markersize = 4)
    end
    if col == 1 && norm(shock_xz_dir) > 1e-10
        L = Float64(plot_range[2]) * 0.9
        lines!(ax_xz, [-L*shock_xz_dir[1], L*shock_xz_dir[1]],
                    [-L*shock_xz_dir[2], L*shock_xz_dir[2]],
                    color=:black, linestyle=:dash, linewidth=2)
    end

    # ---- yz 平面 (第2行) ----
    f_keep_yz = (eflux_fi .> flux_threshold) .& e_keep .&
        (proj_fi[:, 2] .> plot_range[1]) .& (proj_fi[:, 2] .< plot_range[2]) .&
        (proj_fi[:, 3] .> plot_range[1]) .& (proj_fi[:, 3] .< plot_range[2])

    ax_yz = Axis(fig[2, col],
        xlabel = L"$v \cdot (\hat{B} \times \hat{V}_{\mathrm{sw}})$ (km/s)",
        ylabel = col == 1 ? L"$v \cdot (\hat{B} \times (\hat{B} \times \hat{V}_{\mathrm{sw}}))$ (km/s)" : "",
        xticks = ticks, yticks = ticks)
    xlims!(ax_yz, plot_range[1], plot_range[2])
    ylims!(ax_yz, plot_range[1], plot_range[2])
    hlines!(ax_yz, 0, color = :white, linestyle = :dash)
    vlines!(ax_yz, 0, color = :white, linestyle = :dash)
    if count(f_keep_yz) > 0
        scatter!(ax_yz, proj_fi[f_keep_yz, 2], proj_fi[f_keep_yz, 3];
            color = eflux_fi[f_keep_yz], colormap = :gist_earth, colorscale = log10,
            colorrange = crange, markersize = 4)
    end
    if col == 1 && norm(shock_yz_dir) > 1e-10
        L = Float64(plot_range[2]) * 0.9
        lines!(ax_yz, [-L*shock_yz_dir[1], L*shock_yz_dir[1]],
                    [-L*shock_yz_dir[2], L*shock_yz_dir[2]],
                    color=:black, linestyle=:dash, linewidth=2)
    end

    # ---- xy 平面 (第3行) ----
    f_keep_xy = (eflux_fi .> flux_threshold) .& e_keep .&
        (proj_fi[:, 1] .> plot_range[1]) .& (proj_fi[:, 1] .< plot_range[2]) .&
        (proj_fi[:, 2] .> plot_range[1]) .& (proj_fi[:, 2] .< plot_range[2])

    ax_xy = Axis(fig[3, col],
        xlabel = col == 1 ? L"$v \cdot \hat{B}$ (km/s)" : "",
        ylabel = col == 1 ? L"$v \cdot (\hat{B} \times \hat{V}_{\mathrm{sw}})$ (km/s)" : "",
        xticks = ticks, yticks = ticks)
    xlims!(ax_xy, plot_range[1], plot_range[2])
    ylims!(ax_xy, plot_range[1], plot_range[2])
    hlines!(ax_xy, 0, color = :white, linestyle = :dash)
    vlines!(ax_xy, 0, color = :white, linestyle = :dash)
    if count(f_keep_xy) > 0
        scatter!(ax_xy, proj_fi[f_keep_xy, 1], proj_fi[f_keep_xy, 2];
            color = eflux_fi[f_keep_xy], colormap = :gist_earth, colorscale = log10,
            colorrange = crange, markersize = 4)
    end
end

Colorbar(fig[1:3, n_selected+1], limits = crange, colormap = :gist_earth, scale = log10,
    label = L"H$^{+}$ Differential Energy Flux (eV/cm²/s/sr/eV)")

save(joinpath(VDF_path, fig_name), fig)
println("Saved figure: ", joinpath(VDF_path, fig_name))

# ============================================================
# 热力图: 球面样条插值 — 逐帧版
# 每个选中的时间点单独构建球面样条 → 单独插值 → 一格一个时间点
# (与散点图一致: 2 行 × n_selected 列, 每列一个时间点)
# ============================================================

n_energy_ch = size(MVN_staticeflux, 4)  # 32 能道
n_dir = size(MVN_staticeflux, 3)        # 64 方向

# 在 (B̂, V̂sw, B̂×V̂sw) 坐标系下构建球面样条插值热力图 (单帧版本)
# plane: "xz" → (v·B̂, v·(B̂×V̂sw)), "yz" → (v·V̂sw, v·(B̂×V̂sw))
# splines 已按 v_ch 升序排列, 与 build_vdf_grid_bvsw 用法一致
function build_vdf_grid_bvsw(splines, v_ch; n_grid=80, plot_range=(-1000, 1000), plane="xz",
                              e_x_mso=nothing, e_y_mso=nothing, e_z_mso=nothing)
    n_energy = length(v_ch)
    v_min, v_max = extrema(v_ch)
    grid_log = fill(NaN, n_grid, n_grid)

    # 网格坐标
    v_edges_grid = range(plot_range[1], plot_range[2], length=n_grid + 1)
    v_centers = collect((v_edges_grid[1:end-1] .+ v_edges_grid[2:end]) ./ 2)

    # 逐网格点插值 (参考 sta_VDF_eflux_interp.jl 第 235-271 行)
    for ix in 1:n_grid, iy in 1:n_grid
        v1 = v_centers[ix]
        v2 = v_centers[iy]

        # 构建 3D 方向 (第三分量=0, 在 2D 平面上)
        dir_mso = zeros(3)
        if plane == "xz"
            dir_mso[1] = v1   # v·B̂
            dir_mso[3] = v2   # v·(B̂×(B̂×V̂sw))
        elseif plane == "yz"
            dir_mso[2] = v1   # v·(B̂×V̂sw)
            dir_mso[3] = v2   # v·(B̂×(B̂×V̂sw))
        else  # "xy"
            dir_mso[1] = v1   # v·B̂
            dir_mso[2] = v2   # v·(B̂×V̂sw)
        end

        speed_pt = norm(dir_mso)
        speed_pt > 0 || continue
        dir_mso ./= speed_pt  # 归一化方向 (当前在 B̂,V̂sw,B̂×V̂sw 坐标系)

        # 坐标转换: (B̂,V̂sw,B̂×V̂sw) → MSO
        if e_x_mso !== nothing
            dir_bvsw = copy(dir_mso)
            dir_mso[:] = dir_bvsw[1] * e_x_mso + dir_bvsw[2] * e_y_mso + dir_bvsw[3] * e_z_mso
            dir_mso ./= norm(dir_mso)  # 重新归一化
        end

        # 速度范围检查 (必须在数据速度范围内)
        (v_min ≤ speed_pt ≤ v_max) || continue

        # 找到包围 speed_pt 的两个速度通道
        i_high = searchsortedfirst(v_ch, speed_pt)
        i_low = i_high - 1
        (i_low ≥ 1 && i_high ≤ n_energy) || continue

        S_lo = splines[i_low]
        S_hi = splines[i_high]
        (S_lo === nothing || S_hi === nothing) && continue

        # MSO 球面样条插值: 输入 MSO 方向
        log_flux_lo = S_lo(Vector(dir_mso))
        log_flux_hi = S_hi(Vector(dir_mso))
        (isnan(log_flux_lo) || isnan(log_flux_hi)) && continue

        # log-log 速度通道插值 (与 sta_VDF_eflux_interp.jl 第 263-270 行一致)
        log_v = log10(speed_pt)
        log_v_lo = log10(v_ch[i_low])
        log_v_hi = log10(v_ch[i_high])
        t = (log_v - log_v_lo) / (log_v_hi - log_v_lo)

        grid_log[ix, iy] = log_flux_lo + t * (log_flux_hi - log_flux_lo)
    end

    return v_centers, grid_log
end

# 单帧: 构建球面样条 (该帧各能道的方向旋转到 MSO 后 → log10(eflux)),
# 返回按 v_ch 升序排列的 splines 与 v_ch
function build_frame_splines(t_idx)
    swp = MVN_staticswp[t_idx]
    Q_f = QuatRotation(MVN_staticquat_mso[t_idx, :])

    # 该帧各能道平均速度
    v_ch_f = zeros(n_energy_ch)
    for k in 1:n_energy_ch
        speed_k = MVN_staticvelocity[1, :, k, swp, 1]
        svalid = .!isnan.(speed_k)
        v_ch_f[k] = svalid == 0 ? 0.0 : mean(speed_k[svalid])
    end

    # 逐能道构建球面样条 (方向在 MSO 系)
    splines_f = Vector{Any}(undef, n_energy_ch)
    for k in 1:n_energy_ch
        splines_f[k] = nothing
        theta_k = MVN_statictheta[1, :, k, swp, 1]
        phi_k = MVN_staticphi[1, :, k, swp, 1]
        flux_k = MVN_staticeflux[t_idx, 1, :, k]

        valid = .!isnan.(theta_k) .& .!isnan.(phi_k) .& .!isnan.(flux_k) .& (flux_k .> 0)
        nv = count(valid)
        nv ≥ 4 || continue

        ct = cosd.(theta_k[valid]); st = sind.(theta_k[valid])
        cp = cosd.(phi_k[valid]);   sp = sind.(phi_k[valid])
        dirs_SC = hcat(ct.*cp, ct.*sp, st)

        # SC → MSO
        dirs_MSO = (Q_f * dirs_SC')'          # 每行一个 MSO 方向
        dirs_MSO ./= vec(sqrt.(sum(abs2, dirs_MSO, dims=2)))

        splines_f[k] = InterpolationSpline(dirs_MSO, log10.(flux_k[valid]))
    end

    # 按速度排序 (低→高), 与 build_vdf_grid_bvsw 的 i_low/i_high 用法一致
    sort_idx = sortperm(v_ch_f)
    v_ch_f = v_ch_f[sort_idx]
    splines_f = splines_f[sort_idx]
    return splines_f, v_ch_f
end

println("逐帧球面样条插值生成热力图...")

# 先对每个选中帧计算网格, 并收集全局颜色范围 (便于帧间对比)
all_valid = Float64[]
frame_grids = Vector{Tuple{Vector{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64},Vector{Float64},Matrix{Float64}}}(undef, n_selected)
for (col, fi) in enumerate(selected_fi)
    t_idx = MVN_staticind[fi]
    splines_f, v_ch_f = build_frame_splines(t_idx)
    xzc, xzg = build_vdf_grid_bvsw(splines_f, v_ch_f; n_grid=heatmap_nbin,
        plot_range=plot_range, plane="xz", e_x_mso=e_x_mso, e_y_mso=e_y_mso, e_z_mso=e_z_mso)
    yzc, yzg = build_vdf_grid_bvsw(splines_f, v_ch_f; n_grid=heatmap_nbin,
        plot_range=plot_range, plane="yz", e_x_mso=e_x_mso, e_y_mso=e_y_mso, e_z_mso=e_z_mso)
    xyc, xyg = build_vdf_grid_bvsw(splines_f, v_ch_f; n_grid=heatmap_nbin,
        plot_range=plot_range, plane="xy", e_x_mso=e_x_mso, e_y_mso=e_y_mso, e_z_mso=e_z_mso)
    frame_grids[col] = (xzc, xzg, yzc, yzg, xyc, xyg)
    append!(all_valid, filter(isfinite, vec(xzg)))
    append!(all_valid, filter(isfinite, vec(yzg)))
    append!(all_valid, filter(isfinite, vec(xyg)))
end

if length(all_valid) > 10
    lo, hi = quantile(all_valid, [0.05, 0.95])
    cr_log = (lo, hi)
else
    cr_log = (6.0, 8.0)  # fallback
end

# ---- 太阳风方向 V̂sw 在 MVA 基 (B̂, B̂×V̂sw, B̂×(B̂×V̂sw)) 下的分量 ----
vsw_hat = MVN_staticvswrough ./ norm(MVN_staticvswrough)   # V̂sw (STATIC 系)
vsw_x = dot(vsw_hat, e_x)   # V̂sw · B̂
vsw_y = dot(vsw_hat, e_y)   # V̂sw · (B̂×V̂sw)
vsw_z = dot(vsw_hat, e_z)   # V̂sw · (B̂×(B̂×V̂sw))

ticks = [-1000, -500, 0, 500, 1000]
title_heat = "MAVEN STATIC H⁺ VDF $(datestr) $(Dates.format(ts,"HH:MM:SS"))-$(Dates.format(te,"HH:MM:SS")) (MVA) per-frame spline"
colorbar_label = L"$\log_{10}$ Differential Energy Flux (eV/cm$^2$/s/sr/eV)"

fig_heat_name = "MAVEN_STATIC_H+_VDF_MVA_xz_yz_heatmap_frame_" *
    Dates.format(ts, "yyyymmddTHHMM") * "-" *
    Dates.format(te, "HHMM") * ".pdf"

fig2 = Figure(resolution = (600 * n_selected + 80, 1600), fontsize = 12)

for (col, fi) in enumerate(selected_fi)
    xzc, xzg, yzc, yzg, xyc, xyg = frame_grids[col]

    # ---- xz 平面热力图 (第1行: 该时间点) ----
    ax_xz = Axis(fig2[1, col], aspect = 1, yticks = ticks, xticks = ticks,
        limits = (plot_range, plot_range), alignmode = Inside(),
        title = time_labels[fi])
    heatmap!(ax_xz, xzc, xzc, xzg; colormap = :viridis,
        colorrange = cr_log, nan_color = :transparent)
    hlines!(ax_xz, 0, linestyle = :dash, color = :white)
    vlines!(ax_xz, 0, linestyle = :dash, color = :white)
    ax_xz.xlabel = col == 1 ? L"$v \cdot \hat{B}$ (km/s)" : ""
    ax_xz.ylabel = L"$v \cdot (\hat{B} \times (\hat{B} \times \hat{V}_{\mathrm{sw}}))$ (km/s)"

    # 激波面方向 (黑色虚线)
    if norm(shock_xz_dir) > 1e-10
        Lz = Float64(plot_range[2]) * 0.9
        lines!(ax_xz, [-Lz*shock_xz_dir[1], Lz*shock_xz_dir[1]],
                    [-Lz*shock_xz_dir[2], Lz*shock_xz_dir[2]],
                    color=:black, linestyle=:dash, linewidth=2)
    end
    # # 太阳风方向 V̂sw (绿色实线) — 投影到 xz 平面: (vsw_x, vsw_z)
    # Lv = Float64(plot_range[2]) * 0.9
    # lines!(ax_xz, [-Lv*vsw_x, Lv*vsw_x], [-Lv*vsw_z, Lv*vsw_z],
    #     color = :green, linewidth = 2)

    # ---- yz 平面热力图 (第2行: 该时间点) ----
    ax_yz = Axis(fig2[2, col], aspect = 1, yticks = ticks, xticks = ticks,
        limits = (plot_range, plot_range), alignmode = Inside())
    heatmap!(ax_yz, yzc, yzc, yzg; colormap = :viridis,
        colorrange = cr_log, nan_color = :transparent)
    hlines!(ax_yz, 0, linestyle = :dash, color = :white)
    vlines!(ax_yz, 0, linestyle = :dash, color = :white)
    ax_yz.xlabel = L"$v \cdot (\hat{B} \times \hat{V}_{\mathrm{sw}})$ (km/s)"
    ax_yz.ylabel = col == 1 ? L"$v \cdot (\hat{B} \times (\hat{B} \times \hat{V}_{\mathrm{sw}}))$ (km/s)" : ""

    # 激波面方向 (黑色虚线)
    if norm(shock_yz_dir) > 1e-10
        Lz = Float64(plot_range[2]) * 0.9
        lines!(ax_yz, [-Lz*shock_yz_dir[1], Lz*shock_yz_dir[1]],
                    [-Lz*shock_yz_dir[2], Lz*shock_yz_dir[2]],
                    color=:black, linestyle=:dash, linewidth=2)
    end
    # # 太阳风方向 V̂sw (绿色实线) — 投影到 yz 平面: (vsw_y, vsw_z)
    # lines!(ax_yz, [-Lv*vsw_y, Lv*vsw_y], [-Lv*vsw_z, Lv*vsw_z],
    #     color = :green, linewidth = 2)

    # ---- xy 平面热力图 (第3行: 该时间点) ----
    ax_xy = Axis(fig2[3, col], aspect = 1, yticks = ticks, xticks = ticks,
        limits = (plot_range, plot_range), alignmode = Inside())
    heatmap!(ax_xy, xyc, xyc, xyg; colormap = :viridis,
        colorrange = cr_log, nan_color = :transparent)
    hlines!(ax_xy, 0, linestyle = :dash, color = :white)
    vlines!(ax_xy, 0, linestyle = :dash, color = :white)
    ax_xy.xlabel = col == 1 ? L"$v \cdot \hat{B}$ (km/s)" : ""
    ax_xy.ylabel = col == 1 ? L"$v \cdot (\hat{B} \times \hat{V}_{\mathrm{sw}})$ (km/s)" : ""
end

Colorbar(fig2[1:3, n_selected+1], limits = cr_log, colormap = :viridis,
    label = colorbar_label)
Label(fig2[0, :], title_heat, fontsize = 18, tellwidth = false)

save(joinpath(VDF_path, fig_heat_name), fig2)
println("Saved heatmap: ", joinpath(VDF_path, fig_heat_name))

# ============================================================
# Voronoi 图: (B̂, V̂sw, B̂×V̂sw) 坐标系
# 参考 sta_VDF.jl 的 VDF_2d_slip 绘图逻辑
# ============================================================
println("\n生成 Voronoi 图...")

# Voronoi 图参数
voronoi_colorrange = (10^3, 10^9)  # eflux 颜色范围

# Gram-Schmidt 正交化: 从 e_x(B̂) 和 e_y(V̂sw) 构建正交基
# 参考 sta_VDF.jl 的 slice2d_cal_rot 函数
function build_orthonormal_basis(v1, v2)
    a = normalize(v1)
    c = cross(a, normalize(v2))
    c = normalize(c)
    b = cross(c, a)
    b = normalize(b)
    return a, b, c  # 正交基: (B̂_orth, perp, B̂×V̂sw)
end

# xz 平面基: x̂=B̂, ẑ=B̂×(B̂×V̂sw)
basis_xz_a, basis_xz_b, basis_xz_c = build_orthonormal_basis(e_x, e_z)
rot_xz = hcat(basis_xz_a, basis_xz_b, basis_xz_c)'

# yz 平面基: ŷ=B̂×V̂sw, ẑ=B̂×(B̂×V̂sw)
basis_yz_a, basis_yz_b, basis_yz_c = build_orthonormal_basis(e_y, e_z)
rot_yz = hcat(basis_yz_a, basis_yz_b, basis_yz_c)'

# 合并所有帧数据用于 Voronoi
V_all = reduce(vcat, V_frames)
eflux_all = reduce(vcat, eflux_frames)

# 旋转到目标坐标系
V_rot_xz = V_all * rot_xz'
V_rot_yz = V_all * rot_yz'

# 颜色映射 (eflux → 1-256 索引, 处理零/负值)
function color_map_voronoi(c, cr)
    c_safe = copy(c)
    c_safe[c_safe .<= 0] .= 1e-20  # 零/负值设为极小正值
    c_log = log10.(c_safe)
    cr_log = log10.(cr)
    colors = round.(Int, ((c_log .- cr_log[1]) ./ (cr_log[2] - cr_log[1])) .* 255 .+ 1)
    colors[colors .> 256] .= 256
    colors[colors .< 1] .= 1
    return colors
end

# 创建 Voronoi 图
fig_voronoi_name = "MAVEN_STATIC_H+_VDF_MVA_voronoi_" *
    Dates.format(ts, "yyyymmddTHHMM") * "-" *
    Dates.format(te, "HHMM") * ".pdf"

fig3 = Figure(resolution = (1200, 550), fontsize = 15)

# ---- xz 平面 Voronoi: B̂ vs B̂×V̂sw ----
# 无过滤, 使用全部数据点
vx_xz = V_rot_xz[:, 1]  # v 在 B̂_orth 方向的分量
vz_xz = V_rot_xz[:, 3]  # v 在 B̂×V̂sw 方向的分量
eflux_xz = eflux_all

ax_vxz = Axis(fig3[1, 1],
    xlabel = L"$v \cdot \hat{B}$ (km/s)",
    ylabel = L"$v \cdot (\hat{B} \times (\hat{B} \times \hat{V}_{\mathrm{sw}}))$ (km/s)",
    title = "xz plane (B̂ vs B̂×(B̂×V̂sw))",
    xticks = [-1000, -500, 0, 500, 1000], yticks = [-1000, -500, 0, 500, 1000])
xlims!(ax_vxz, plot_range[1], plot_range[2])
ylims!(ax_vxz, plot_range[1], plot_range[2])
hlines!(ax_vxz, 0, color = :white, linestyle = :dash)
vlines!(ax_vxz, 0, color = :white, linestyle = :dash)

if length(vx_xz) > 3
    pts_xz = [Point2f(vx_xz[i], vz_xz[i]) for i in eachindex(vx_xz)]
    tri_xz = triangulate(pts_xz)
    vor_xz = voronoi(tri_xz)
    vor_colors_xz = color_map_voronoi(eflux_xz, voronoi_colorrange)
    voronoiplot!(ax_vxz, vor_xz, color = vor_colors_xz, colormap = :viridis,
        strokewidth = 0, markersize = 0)
end

# 画激波面方向 (黑色虚线)
if norm(shock_xz_dir) > 1e-10
    L = Float64(plot_range[2]) * 0.9
    lines!(ax_vxz, [-L*shock_xz_dir[1], L*shock_xz_dir[1]],
                [-L*shock_xz_dir[2], L*shock_xz_dir[2]],
                color=:black, linestyle=:dash, linewidth=2)
end

# ---- yz 平面 Voronoi: V̂sw vs B̂×V̂sw ----
# 无过滤, 使用全部数据点
vy_yz = V_rot_yz[:, 1]  # v 在 V̂sw 方向的分量
vz_yz = V_rot_yz[:, 3]  # v 在 B̂×V̂sw 方向的分量
eflux_yz = eflux_all

ax_vyz = Axis(fig3[1, 2],
    xlabel = L"$v \cdot (\hat{B} \times \hat{V}_{\mathrm{sw}})$ (km/s)",
    ylabel = L"$v \cdot (\hat{B} \times (\hat{B} \times \hat{V}_{\mathrm{sw}}))$ (km/s)",
    title = "yz plane (B̂×V̂sw vs B̂×(B̂×V̂sw))",
    xticks = [-1000, -500, 0, 500, 1000], yticks = [-1000, -500, 0, 500, 1000])
xlims!(ax_vyz, plot_range[1], plot_range[2])
ylims!(ax_vyz, plot_range[1], plot_range[2])
hlines!(ax_vyz, 0, color = :white, linestyle = :dash)
vlines!(ax_vyz, 0, color = :white, linestyle = :dash)

if length(vy_yz) > 3
    pts_yz = [Point2f(vy_yz[i], vz_yz[i]) for i in eachindex(vy_yz)]
    tri_yz = triangulate(pts_yz)
    vor_yz = voronoi(tri_yz)
    vor_colors_yz = color_map_voronoi(eflux_yz, voronoi_colorrange)
    voronoiplot!(ax_vyz, vor_yz, color = vor_colors_yz, colormap = :viridis,
        strokewidth = 0, markersize = 0)
end

# 画激波面方向 (黑色虚线)
if norm(shock_yz_dir) > 1e-10
    L = Float64(plot_range[2]) * 0.9
    lines!(ax_vyz, [-L*shock_yz_dir[1], L*shock_yz_dir[1]],
                [-L*shock_yz_dir[2], L*shock_yz_dir[2]],
                color=:black, linestyle=:dash, linewidth=2)
end

Colorbar(fig3[1, 3], limits = voronoi_colorrange, colormap = :viridis, scale = log10,
    label = L"H$^{+}$ Differential Energy Flux (eV/cm²/s/sr/eV)")

Label(fig3[0, :], "MAVEN STATIC H⁺ VDF $(datestr) $(Dates.format(ts,"HH:MM:SS"))-$(Dates.format(te,"HH:MM:SS")) (MVA) Voronoi",
    fontsize = 18, tellwidth = false)

save(joinpath(VDF_path, fig_voronoi_name), fig3)
println("Saved Voronoi: ", joinpath(VDF_path, fig_voronoi_name))
