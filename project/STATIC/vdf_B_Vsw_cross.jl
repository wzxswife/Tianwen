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
ts = DateTime(date) + Dates.Hour(7) + Dates.Minute(25)
te = ts + Dates.Minute(1)

# --- 太阳风速度估算参数 (eflux/E² 峰值搜寻 bin) ---
sw_i = 2
sw_j = 14

# --- 过滤参数 ---
flux_threshold = 10^7.5   # eflux 阈值 (eV/cm²/s/sr/eV)
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
println("太阳风速度 (STATIC系): ", round.(MVN_staticvswrough, digits = 1), " km/s, |Vsw| = ", round(norm(MVN_staticvswrough), digits = 1), " km/s")

# ============================================================
# 构建坐标系: x̂=B̂, ŷ=V̂sw, ẑ=B̂×V̂sw (STATIC 系)
# ============================================================
e_x = MVN_parastatic
e_y = MVN_staticvswrough ./ norm(MVN_staticvswrough)
e_z = cross(e_x, e_y)
e_z = e_z ./ norm(e_z)

println("\n基矢量 (STATIC系):")
println("  x̂ = B̂      = ($(round.(e_x, digits = 4)))")
println("  ŷ = V̂sw    = ($(round.(e_y, digits = 4)))")
println("  ẑ = B̂×V̂sw  = ($(round.(e_z, digits = 4)))")
println("  x̂·ŷ = $(round(dot(e_x, e_y), digits = 3)) (非正交)")

# 基矢量转换到 MSO 坐标系 (用于热力图坐标转换)
e_x_mso = MVN_static2msomean * e_x  # B̂ in MSO
e_y_mso = MVN_static2msomean * e_y  # V̂sw in MSO
e_z_mso = MVN_static2msomean * e_z  # B̂×V̂sw in MSO

# ============================================================
# 逐帧遍历: 每个时刻 × 每个 bin → 速度矢量 + eflux
# ============================================================
n_frames = length(MVN_staticind)
n_total = n_frames * 64 * 32
V_all = Array{Float64}(undef, n_total, 3)
eflux_all = Array{Float64}(undef, n_total)
energy_all = Array{Float64}(undef, n_total)  # 每个 bin 的离子能量 (eV)

println("\n处理 $n_frames 帧 × $(64 * 32) bin = $n_total 点...")
for (fi, t_idx) in enumerate(MVN_staticind)
    swp = MVN_staticswp[t_idx]
    base_row = (fi - 1) * 64 * 32
    for (bin_idx, (ij, k)) in enumerate(Iterators.product(1:64, 1:32))
        row = base_row + bin_idx
        dir = CartesianFromSpherical()(
            CoordinateTransformations.Spherical(1.0,
                deg2rad(MVN_staticphi[1, ij, k, swp, 1]),
                deg2rad(MVN_statictheta[1, ij, k, swp, 1])))
        V_all[row, :] = MVN_staticvelocity[1, ij, k, swp, 1] .* Array(dir)
        eflux_all[row] = MVN_staticeflux[t_idx, 1, ij, k]
        energy_all[row] = MVN_staticenergy[1, ij, k, swp, 1]
    end
    fi % 10 == 0 && println("  帧 $fi / $n_frames")
end

# ============================================================
# 投影到新基: (v·x̂, v·ŷ, v·ẑ)
# ============================================================
proj = V_all * hcat(e_x, e_y, e_z) # (n_total, 3)

# 能量排除掩码: true 表示保留
energy_keep = trues(n_total)
if energy_exclude_min !== nothing && energy_exclude_max !== nothing
    energy_keep .&= .!(energy_exclude_min .<= energy_all .<= energy_exclude_max)
    println("排除能量范围: $energy_exclude_min ~ $energy_exclude_max eV, 排除点数: $(count(.!energy_keep))")
elseif energy_exclude_min !== nothing
    energy_keep .&= .!(energy_all .>= energy_exclude_min)
    println("排除能量 >= $energy_exclude_min eV, 排除点数: $(count(.!energy_keep))")
elseif energy_exclude_max !== nothing
    energy_keep .&= .!(energy_all .<= energy_exclude_max)
    println("排除能量 <= $energy_exclude_max eV, 排除点数: $(count(.!energy_keep))")
end

keep = (eflux_all .> flux_threshold) .&
       energy_keep .&
       (proj[:, 1] .> plot_range[1]) .& (proj[:, 1] .< plot_range[2]) .&
       (proj[:, 3] .> plot_range[1]) .& (proj[:, 3] .< plot_range[2])
println("有效点: $(count(keep)) / $n_total")

# ============================================================
# 绘图 (CairoMakie 散点)
# ============================================================
CairoMakie.activate!()

fig_name = "MAVEN_STATIC_H+_VDF_BVsw_xz_yz_frame_" *
    Dates.format(ts, "yyyymmddTHHMM") * "-" *
    Dates.format(te, "HHMM") * ".pdf"

crange = (flux_threshold, 1e9)

fig = Figure(resolution = (1200, 550), fontsize = 15)

# ---- xz 平面: v·B̂ vs v·B̂×V̂sw ----
ax1 = Axis(fig[1, 1],
    xlabel = L"$v \cdot \hat{B}$ (km/s)",
    ylabel = L"$v \cdot (\hat{B} \times \hat{V}_{\mathrm{sw}})$ (km/s)",
    title = "xz plane (B̂ vs B̂×V̂sw)",
    xticks = [-1000, -500, 0, 500, 1000], yticks = [-1000, -500, 0, 500, 1000])
xlims!(ax1, -1100, 1100)
ylims!(ax1, -1100, 1100)
hlines!(ax1, 0, color = :white, linestyle = :dash)
vlines!(ax1, 0, color = :white, linestyle = :dash)
scatter!(ax1, proj[keep, 1], proj[keep, 3];
    color = eflux_all[keep], colormap = :gist_earth, colorscale = log10,
    colorrange = crange, markersize = 8)

# ---- yz 平面: v·V̂sw vs v·B̂×V̂sw ----
ax2 = Axis(fig[1, 2],
    xlabel = L"$v \cdot \hat{V}_{\mathrm{sw}}$ (km/s)",
    ylabel = L"$v \cdot (\hat{B} \times \hat{V}_{\mathrm{sw}})$ (km/s)",
    title = "yz plane (V̂sw vs B̂×V̂sw)",
    xticks = [-1000, -500, 0, 500, 1000], yticks = [-1000, -500, 0, 500, 1000])
xlims!(ax2, -1100, 1100)
ylims!(ax2, -1100, 1100)
hlines!(ax2, 0, color = :white, linestyle = :dash)
vlines!(ax2, 0, color = :white, linestyle = :dash)
scatter!(ax2, proj[keep, 2], proj[keep, 3];
    color = eflux_all[keep], colormap = :gist_earth, colorscale = log10,
    colorrange = crange, markersize = 8)

Colorbar(fig[1, 3], limits = crange, colormap = :gist_earth, scale = log10,
    label = L"H$^{+}$ Differential Energy Flux (eV/cm²/s/sr/eV)")

save(joinpath(VDF_path, fig_name), fig)
println("Saved figure: ", joinpath(VDF_path, fig_name))

# ============================================================
# 热力图: 球面样条插值 (参考 sta_VDF_eflux_interp.jl)
# ============================================================

# 热力图数据: 仅过滤能量排除 + 速度范围, 不过滤 eflux
keep_heat = energy_keep .&
    (proj[:, 1] .> plot_range[1]) .& (proj[:, 1] .< plot_range[2]) .&
    (proj[:, 3] .> plot_range[1]) .& (proj[:, 3] .< plot_range[2])
println("\n热力图有效点: $(count(keep_heat)) / $n_total")

# ============================================================
# 步骤1: 时间平均 (参考 sta_VDF_eflux_interp.jl)
# 逐帧构建球面样条 → 插值到公共方向 → 累加平均
# ============================================================
println("沿时间平均数据...")

n_energy_ch = size(MVN_staticeflux, 4)  # 32 能道
n_dir = size(MVN_staticeflux, 3)        # 64 方向

# 使用窗口首帧的 swp 和方向作为公共方向
swp_ref = MVN_staticswp[MVN_staticind[1]]
theta_ref = MVN_statictheta[1, :, :, swp_ref, 1]  # 64×32
phi_ref = MVN_staticphi[1, :, :, swp_ref, 1]      # 64×32
energy_ref = MVN_staticenergy[1, :, :, swp_ref, 1] # 64×32
velocity_ref = MVN_staticvelocity[1, :, :, swp_ref, 1] # 64×32

# 公共方向: SC → (B̂, V̂sw, B̂×V̂sw)
dirs_common = zeros(n_dir * n_energy_ch, 3)
speeds_common = zeros(n_dir * n_energy_ch)
global row_common = 0
for (ij, k) in Iterators.product(1:n_dir, 1:n_energy_ch)
    global row_common += 1
    dir_sc = CartesianFromSpherical()(
        CoordinateTransformations.Spherical(1.0,
            deg2rad(phi_ref[ij, k]),
            deg2rad(theta_ref[ij, k])))
    v_sc = velocity_ref[ij, k] .* Array(dir_sc)
    # SC → MSO → (B̂, V̂sw, B̂×V̂sw)
    Q_first = QuatRotation(MVN_staticquat_mso[MVN_staticind[1], :])
    v_mso = Q_first * v_sc
    dirs_common[row_common, :] = v_mso ./ norm(v_mso)
    speeds_common[row_common] = norm(v_sc)
end

# 逐能道: 公共方向 + 累加器
common = Vector{Dict{Symbol,Any}}(undef, n_energy_ch)
for k in 1:n_energy_ch
    mask_k = ((k-1)*n_dir+1):(k*n_dir)
    common[k] = Dict(
        :dirs_MSO => dirs_common[mask_k, :],
        :speeds => speeds_common[mask_k],
        :n_total => n_dir,
        :flux_sum => zeros(Float64, n_dir),
        :cnt => zeros(Int, n_dir),
    )
end

# 逐帧: 构建球面样条 → 插值到公共方向 → 累加
global n_frames_ok = 0
for t_idx in MVN_staticind
    swp = MVN_staticswp[t_idx]
    Q_f = QuatRotation(MVN_staticquat_mso[t_idx, :])

    for k in 1:n_energy_ch
        theta_k = MVN_statictheta[1, :, k, swp, 1]
        phi_k = MVN_staticphi[1, :, k, swp, 1]
        flux_k = MVN_staticeflux[t_idx, 1, :, k]

        # 筛选有效数据
        valid = .!isnan.(theta_k) .& .!isnan.(phi_k) .& .!isnan.(flux_k) .& (flux_k .> 0)
        nv = count(valid)
        nv ≥ 4 || continue

        ct = cosd.(theta_k[valid])
        st = sind.(theta_k[valid])
        cp = cosd.(phi_k[valid])
        sp = sind.(phi_k[valid])
        dirs_SC_this = hcat(ct.*cp, ct.*sp, st)

        # 球面样条: SC 方向 → log10(eflux)
        spline = InterpolationSpline(dirs_SC_this, log10.(flux_k[valid]))

        # 对每个公共方向, 插值该帧的 eflux
        c = common[k]
        for i in 1:c[:n_total]
            d_mso = c[:dirs_MSO][i, :]
            d_sc = normalize(Q_f' * d_mso)  # MSO → SC (该帧姿态)
            log_flux_val = spline(Vector(d_sc))
            if !isnan(log_flux_val)
                c[:flux_sum][i] += 10.0^log_flux_val  # 还原线性 eflux
                c[:cnt][i] += 1
            end
        end
    end

    global n_frames_ok += 1
end
println("时间平均: $n_frames_ok 帧")

# 步骤2: 平均 → 构建最终球面样条
println("构建平均后的球面样条...")

# 平均 eflux
eflux_avg = zeros(n_dir, n_energy_ch)
for k in 1:n_energy_ch
    c = common[k]
    for i in 1:c[:n_total]
        if c[:cnt][i] > 0
            eflux_avg[i, k] = c[:flux_sum][i] / c[:cnt][i]
        end
    end
end

# 速度通道 (公共方向的3D全速)
v_ch = zeros(n_energy_ch)
for k in 1:n_energy_ch
    v_ch[k] = mean(common[k][:speeds])
end

# 按速度排序 (低→高), 与 sta_VDF_eflux_interp.jl 一致
sort_idx = sortperm(v_ch)
v_ch = v_ch[sort_idx]
eflux_avg = eflux_avg[:, sort_idx]

# 逐能道构建球面样条 (公共 MSO 方向 → log10(eflux_avg))
splines = Vector{Any}(undef, n_energy_ch)
for k in 1:n_energy_ch
    splines[k] = nothing
    flux_k = eflux_avg[:, k]
    valid = .!isnan.(flux_k) .& (flux_k .> 0)
    nv = count(valid)
    nv ≥ 4 || continue

    dirs_k = common[sort_idx[k]][:dirs_MSO][valid, :]
    splines[k] = InterpolationSpline(dirs_k, log10.(flux_k[valid]))
end
n_spline_ok = count(s -> s !== nothing, splines)
println("构建 $n_spline_ok / $n_energy_ch 个球面样条")

# 在 (B̂, V̂sw, B̂×V̂sw) 坐标系下构建球面样条插值热力图
# 参考 sta_VDF_eflux_interp.jl 的 build_vdf_grid 方法
# plane: "xz" → (v·B̂, v·(B̂×V̂sw)), "yz" → (v·V̂sw, v·(B̂×V̂sw))
# 使用时间平均后的 splines 和 v_ch
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
            dir_mso[3] = v2   # v·(B̂×V̂sw)
        else  # "yz"
            dir_mso[2] = v1   # v·V̂sw
            dir_mso[3] = v2   # v·(B̂×V̂sw)
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

println("球面样条插值生成热力图...")
xz_centers, xz_grid = build_vdf_grid_bvsw(splines, v_ch; n_grid=heatmap_nbin, plot_range=plot_range, plane="xz",
    e_x_mso=e_x_mso, e_y_mso=e_y_mso, e_z_mso=e_z_mso)
yz_centers, yz_grid = build_vdf_grid_bvsw(splines, v_ch; n_grid=heatmap_nbin, plot_range=plot_range, plane="yz",
    e_x_mso=e_x_mso, e_y_mso=e_y_mso, e_z_mso=e_z_mso)

fig_heat_name = "MAVEN_STATIC_H+_VDF_BVsw_xz_yz_heatmap_" *
    Dates.format(ts, "yyyymmddTHHMM") * "-" *
    Dates.format(te, "HHMM") * ".pdf"

# 颜色范围: 从有效数据的 5%~95% 分位数计算
all_valid = filter(isfinite, vec(xz_grid))
append!(all_valid, filter(isfinite, vec(yz_grid)))
if length(all_valid) > 10
    lo, hi = quantile(all_valid, [0.05, 0.95])
    cr_log = (lo, hi)
else
    cr_log = (6.0, 8.0)  # fallback
end

ticks = [-1000, -500, 0, 500, 1000]
title_heat = "MAVEN STATIC H⁺ VDF $(datestr) $(Dates.format(ts,"HH:MM:SS"))-$(Dates.format(te,"HH:MM:SS")) (B̂, V̂sw)"
colorbar_label = L"$\log_{10}$ Differential Energy Flux (eV/cm$^2$/s/sr/eV)"

fig2 = Figure(resolution = (1200, 550), fontsize = 15)

ax_xz = Axis(fig2[1, 1], aspect = 1, yticks = ticks, xticks = ticks,
    limits = (plot_range, plot_range), alignmode = Inside())
ax_yz = Axis(fig2[1, 2], aspect = 1, yticks = ticks, xticks = ticks,
    limits = (plot_range, plot_range), alignmode = Inside())
linkxaxes!(ax_xz, ax_yz)

# ---- xz 平面热力图 ----
heatmap!(ax_xz, xz_centers, xz_centers, xz_grid; colormap = :viridis,
    colorrange = cr_log, nan_color = :transparent)
hlines!(ax_xz, 0, linestyle = :dash, color = :white)
vlines!(ax_xz, 0, linestyle = :dash, color = :white)
ax_xz.xlabel = L"$v \cdot \hat{B}$ (km/s)"
ax_xz.ylabel = L"$v \cdot (\hat{B} \times \hat{V}_{\mathrm{sw}})$ (km/s)"

# ---- yz 平面热力图 ----
heatmap!(ax_yz, yz_centers, yz_centers, yz_grid; colormap = :viridis,
    colorrange = cr_log, nan_color = :transparent)
hlines!(ax_yz, 0, linestyle = :dash, color = :white)
vlines!(ax_yz, 0, linestyle = :dash, color = :white)
ax_yz.xlabel = L"$v \cdot \hat{V}_{\mathrm{sw}}$ (km/s)"
ax_yz.ylabel = L"$v \cdot (\hat{B} \times \hat{V}_{\mathrm{sw}})$ (km/s)"

Colorbar(fig2[1, 3], limits = cr_log, colormap = :viridis,
    label = colorbar_label)
Label(fig2[0, :], title_heat, fontsize = 18, tellwidth = false)

save(joinpath(VDF_path, fig_heat_name), fig2)
println("Saved heatmap: ", joinpath(VDF_path, fig_heat_name))
