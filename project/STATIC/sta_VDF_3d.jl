using Glob
using Dates
using Statistics, LinearAlgebra
using Rotations
using CoordinateTransformations
using WENO4
include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/SphericalSplines.jl")
using .MAVEN_load
using .SphericalSplines

# ═══════════════════════════════════════════════════════
#  参数设置
# ═══════════════════════════════════════════════════════
date = Date(2022, 8, 24)
datestr = Dates.format.(date, "yyyymmdd")
year_num = Dates.format(date, "yyyy")
month_num = Dates.format(date, "mm")

ts = DateTime(date) + Dates.Hour(7) + Dates.Minute(27)
te = ts + Dates.Minute(3)

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

# ═══════════════════════════════════════════════════════
#  读取数据 (ScienceReadData.jl 方式)
# ═══════════════════════════════════════════════════════
mag_data = MAVEN_load.load_mag_l2(MAG_file[1])
sta_data = MAVEN_load.load_cdf(STA_file[1])
println("Data loaded.")

# load_cdf 返回形状:
#   eflux:     (ntime, 8, 64, 32) — record-varying, mass/bin/energy
#   energy:    (8, 64, 32, nswp)  — mass/bin/energy/swp
#   phi/theta: (8, 64, 32, nswp)
#   mass_arr:  (8, 64, 32, nswp)
#   swp_ind:   (ntime,)
#   quat_mso:  (ntime, 4)
#   epoch:     (ntime,)

MVN_staswp = sta_data[:swp_ind] .+ 1   # +1 for 1-based indexing
MVN_statime = sta_data[:epoch]
MVN_staeflux = sta_data[:eflux]         # (ntime, 8, 64, 32)
MVN_staenergy = sta_data[:energy]       # (8, 64, 32, nswp)
MVN_statheta = sta_data[:theta]         # (8, 64, 32, nswp)
MVN_staphi = sta_data[:phi]             # (8, 64, 32, nswp)
MVN_stamass = sta_data[:mass_arr]       # (8, 64, 32, nswp)

# 速度计算: v = c * sqrt(2E / (m_p*c^2 * mass_amu))
MVN_stavelocity = sqrt.(MVN_staenergy .* (2.0 / 938.272e6) ./ MVN_stamass) .* 299792.458  # km/s, (8,64,32,nswp)

# 四元数 (SC→MSO)
MVN_staat_mso = sta_data[:quat_mso]     # (ntime, 4)
MVN_sta_msotrans = [QuatRotation(MVN_staat_mso[i, :]) for i in 1:size(MVN_staat_mso, 1)]

# ═══════════════════════════════════════════════════════
#  时间筛选与平均
# ═══════════════════════════════════════════════════════
MVN_swpind = MVN_staswp[1]
time_inds = findall(ts .<= MVN_statime .<= te)
if isempty(time_inds)
    error("No STATIC data found in [ts, te]")
end
@info "Found $(length(time_inds)) STATIC time frames"

# eflux 时间平均: (ntime, 8, 64, 32) → mean over time → (8, 64, 32)
MVN_staefluxmean = dropdims(mean(MVN_staeflux[time_inds, :, :, :], dims=1), dims=1)
@info "eflux mean shape: $(size(MVN_staefluxmean))"
@info "eflux mean range: $(extrema(MVN_staefluxmean[MVN_staefluxmean .> 0]))"

# 磁场时间平均 (MSO)
MVN_Bmeanmso = vec(mean(mag_data[:B][findall(ts .<= mag_data[:epoch] .<= te), :], dims=1))

# 平均时间索引 → 坐标变换
MVN_staindmean = round(Int64, mean(time_inds))
MVN_stamsomean = MVN_sta_msotrans[MVN_staindmean]

# 磁场 → STATIC 坐标系
MVN_Bmeanstatic = MVN_stamsomean' * MVN_Bmeanmso
MVN_parastatic = MVN_Bmeanstatic ./ norm(MVN_Bmeanstatic)

# ═══════════════════════════════════════════════════════
#  太阳风速度方向 (ScienceReadData.jl 方式)
# ═══════════════════════════════════════════════════════
# 找 eflux/energy^2 最大值对应的 bin → 太阳风速度方向
i_sw, j_sw = 2, 14
bin_sw = (j_sw - 1) * 4 + i_sw   # = 54
ind_max = argmax(vec(MVN_staefluxmean[1, bin_sw, :]) ./
    vec(MVN_staenergy[1, bin_sw, :, MVN_swpind]).^2)
MVN_stavswrough = Array(CartesianFromSpherical()(
    CoordinateTransformations.Spherical(
        MVN_stavelocity[1, bin_sw, ind_max, MVN_swpind],
        deg2rad(MVN_staphi[1, bin_sw, ind_max, MVN_swpind]),
        deg2rad(MVN_statheta[1, bin_sw, ind_max, MVN_swpind])
    )))
@info "Vsw direction from STATIC: $(round.(MVN_stavswrough, digits=3))"

# ═══════════════════════════════════════════════════════
#  构建 MFA 坐标系 (ScienceReadData.jl 方式)
# ═══════════════════════════════════════════════════════
# e1 (perp1) = B × Vsw, e2 (para) = B̂, e3 (perp2) = e1 × B̂
MVN_perpstatic1 = cross(MVN_parastatic, MVN_stavswrough)
MVN_perpstatic1 = MVN_perpstatic1 ./ norm(MVN_perpstatic1)
MVN_perpstatic2 = cross(MVN_perpstatic1, MVN_parastatic)
MVN_perpstatic2 = MVN_perpstatic2 ./ norm(MVN_perpstatic2)

# 变换矩阵: v_mfa = MVN_static2mfa * v_static
MVN_static2mfa = transpose(hcat(MVN_perpstatic1, MVN_parastatic, MVN_perpstatic2))

@info "MFA frame: B_hat = $(round.(MVN_parastatic, digits=3))"
@info "MFA frame: perp1 = $(round.(MVN_perpstatic1, digits=3))"
@info "MFA frame: perp2 = $(round.(MVN_perpstatic2, digits=3))"

# ═══════════════════════════════════════════════════════
#  球面样条 (ScienceReadData.jl 方式: phi/theta → 方向)
# ═══════════════════════════════════════════════════════
n_energy = 32
n_bins = 64

# 能量网格 → 速度 (km/s), 取 mass=1 (H+)
MVN_staenergygrid = vec(MVN_stavelocity[1, 1, :, MVN_swpind])
@info "Energy grid (km/s): $(round.(MVN_staenergygrid[1:5], digits=1)) ... $(round.(MVN_staenergygrid[end-4:end], digits=1))"

# 方向向量 (从 phi/theta 直接计算, 用参考能道)
ie_ref = 10
MVN_stadirections = [CartesianFromSpherical()(
    CoordinateTransformations.Spherical(
        1.0,
        deg2rad(MVN_staphi[1, ij, ie_ref, MVN_swpind]),
        deg2rad(MVN_statheta[1, ij, ie_ref, MVN_swpind])
    )) for ij in 1:n_bins]

# 按能道构建球面样条 (ScienceReadData.jl: 逐能道, mass=1)
Snan = SphericalSplines.InterpolationSpline(
    stack(MVN_stadirections; dims=1), fill(-30.0, n_bins))
MVN_stass = []
for k in 1:n_energy
    ind = findall(vec(MVN_staefluxmean[1, :, k]) .> 0.0)
    if !isempty(ind)
        S = SphericalSplines.InterpolationSpline(
            stack(MVN_stadirections[ind]; dims=1),
            log10.(vec(MVN_staefluxmean[1, ind, k])))
        push!(MVN_stass, S)
    else
        push!(MVN_stass, Snan)
    end
end
n_spline_ok = count(S -> S !== Snan, MVN_stass)
@info "Built $n_spline_ok/$n_energy energy channel splines"

# ═══════════════════════════════════════════════════════
#  WENO4 跨能道插值函数 (ScienceReadData.jl 方式)
#  只用有效能道 (非Snan) 做 WENO4 插值, 避免 -30 污染
# ═══════════════════════════════════════════════════════
# 预计算: 哪些能道(6~32)有有效样条
valid_channels = [k for k in 6:n_energy if MVN_stass[k] !== Snan]
valid_logv = log10.(MVN_staenergygrid[valid_channels])
@info "Valid WENO4 channels: $(length(valid_channels))/27 (indices $(valid_channels[1])..$(valid_channels[end]))"

function MVN_STATIC_eflux_interp(v3)
    vs = SphericalFromCartesian()(v3)
    vtot = vs.r
    lat = vs.ϕ
    v_lo = MVN_staenergygrid[valid_channels[1]]   # 最高速度 (最低能道)
    v_hi = MVN_staenergygrid[valid_channels[end]]  # 最低速度 (最高能道)
    if (abs(rad2deg(lat)) <= 45.0) & (v_lo >= vtot >= v_hi)
        dir = Vector{Float64}(v3 ./ vtot)
        nvalid = length(valid_channels)
        logeflux = Vector{Float64}(undef, nvalid)
        for (idx, k) in enumerate(valid_channels)
            S = MVN_stass[k]
            logeflux[idx] = S(dir)
        end
        x = reverse(valid_logv)
        y = reverse(logeflux)
        return 10.0^interpolate_weno4([log10(vtot)], x, y)[1]
    else
        return NaN64
    end
end

# ═══════════════════════════════════════════════════════
#  3D 速度网格插值
# ═══════════════════════════════════════════════════════
v_min = MVN_staenergygrid[32]
v_max = MVN_staenergygrid[1]
@info "Velocity range: $v_min - $v_max km/s"

n_grid = 80
v_range = range(-1000, 1000, length=n_grid)
v_para  = collect(v_range)
v_perp1 = collect(v_range)
v_perp2 = collect(v_range)

grid_3d = fill(NaN, n_grid, n_grid, n_grid)

@info "Interpolating on $(n_grid)×$(n_grid)×$(n_grid) grid..."

n_interp_ok = 0
for k in 1:n_grid, j in 1:n_grid, i in 1:n_grid
    v3 = v_para[j] .* MVN_parastatic .+ v_perp1[i] .* MVN_perpstatic1 .+ v_perp2[k] .* MVN_perpstatic2
    grid_3d[i, j, k] = MVN_STATIC_eflux_interp(v3)
    if isfinite(grid_3d[i, j, k])
        global n_interp_ok += 1
    end
end

@info "Interpolation done: $n_interp_ok / $(n_grid^3) grid points filled"

# 诊断: 打印 grid_3d 值分布
finite_vals = grid_3d[isfinite.(grid_3d)]
if !isempty(finite_vals)
    log_vals = log10.(finite_vals)
    @info "Grid finite values: n=$(length(finite_vals)), min=$(minimum(finite_vals)), max=$(maximum(finite_vals))"
    @info "Log10 range: $(minimum(log_vals)) - $(maximum(log_vals))"
    # 分位数
    for q in [0.1, 0.25, 0.5, 0.75, 0.9]
        @info "  log10 quantile $q: $(quantile(log_vals, q))"
    end
end

# 诊断: 每个 v_para 切片有多少有限值
for j in 1:n_grid
    slice = @view grid_3d[:, j, :]
    cnt = count(isfinite, slice)
    if cnt > 0
        sv = slice[isfinite.(slice)]
        @info "  slice ipara=$j: $cnt finite, log10=[$(round(log10(minimum(sv)),digits=1)), $(round(log10(maximum(sv)),digits=1))]"
    end
end

# ═══════════════════════════════════════════════════════
#  3D 绘图 (GLMakie) — 匹配 MINPA_read.jl 方式
#  传入线性值, colorrange 用线性范围, colorscale=log10
# ═══════════════════════════════════════════════════════
using GLMakie
GLMakie.activate!()

# 找到 v_perp1 方向数据最多的切片 (用于右图)
global best_iperp1 = 0
global best_cnt = 0
for i in 1:n_grid
    cnt = count(isfinite, @view grid_3d[i, :, :])
    if cnt > best_cnt
        global best_cnt = cnt
        global best_iperp1 = i
    end
end
@info "Right plot: best_iperp1=$best_iperp1 with $best_cnt finite points"

fig_name = "MAVEN_STATIC_H+_VDF_3D_" * Dates.format(ts, "yyyymmddTHHMMSS") *
    "-" * Dates.format(te, "HHMMSS") * ".png"

size_inches = (16, 9)
size_pt = 72 .* size_inches
fig = Figure(size = size_pt, fontsize = 15,
    figure_padding = (5, 10, 9, 8))

title_str = "MAVEN STATIC H⁺ VDF 3D " *
    Dates.format(ts, "yyyymmdd HH:MM:SS") * "-" *
    Dates.format(te, "HH:MM:SS")
Label(fig[1, 1, Top()], title_str, fontsize = 18, tellheight = true)

# 左图: 沿 v_para 切片 — 使用线性值 + colorscale=log10
ls1 = LScene(fig[1, 1])
colorrange = (minimum(finite_vals), maximum(finite_vals))
@info "Linear colorrange: $(colorrange[1]) - $(colorrange[2])"

# 先画第1个有数据的切片作为 Colorbar 引用
hm_ref = nothing
for ipara in 1:n_grid
    cnt = count(isfinite, @view grid_3d[:, ipara, :])
    if cnt >= 4
        hm_ref = heatmap!(ls1, v_perp1, v_perp2, grid_3d[:, ipara, :],
            transformation = (:xz, v_para[ipara]),
            colormap = :viridis, colorrange = colorrange, colorscale = log10)
        break
    end
end
# 再画其余有数据的切片
for ipara in 2:n_grid
    cnt = count(isfinite, @view grid_3d[:, ipara, :])
    if cnt >= 4
        heatmap!(ls1, v_perp1, v_perp2, grid_3d[:, ipara, :],
            transformation = (:xz, v_para[ipara]),
            colormap = :viridis, colorrange = colorrange, colorscale = log10)
    end
end

# 右图: 沿 v_perp1 切片
Label(fig[1, 2, Top()], title_str, fontsize = 18, tellheight = true)
ls2 = LScene(fig[1, 2])

heatmap!(ls2, v_para, v_perp2, grid_3d[best_iperp1, :, :],
    transformation = (:yz, v_perp1[best_iperp1]),
    colormap = :viridis, colorrange = colorrange, colorscale = log10)

# 磁场方向箭头 (MFA y轴)
ps = Point3f(0.0, 0.0, 0.0)
pe = Point3f(0.0, v_max * 0.8, 0.0)
arrows3d!(ls1, ps, pe, argmode = :endpoint,
    color = :red, shaftradius = 0.02, tipradius = 0.05, tailradius = 0.02)
arrows3d!(ls2, ps, pe, argmode = :endpoint,
    color = :red, shaftradius = 0.02, tipradius = 0.05, tailradius = 0.02)

if hm_ref !== nothing
    Colorbar(fig[1, 3], hm_ref, label = "Energy Flux (eV/cm²/s/sr/eV)")
else
    Colorbar(fig[1, 3], limits=colorrange, colormap=:viridis,
        label = "Energy Flux")
end

save(joinpath(VDF_path, fig_name), fig)
println("Saved figure: ", joinpath(VDF_path, fig_name))

# 持续显示窗口
display(fig)
while isopen(fig.scene)
    sleep(0.1)
end
