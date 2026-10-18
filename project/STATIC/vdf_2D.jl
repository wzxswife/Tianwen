using Glob 
using DelimitedFiles, DataFrames, EzXML
using Dates
using LinearAlgebra
using ForwardDiff
using Statistics
using Interpolations
using Rotations
using Unitful
using ThreadsX
using WENO4
using Polyester
using CoordinateTransformations
import PhysicalConstants.CODATA2018: m_p, m_e, e, ε_0, μ_0, c_0, G, k_B
const M_Mars = 0.64169e24u"kg"
const R_Mars = 3389.5u"km"  #m
const au = 149597870.7u"km"

include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_STATIC.jl")
include("../../src/scripts/SphericalSplines.jl")
using .MAVEN_STATIC
using .MAVEN_load
using .SphericalSplines

date = Date(2022, 8, 24)
datestr = Dates.format.(date, "yyyymmdd")
year_num = Dates.format(date, "yyyy")
month_num = Dates.format(date, "mm")

ts = DateTime(date) + Dates.Hour(7) + Dates.Minute(27)
te = ts + Dates.Minute(1)

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

# MVN MAG
datejul = datetime2julian(DateTime(date))
MVN_magut1hz = MVN_mag2c1hz[:JulUTtime].-datejul
MVN_BMSO1hz = MVN_mag2c1hz[:B]
ind = ThreadsX.findall(ts .<= MVN_mag2c1hz[:epoch] .<= te)
MVN_Bmeanmso = dropdims(mean(Matrix(MVN_BMSO1hz[ind, :]), dims=(1,)), dims=(1,))

# MVN STATIC
MVN_staticswp = MVN_static[:swp_ind] |> Array 
MVN_staticswp .+= 1
MVN_statictime = MVN_static[:epoch] |> Array
MVN_staticjulian=datetime2julian.(DateTime.(MVN_statictime))
MVN_staticeflux = MVN_static[:eflux] |> Array
MVN_staticmass = MVN_static[:mass_arr] |> Array
MVN_staticenergy = MVN_static[:energy] |> Array
MVN_statictheta = MVN_static[:theta] |> Array
MVN_staticphi = MVN_static[:phi] |> Array
MVN_staticdenergy = MVN_static[:denergy] |> Array
MVN_staticdtheta = MVN_static[:dtheta] |> Array
MVN_staticdphi = MVN_static[:dphi] |> Array
MVN_staticdomega = MVN_static[:domega] |> Array
MVN_staticvelocity = sqrt.(MVN_staticenergy.*(2.0/ustrip(u"eV", m_p*c_0^2))./MVN_staticmass).*ustrip(u"km/s", c_0)
MVN_staticpos_sc_mso = MVN_static[:pos_sc_mso] |> Array
MVN_staticpos_sc_mso /= ustrip(u"km", R_Mars)

MVN_staticquat_mso = MVN_static[:quat_mso] |> Array
MVN_staticquat_sc = MVN_static[:quat_sc] |> Array
MVN_staticmsotrans = [QuatRotation(quat) for quat in eachcol(MVN_staticquat_mso')]
MVN_staticind = ThreadsX.findall(ts .<= MVN_statictime .<= te)
MVN_staticswpind = MVN_staticswp[MVN_staticind[1]]
MVN_staticefluxmean = dropdims(mean(MVN_staticeflux[MVN_staticind, :, :, :], dims = (1, ) ), dims=(1, ))
MVN_staticpsdmean = ustrip(u"cm^-3 * km^-3 * s^3", 0.5 * 1u"cm^-2 * s^-1" * (m_p/1u"eV")^2)  .* MVN_staticefluxmean .* (MVN_staticmass[:,:,:,MVN_staticswpind,1]./ MVN_staticenergy[:,:,:,MVN_staticswpind,1]).^2

MVN_staticindmean = round(Int64, mean(MVN_staticind))
MVN_statictimemean = MVN_statictime[MVN_staticindmean]
MVN_staticjulianmean = MVN_staticjulian[MVN_staticindmean]
MVN_static2msomean = QuatRotation(MVN_staticquat_mso[MVN_staticindmean, :])
MVN_Bmeanstatic = MVN_static2msomean' * MVN_Bmeanmso
MVN_parastatic = MVN_Bmeanstatic ./ norm(MVN_Bmeanstatic)
swi, swj = 2, 12
ind = argmax(MVN_staticefluxmean[1, (swj-1)*4+swi, :] ./ MVN_staticenergy[1, 1, :, MVN_staticswpind, 1].^2)
MVN_staticvswrough = CartesianFromSpherical()(CoordinateTransformations.Spherical(
    MVN_staticvelocity[1, (swj-1)*4+swi, ind, MVN_staticswpind, 1], 
    deg2rad(MVN_staticphi[1, (swj-1)*4+swi, ind, MVN_staticswpind, 1]), 
    deg2rad(MVN_statictheta[1, (swj-1)*4+swi, ind, MVN_staticswpind, 1])))
MVN_staticdirections = [CartesianFromSpherical()(CoordinateTransformations.Spherical(1.0, deg2rad(MVN_staticphi[1, ij, ind, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, ij, ind, MVN_staticswpind, 1]))) for ij in 1:64]
MVN_staticss = []
Snan = SphericalSplines.InterpolationSpline(stack(MVN_staticdirections; dims=1), fill(-30.0, 64))
for k in 1:32
    ind = findall(MVN_staticefluxmean[1, :, k] .> 0.0)
    if !isempty(ind)
        S = SphericalSplines.InterpolationSpline(stack(MVN_staticdirections[ind]; dims=1), log10.(MVN_staticefluxmean[1, ind, k]))
        push!(MVN_staticss, S)
    else
        push!(MVN_staticss, Snan)
    end  
end
MVN_perpstatic1 = cross(MVN_parastatic, MVN_staticvswrough)
MVN_perpstatic1 = MVN_perpstatic1 ./ norm(MVN_perpstatic1)
MVN_perpstatic2 = cross(MVN_perpstatic1, MVN_parastatic)
MVN_static2mfa = transpose(hcat(MVN_perpstatic1, MVN_parastatic, MVN_perpstatic2))

# 预计算有效能道（排除Snan）
valid_channels = [k for k in 1:32 if MVN_staticss[k] !== Snan]
valid_logv = log10.(MVN_staticvelocity[1, 1, valid_channels, MVN_staticswpind, 1])
println("Valid WENO4 channels: $(length(valid_channels))/27")

# ═══════════════════════════════════════════════════════
#  插值函数
# ═══════════════════════════════════════════════════════
function MVN_STATIC_eflux_interp(v3)
    vs = SphericalFromCartesian()(v3) 
    vtot = vs.r
    lon = vs.θ
    lat = vs.ϕ
    v_lo = MVN_staticvelocity[1, 1, valid_channels[1], MVN_staticswpind, 1]
    v_hi = MVN_staticvelocity[1, 1, valid_channels[end], MVN_staticswpind, 1]
    if (abs(rad2deg(lat)) <= 45.0) & (v_lo ≥ vtot ≥ v_hi)
        dir = v3 ./ vtot
        nvalid = length(valid_channels)
        logeflux = Vector{Float64}(undef, nvalid)
        for (idx, k) in enumerate(valid_channels)
            S = MVN_staticss[k]
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
#  3D直角网格插值 (与vdf_3D.jl相同)
# ═══════════════════════════════════════════════════════
vpara = collect(range(-1000.0, 1000.0, step = 10.0))
vperp1 = collect(range(-1000.0, 1000.0, step = 10.0))
vperp2 = collect(range(-1000.0, 1000.0, step = 10.0))
MVN_static_eflux_rect = Array{Float64, 3}(undef, length(vperp1), length(vpara), length(vperp2))
@batch for k in 1:length(vperp2)
    for j in 1:length(vpara)
        for i in 1:length(vperp1)
            v3 = Array(vpara[j] .* MVN_parastatic .+ vperp1[i] .* MVN_perpstatic1 .+ vperp2[k] .* MVN_perpstatic2)
            MVN_static_eflux_rect[i, j, k] = MVN_STATIC_eflux_interp(v3)
        end
    end
end
println("3D rectilinear grid interpolation completed.")

function nanmean(arr::AbstractArray{T,3}, dim::Int) where T
    # Count non-NaN along dim
    count_valid = sum(.!isnan.(arr), dims=dim)
    # Sum ignoring NaN
    sum_valid = sum(ifelse.(isnan.(arr), zero(T), arr), dims=dim)
    # Divide, keeping NaN where all were NaN
    result = sum_valid ./ count_valid
    result[count_valid .== 0] .= NaN
    return dropdims(result, dims=dim)
end

# Usage: average along dimension 3
result = nanmean(MVN_static_eflux_rect, 1)

using CairoMakie
CairoMakie.activate!()

fig = Figure(size = (800, 700))
ax = Axis(fig[1, 1], 
    xlabel = "v_perp2 (km/s)", 
    ylabel = "v_para (km/s)",
    aspect = 1)  # square cells

hm = heatmap!(ax, vperp2, vpara, result,
    colormap = :gist_earth,    # same as codebase
    colorscale = log10,        # log scale (common for flux)
    colorrange = (1e4, 1e8),   # adjust to your data range
    nan_color = :transparent   # NaN shows as transparent
)

Colorbar(fig[1, 2], hm, label = "Flux")
println("2D VDF heatmap generated.")

# ═══════════════════════════════════════════════════════
#  叠加标注: 太阳风方向 / 激波面投影 / 预期反射离子方向
#
#  本图所在平面 (MFA / 磁场坐标系):
#    横轴 x = vperp2 → 沿 MVN_perpstatic2 方向
#    纵轴 y = vpara  → 沿 MVN_parastatic (背景磁场 B) 方向
#    被积分掉的维度 = vperp1 (MVN_perpstatic1)
#
#  注意: MVN_perpstatic1 = normalize(cross(B, vsw))，即 perp1 ⊥ vsw。
#  因此太阳风/束流在 perp1 上的分量恒为 0，只有在 (perp2, para)
#  平面内才能完整看到太阳风方向。故本图取 (perp2, para) 平面。
#
#  任意 3D 向量 v (STATIC 坐标系) 在该平面上的投影分量:
#    (x, y) = (dot(v, MVN_perpstatic2), dot(v, MVN_parastatic))
# ═══════════════════════════════════════════════════════

# 把 STATIC 坐标系下的 3D 向量投影到 2D 平面坐标 (x=vperp2, y=vpara)
proj2d(v) = (dot(v, MVN_perpstatic2), dot(v, MVN_parastatic))

# ── 1. 激波法向 (MVA 最小方差分析) ──
# 复用 vdf_3D_slice.jl 的做法: 在激波穿越时刻附近取 MAG 数据做 MVA
# 特征值最小的方向即激波法向
function MVA(b::AbstractMatrix)
    nb = size(b, 1)
    bm = mean(b, dims = 1)
    muv = b' * b ./ nb - bm' * bm
    return eigen(muv)
end

# 激波穿越时刻 (按需修改). 默认取平均时刻
shock_time = DateTime(2022, 8, 24, 7, 27, 20)
shock_window = Second(30)
ind_shock = ThreadsX.findall(shock_time - shock_window .<= MVN_mag2c1hz[:epoch] .<= shock_time + shock_window)

n_shock_mso = zeros(3)
if length(ind_shock) >= 3
    B_shock = Matrix(MVN_BMSO1hz[ind_shock, :])
    eig_shock = MVA(B_shock)
    n_shock_mso = eig_shock.vectors[:, 1]         # 最小方差方向 = 法向
    if n_shock_mso[1] < 0
        n_shock_mso = -n_shock_mso                # 统一取指向太阳侧
    end
    println("Shock normal (MSO) = ", round.(n_shock_mso, digits = 3))
    println("MVA eigenvalues    = ", round.(eig_shock.values, digits = 2))
else
    @warn "MAG 数据在激波时刻附近不足, 跳过激波面/反射方向标注" n_mag = length(ind_shock)
end

# 法向由 MSO → STATIC. 用平均时刻的 STATIC→MSO 旋转求逆
# 之后所有矢量都统一在 STATIC 坐标系下, 与 proj2d 的输入保持一致
n_shock_static = MVN_static2msomean' * n_shock_mso
n_shock_static = n_shock_static ./ norm(n_shock_static)
n_shock_mfa = MVN_static2mfa * n_shock_static

# ── 2. 太阳风方向 (STATIC 坐标) ──
# MVN_staticvswrough 本身就是 STATIC 坐标系下、由 eflux 峰值估计的太阳风速度矢量
vsw_static = MVN_staticvswrough
vsw_mfa = MVN_static2mfa * vsw_static
vsw_hat_static = vsw_static ./ norm(vsw_static)
vsw_hat_mfa = vsw_mfa ./ norm(vsw_mfa)

# ── 3. 预期反射离子方向: v_ref = v - 2 (v·n̂) n̂ (均为 STATIC 坐标) ──
v_ref_static = vsw_static .- 2.0 .* dot(vsw_static, n_shock_static) .* n_shock_static
v_ref_mfa = vsw_hat_mfa .- 2.0 .* dot(vsw_mfa, n_shock_mfa) .* n_shock_mfa
v_ref_mfa = v_ref_mfa ./ norm(v_ref_mfa)
println("Solar wind (STATIC) = ", round.(vsw_static, digits = 1))
println("Shock normal(STATIC) = ", round.(n_shock_static, digits = 3))
println("Reflected   (STATIC) = ", round.(v_ref_static, digits = 1))
println("  vsw  in-plane (perp2, para) = ", round.(collect(proj2d(vsw_hat_static)), digits = 3))
println("  vref in-plane (perp2, para) = ", round.(collect(proj2d(v_ref_static ./ norm(v_ref_static))), digits = 3))

# ── 4. 绘图 ──
arrow_len = 900.0   # 标注箭头长度 (km/s), 可调整

# 太阳风方向 (红色箭头)
vsw_x, vsw_y = vsw_hat_mfa[2], vsw_hat_mfa[3]
arrows2d!(ax, [0.0], [0.0], [vsw_x * arrow_len], [vsw_y * arrow_len];
    color = :red, shaftwidth = 6, tipwidth = 20, tiplength = 20)
text!(ax, vsw_x * arrow_len, vsw_y * arrow_len; text = " Vsw",
    color = :red, fontsize = 16, align = (:left, :center))
    
# 预期反射离子方向 (橙色箭头)
# ref_x, ref_y = proj2d(v_ref_static ./ norm(v_ref_static))
ref_x, ref_y = v_ref_mfa[2], v_ref_mfa[3]
arrows2d!(ax, [0.0], [0.0], [ref_x * arrow_len], [ref_y * arrow_len];
    color = :orange, shaftwidth = 6, tipwidth = 20, tiplength = 20)
text!(ax, ref_x * arrow_len, ref_y * arrow_len; text = " Vref",
    color = :orange, fontsize = 16, align = (:left, :center))

# 激波面在 (vperp2, vpara) 平面上的投影 = 一条过原点的直线
# 平面: dot(v, n̂) = 0, v = x·e_perp2 + y·e_para  ⇒  x·n_perp2 + y·n_para = 0
# 取平面内的方向 (dx, dy) = (-n_para, n_perp2), 沿该方向画线
# np2 = dot(n_shock_static, MVN_perpstatic2)
# npa = dot(n_shock_static, MVN_parastatic)
np2, npa = n_shock_mfa[2], n_shock_mfa[3]
shock_dir = (npa, -np2)                   # 垂直于 (np2, npa) 的平面内方向
shock_dir = shock_dir ./ norm(shock_dir)
shock_len = 1000.0                        # 直线半长 (km/s)
lines!(ax,
    [-shock_dir[1], shock_dir[1]] .* shock_len,
    [-shock_dir[2], shock_dir[2]] .* shock_len;
    color = :white, linewidth = 2, linestyle = :dash)
text!(ax, shock_dir[1] * shock_len, shock_dir[2] * shock_len; text = " Shock front",
    color = :white, fontsize = 16, align = (:left, :center))

# 图例 (仅说明标注含义, heatmap 走单独的 colorbar)
legend_elems = [
    LineElement(color = :red, linewidth = 3, linestyle = :solid),
    LineElement(color = :orange, linewidth = 3, linestyle = :solid),
    LineElement(color = :white, linewidth = 2, linestyle = :dash),
]
Legend(fig[1, 2, Bottom()],
    legend_elems, ["Solar wind Vsw", "Reflected Vref", "Shock front"];
    framevisible = false, orientation = :horizontal, fontsize = 14)


vdf_2d_name = "MVN_STATIC_eflux_2D_$(Dates.format(ts, "yyyymmdd_HHMM")).png"
save(joinpath(VDF_path, vdf_2d_name), fig)
println("2D VDF heatmap saved to: ", joinpath(VDF_path, vdf_2d_name))