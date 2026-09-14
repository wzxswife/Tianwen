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

MVN_staticefluxmid = similar(MVN_staticefluxmean)
MVN_staticefluxuquant = similar(MVN_staticefluxmean)
MVN_staticefluxlquant = similar(MVN_staticefluxmean)
for i in axes(MVN_staticeflux, 2), j in axes(MVN_staticeflux, 3), k in axes(MVN_staticeflux, 4)
    data = MVN_staticeflux[MVN_staticind, i, j, k]
    MVN_staticefluxmid[i, j, k] = quantile(data, 0.5)
    MVN_staticefluxuquant[i, j, k] = quantile(data, 0.75)
    MVN_staticefluxlquant[i, j, k] = quantile(data, 0.25)
end
# i, j =2, 14
# MVN_staticenergyl = zeros(32)
# MVN_staticenergyu = zeros(32)
# for ie in 1:32
#     R = MVN_staticdenergy[1, (j-1)*4+i, ie, MVN_staticswpind, 1]./MVN_staticenergy[1, (j-1)*4+i, ie, MVN_staticswpind, 1]
#     ratio = 0.5R+sqrt(1.0+0.25R^2)
#     MVN_staticenergyl[ie] = MVN_staticenergy[1, (j-1)*4+i, ie, MVN_staticswpind, 1]/ratio
#     MVN_staticenergyu[ie] = MVN_staticenergy[1, (j-1)*4+i, ie, MVN_staticswpind, 1]*ratio
# end
# MVN_staticenergygridH = vcat(MVN_staticenergyu, MVN_staticenergyl[end])
# MVN_staticvelocitygridH = reverse(sqrt.(MVN_staticenergygridH.*(2.0/ustrip(u"eV", m_p*c_0^2))/1.03).*ustrip(u"km/s", c_0))
# MVN_staticphigridH = [MVN_staticphi[1, (j-1)*4+1, 7, MVN_staticswpind, 1]-MVN_staticdphi[1, (j-1)*4+1, 7, MVN_staticswpind, 1]*0.5 for j in 1:16]
# MVN_staticphigridH = vcat(MVN_staticphigridH, MVN_staticphigridH[begin]+360.0)
# MVN_staticthetau = 90.0.-[MVN_statictheta[1, i, 7, MVN_staticswpind, 1]-MVN_staticdtheta[1, i, 7, MVN_staticswpind, 1]*0.5 for i in 1:4]
# MVN_staticthetal = 90.0.-[MVN_statictheta[1, i, 7, MVN_staticswpind, 1]+MVN_staticdtheta[1, i, 7, MVN_staticswpind, 1]*0.5 for i in 1:4]
# MVN_staticthetagridH = vcat(MVN_staticthetal, MVN_staticthetau[end])
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
valid_channels = [k for k in 6:32 if MVN_staticss[k] !== Snan]
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
#  任意斜面截面插值函数
#  plane_origin: 平面上的一点 (3D向量, MFA坐标)
#  plane_normal: 平面法向量 (3D向量, 会被归一化)
#  range_s, range_t: 截面两个方向的范围 (km/s)
#  n_s, n_t: 两个方向的采样点数
#  返回: (grid_s, grid_t, eflux_2d)
#    grid_s, grid_t: 2D网格坐标 (用于绘图)
#    eflux_2d: 2D截面数据 (n_s × n_t)
# ═══════════════════════════════════════════════════════
function slice_plane_interp(plane_origin, plane_normal;
        range_s = (-1000.0, 1000.0), range_t = (-1000.0, 1000.0),
        n_s = 200, n_t = 200)

    # 归一化法向量
    n = plane_normal ./ norm(plane_normal)

    # 构建平面内正交基 (Gram-Schmidt)
    # 选一个不平行于n的辅助向量
    if abs(n[1]) < 0.9
        aux = [1.0, 0.0, 0.0]
    else
        aux = [0.0, 1.0, 0.0]
    end
    u = cross(aux, n)
    u = u ./ norm(u)
    v = cross(n, u)  # v = n × u, 保证右手系

    # 2D网格
    s_vals = collect(range(range_s[1], range_s[2], length = n_s))
    t_vals = collect(range(range_t[1], range_t[2], length = n_t))

    # 插值
    eflux_2d = Array{Float64, 2}(undef, n_s, n_t)
    @batch for it in 1:n_t
        for is in 1:n_s
            v3 = plane_origin .+ s_vals[is] .* u .+ t_vals[it] .* v
            eflux_2d[is, it] = MVN_STATIC_eflux_interp(v3)
        end
    end

    return s_vals, t_vals, eflux_2d, u, v
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

# ═══════════════════════════════════════════════════════
#  激波法向计算 (MVA最小方差分析)
# ═══════════════════════════════════════════════════════
function MVA(b::AbstractMatrix)
    nb = size(b, 1)
    bm = mean(b, dims=1)
    muv = b' * b ./ nb - bm' * bm
    return eigen(muv)
end

shock_time = DateTime(2022, 8, 24, 7, 27, 0)
shock_window = Second(10)
t_shock_start = shock_time - shock_window
t_shock_end = shock_time + shock_window

ind_shock = ThreadsX.findall(t_shock_start .<= MVN_mag2c1hz[:epoch] .<= t_shock_end)
B_shock = MVN_mag2c1hz[:B][ind_shock, :]

eig_shock = MVA(B_shock)
n_shock_mso = eig_shock.vectors[:, 1]
if n_shock_mso[1] < 0
    n_shock_mso = -n_shock_mso
end
n_shock_mfa = MVN_static2mfa * n_shock_mso

# 太阳风速度 (MFA坐标)
MVN_vsw_mfa = MVN_static2mfa * MVN_staticvswrough

# 反射方向: v_ref = v - 2(v·n)n
dot_vn = dot(MVN_vsw_mfa, n_shock_mfa)
v_reflected = MVN_vsw_mfa - 2.0 * dot_vn .* n_shock_mfa

println("Shock normal (MSO): $(round.(n_shock_mso, digits=3))")
println("Shock normal (MFA): $(round.(n_shock_mfa, digits=3))")
println("Eigenvalues: $(round.(eig_shock.values, digits=2))")
println("Solar wind (MFA): $(round.(MVN_vsw_mfa, digits=1))")
println("Reflected (MFA): $(round.(v_reflected, digits=1))")

# ═══════════════════════════════════════════════════════
#  截面定义: 太阳风与反射方向所在平面
# ═══════════════════════════════════════════════════════
# 截面法向 = 太阳风 × 反射方向 (垂直于入射-反射平面)
plane_origin = [0.0, 0.0, 0.0]
plane_normal = normalize(cross(MVN_vsw_mfa, v_reflected))

println("Computing oblique slice...")
println("  Plane normal (MFA): $(round.(plane_normal, digits=3))")
s_vals, t_vals, eflux_slice, u_vec, v_vec = slice_plane_interp(
    plane_origin, plane_normal;
    range_s = (-1000.0, 1000.0), range_t = (-1000.0, 1000.0),
    n_s = 200, n_t = 200)

println("  Slice shape: $(size(eflux_slice))")
finite_count = count(isfinite, eflux_slice)
println("  Finite values: $finite_count / $(length(eflux_slice))")

# ═══════════════════════════════════════════════════════
#  绘图: 左=3D体渲染 + 截面, 右=2D截面展开图
# ═══════════════════════════════════════════════════════
using GLMakie
GLMakie.activate!()

size_inches = (20, 9)
size_pt = 72 .* size_inches
fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 10, 9, 8))
title_str = "MAVEN STATIC H⁺ VDF " * Dates.format(ts, "yyyymmdd HH:MM") * "-" * Dates.format(te, "HH:MM")
fig_name = "MAVEN_STATIC_H+_VDF_slice_" * Dates.format(ts, "yyyymmddTHHMMSS") *
    "-" * Dates.format(te, "HHMMSS") * ".png"

# ── 左图: 3D场景 + 斜面截面 ──
Label(fig[1, 1, Top()], title_str * "\n3D + Oblique Slice", fontsize = 16, tellheight = true)
ls1 = LScene(fig[1, 1])

# 画几层 vpara 切片作为参考背景 (半透明)
for ipara in 1:2:length(vpara)
    heatmap!(ls1, vperp1, vperp2, MVN_static_eflux_rect[:, ipara, :],
        transformation = (:xz, vpara[ipara]),
        colormap = :gist_earth, colorrange = (1e4, 1e8), colorscale = log10)
end

# 在3D场景中画斜面截面
# 构建截面的4个角点 (在MFA坐标中)
corner_s = [-1000.0, 1000.0, 1000.0, -1000.0]
corner_t = [-1000.0, -1000.0, 1000.0, 1000.0]
# 转为3D坐标 (用于mesh!的vertices)
function plane_point(origin, u, v, s, t)
    return origin .+ s .* u .+ t .* v
end

# 用mesh画截面平面 (半透明底色)
# mesh! 需要Mesh primitive, 这里用 surface! 替代
# 实际用 lines! 画截面边框
border_s = vcat(range(-1000, 1000, length=50), fill(1000.0, 50), range(1000, -1000, length=50), fill(-1000.0, 50))
border_t = vcat(fill(-1000.0, 50), range(-1000, 1000, length=50))  # 简化: 画边框线

# 在3D场景中插入截面 (用 heatmap surface)
# 把截面数据转为 surface! 可用的格式
# surface! 需要 x, y, z 三个2D数组 + color 2D数组
n_s, n_t = size(eflux_slice)
xs = zeros(n_s, n_t)
ys = zeros(n_s, n_t)
zs = zeros(n_s, n_t)
for it in 1:n_t, is in 1:n_s
    pt = plane_origin .+ s_vals[is] .* u_vec .+ t_vals[it] .* v_vec
    xs[is, it] = pt[1]  # perp1
    ys[is, it] = pt[2]  # para
    zs[is, it] = pt[3]  # perp2
end

# 截面边框线 (在3D场景中可视化截面位置)
# 构建截面4条边的3D坐标
function make_plane_corners(origin, u, v, s_range, t_range)
    s_min, s_max = s_range
    t_min, t_max = t_range
    corners = [
        origin .+ s_min .* u .+ t_min .* v,
        origin .+ s_max .* u .+ t_min .* v,
        origin .+ s_max .* u .+ t_max .* v,
        origin .+ s_min .* u .+ t_max .* v,
        origin .+ s_min .* u .+ t_min .* v,  # 闭合
    ]
    return [Point3f(c...) for c in corners]
end

plane_corners = make_plane_corners(plane_origin, u_vec, v_vec, (-1000, 1000), (-1000, 1000))
lines!(ls1, plane_corners, color = :cyan, linewidth = 1.5, linestyle = :dash)

# 太阳风速度方向箭头 (MFA坐标, 固定长度600 km/s)
MVN_vsw_mfa = MVN_static2mfa * MVN_staticvswrough  # 静态坐标→MFA坐标
vsw_dir = MVN_vsw_mfa ./ norm(MVN_vsw_mfa)
vsw_arrow_len = 600.0  # 固定箭头长度 (km/s)
ps_vsw = Point3f(0.0, 0.0, 0.0)
pe_vsw = Point3f(vsw_dir .* vsw_arrow_len...)
linesegments!(ls1, [ps_vsw, pe_vsw], color = :red, linewidth = 8, overdraw = true)
arrows3d!(ls1, ps_vsw, pe_vsw, argmode = :endpoint,
    color = :red, shaftradius = 0.1, tipradius = 0.2, tailradius = 0.06, overdraw = true)

# # 反射方向箭头 (橙色, 固定长度600 km/s)
# ref_dir = v_reflected ./ norm(v_reflected)
# ref_arrow_len = 600.0  # 固定箭头长度 (km/s)
# ps_ref = Point3f(0.0, 0.0, 0.0)
# pe_ref = Point3f(ref_dir .* ref_arrow_len...)
# linesegments!(ls1, [ps_ref, pe_ref], color = :orange, linewidth = 6, overdraw = true)
# arrows3d!(ls1, ps_ref, pe_ref, argmode = :endpoint,
#     color = :orange, shaftradius = 0.08, tipradius = 0.15, tailradius = 0.05, overdraw = true)

# # 激波法向箭头 (白色)
# ps_shock = Point3f(0.0, 0.0, 0.0)
# pe_shock = Point3f(n_shock_mfa .* 600.0)
# linesegments!(ls1, [ps_shock, pe_shock], color = :white, linewidth = 6, overdraw = true)
# arrows3d!(ls1, ps_shock, pe_shock, argmode = :endpoint,
#     color = :white, shaftradius = 0.08, tipradius = 0.15, tailradius = 0.05, overdraw = true)

# 激波面 (垂直于法向的平面, 白色虚线边框)
# 构建平面内正交基
if abs(n_shock_mfa[1]) < 0.9
    aux_shock = [1.0, 0.0, 0.0]
else
    aux_shock = [0.0, 1.0, 0.0]
end
u_shock = cross(aux_shock, n_shock_mfa)
u_shock = u_shock ./ norm(u_shock)
v_shock = cross(n_shock_mfa, u_shock)

# 画激波面边框
shock_corners = [
    -1000.0 .* u_shock .- 1000.0 .* v_shock,
     1000.0 .* u_shock .- 1000.0 .* v_shock,
     1000.0 .* u_shock .+ 1000.0 .* v_shock,
    -1000.0 .* u_shock .+ 1000.0 .* v_shock,
    -1000.0 .* u_shock .- 1000.0 .* v_shock,  # 闭合
]
shock_corners_3f = [Point3f(c...) for c in shock_corners]
lines!(ls1, shock_corners_3f, color = :orange, linewidth = 2, linestyle = :dash, overdraw = true)

# # ── 右图: 2D截面展开图 ──
# Label(fig[1, 2, Top()], title_str * "\nSlice: n=$(round.(plane_normal, digits=2))", fontsize = 16, tellheight = true)
# ax2 = Axis(fig[1, 2],
#     xlabel = "s (km/s) ← $(round.(u_vec, digits=2))",
#     ylabel = "t (km/s) ← $(round.(v_vec, digits=2))",
#     aspect = 1)

# hm2 = heatmap!(ax2, s_vals, t_vals, eflux_slice',
#     colormap = :gist_earth, colorrange = (1e4, 1e8), colorscale = log10)

# # 标注MFA轴方向
# text!(ax2, "n̂ = $(round.(plane_normal, digits=2))", position = (-900, 900),
#     fontsize = 12, color = :cyan, font = :bold)
# text!(ax2, "û = $(round.(u_vec, digits=2))", position = (-900, 800),
#     fontsize = 11, color = :yellow)
# text!(ax2, "v̂ = $(round.(v_vec, digits=2))", position = (-900, 700),
#     fontsize = 11, color = :magenta)
# text!(ax2, "n̂_shock = $(round.(n_shock_mfa, digits=2))", position = (-900, 600),
#     fontsize = 11, color = :orange, font = :bold)

# Colorbar(fig[1, 3], hm2, label = "Differential Energy Flux")

save(joinpath(VDF_path, fig_name), fig)
println("Saved: ", joinpath(VDF_path, fig_name))

# # ═══════════════════════════════════════════════════════
# #  新图: 3D VDF 每个vpara切片投影到截面上 (MIP)
# # ═══════════════════════════════════════════════════════

# # 截面投影坐标轴 (与截面heatmap一致)
# n_proj = 200
# s_range = (-1000.0, 1000.0)
# t_range = (-1000.0, 1000.0)
# s_grid = range(s_range[1], s_range[2], length=n_proj)
# t_grid = range(t_range[1], t_range[2], length=n_proj)

# # 将3D MFA坐标投影到截面2D坐标 (s, t)
# # p = vperp1 * e1 + vpara * e2 + vperp2 * e3
# # s = dot(p, u_vec), t = dot(p, v_vec)
# function project_mfa_to_slice(vperp1, vpara_val, vperp2, u, v, e1, e2, e3)
#     p = vperp1 .* e1 .+ vpara_val .* e2 .+ vperp2 .* e3
#     s = dot(p, u)
#     t = dot(p, v)
#     return s, t
# end

# # 计算投影: 遍历每个vpara切片, 将其所有点投影到截面, 取最大值(MIP)
# println("Computing slice projection (MIP)...")
# n_vp = length(vpara)
# n_vp1 = length(vperp1)
# n_vp2 = length(vperp2)
# proj_mip = fill(-Inf64, n_proj, n_proj)  # Maximum Intensity Projection

# ds = step(s_grid)
# dt = step(t_grid)

# for ipara in 1:n_vp
#     for i1 in 1:n_vp1, i2 in 1:n_vp2
#         val = MVN_static_eflux_rect[i1, ipara, i2]
#         isnan(val) && continue
#         # 该点的3D MFA坐标
#         p = vperp1[i1] .* MVN_perpstatic1 .+ vpara[ipara] .* MVN_parastatic .+ vperp2[i2] .* MVN_perpstatic2
#         s_coord = dot(p, u_vec)
#         t_coord = dot(p, v_vec)
#         # 映射到投影网格索引
#         is_idx = round(Int, (s_coord - s_range[1]) / ds) + 1
#         it_idx = round(Int, (t_coord - t_range[1]) / dt) + 1
#         if 1 <= is_idx <= n_proj && 1 <= it_idx <= n_proj
#             proj_mip[is_idx, it_idx] = max(proj_mip[is_idx, it_idx], val)
#         end
#     end
#     mod(ipara, 5) == 0 && println("  vpara slice $ipara/$n_vp done")
# end

# # 将 -Inf 替换为 NaN (无数据区域)
# proj_mip[.!isfinite.(proj_mip)] .= NaN

# println("Projection done. Finite points: $(count(isfinite, proj_mip))/$(n_proj*n_proj)")

# # 将3D方向投影到截面2D坐标 (s, t)
# function project_to_plane(v3d, u, v)
#     s = dot(v3d, u)
#     t = dot(v3d, v)
#     return s, t
# end

# # 太阳风方向投影到截面
# vsw_proj_s, vsw_proj_t = project_to_plane(MVN_vsw_mfa, u_vec, v_vec)
# # 反射方向投影到截面
# ref_proj_s, ref_proj_t = project_to_plane(v_reflected, u_vec, v_vec)
# # 激波法向投影到截面
# shock_proj_s, shock_proj_t = project_to_plane(n_shock_mfa, u_vec, v_vec)

# # 固定箭头长度 600 km/s
# arrow_scale = 600.0 / max(norm(MVN_vsw_mfa), 1e-6)

# # 新figure
# fig2_name = "MAVEN_STATIC_H+_VDF_slice_proj_" * Dates.format(ts, "yyyymmddTHHMMSS") *
#     "-" * Dates.format(te, "HHMMSS") * ".png"

# fig2 = Figure(size = (800, 700), fontsize = 14, figure_padding = (5, 10, 9, 8))
# Label(fig2[1, 1, Top()],
#     "MAVEN STATIC H⁺ VDF — Slice Projection (MIP)\n" *
#     Dates.format(ts, "yyyymmdd HH:MM") * "-" * Dates.format(te, "HH:MM"),
#     fontsize = 15, tellheight = true)

# ax3 = Axis(fig2[2, 1],
#     xlabel = "û (km/s)",
#     ylabel = "v̂ (km/s)",
#     title = "n̂ = $(round.(plane_normal, digits=3))",
#     aspect = 1)

# # 截面投影 heatmap (MIP)
# # CairoMakie heatmap 需要转置: 数据[i,j] → x=j, y=i
# hm3 = heatmap!(ax3, collect(s_grid), collect(t_grid), proj_mip',
#     colormap = :gist_earth, colorrange = (1e4, 1e8), colorscale = log10)

# # 太阳风方向 (红色箭头, 固定600 km/s)
# vsw_s = vsw_proj_s * arrow_scale
# vsw_t = vsw_proj_t * arrow_scale
# arrows!(ax3, [0.0], [0.0], [vsw_s], [vsw_t],
#     color = :red, linewidth = 3, arrowsize = 20, lengthscale = 1.0)
# text!(ax3, "  Vsw", position = (vsw_s, vsw_t),
#     fontsize = 12, color = :red, font = :bold, align = (:left, :center))

# # 反射方向 (橙色箭头, 固定600 km/s)
# ref_s = ref_proj_s * arrow_scale
# ref_t = ref_proj_t * arrow_scale
# arrows!(ax3, [0.0], [0.0], [ref_s], [ref_t],
#     color = :orange, linewidth = 3, arrowsize = 20, lengthscale = 1.0)
# text!(ax3, "  Vref", position = (ref_s, ref_t),
#     fontsize = 12, color = :orange, font = :bold, align = (:left, :center))

# # 激波法向 (白色虚线, 固定600 km/s)
# shock_s = shock_proj_s * arrow_scale
# shock_t = shock_proj_t * arrow_scale
# if abs(shock_s) > 1.0 || abs(shock_t) > 1.0
#     lines!(ax3, [0.0, shock_s], [0.0, shock_t],
#         color = :white, linewidth = 2, linestyle = :dash)
#     text!(ax3, "  n̂_shock", position = (shock_s, shock_t),
#         fontsize = 11, color = :white, font = :bold, align = (:left, :center))
# end

# # MFA轴标注
# text!(ax3, "û = $(round.(u_vec, digits=2))", position = (-900, -900),
#     fontsize = 10, color = :gray)
# text!(ax3, "v̂ = $(round.(v_vec, digits=2))", position = (-900, -950),
#     fontsize = 10, color = :gray)

# Colorbar(fig2[2, 2], hm3, label = "Differential Energy Flux")

# save(joinpath(VDF_path, fig2_name), fig2)
# println("Saved: ", joinpath(VDF_path, fig2_name))

display(fig)
while isopen(fig.scene)
    sleep(0.1)
end