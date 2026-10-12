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

ts = DateTime(date) + Dates.Hour(7) + Dates.Minute(26)
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
MVN_statictheta = MVN_static[:theta] |> Array #-45 to 45 deg, corresponding to latitude in the STATIC/APP coordinate system
MVN_staticphi = MVN_static[:phi] |> Array #-180 to 180 deg, corresponding to longitude in the STATIC/APP coordinate system
MVN_staticdenergy = MVN_static[:denergy] |> Array
MVN_staticdtheta = MVN_static[:dtheta] |> Array
MVN_staticdphi = MVN_static[:dphi] |> Array
MVN_staticdomega = MVN_static[:domega] |> Array
MVN_staticvelocity = sqrt.(MVN_staticenergy.*(2.0/ustrip(u"eV", m_p*c_0^2))./MVN_staticmass).*ustrip(u"km/s", c_0)
MVN_staticpos_sc_mso = MVN_static[:pos_sc_mso] |> Array
MVN_staticpos_sc_mso /= ustrip(u"km", R_Mars)

MVN_staticquat_mso = MVN_static[:quat_mso] |> Array #Quaternion elements to rotate from STATIC coordinates (same as APP coordinates) to MSO coordinates (NUM_DISTS, 3)
MVN_staticquat_sc = MVN_static[:quat_sc] |> Array #Quaternion elements to rotate from STATIC coordinates (same as APP coordinates) to SC coordinates (NUM_DISTS, 3)   
MVN_staticmsotrans = [QuatRotation(quat) for quat in eachcol(MVN_staticquat_mso')] # active (right handed) rotations where applicable.
MVN_staticind = ThreadsX.findall(ts .<= MVN_statictime .<= te)
# @show minimum(MVN_staticswp[MVN_staticind]), maximum(MVN_staticswp[MVN_staticind])
MVN_staticswpind = MVN_staticswp[MVN_staticind[1]]
MVN_staticefluxmean = dropdims(mean(MVN_staticeflux[MVN_staticind, :, :, :], dims = (1, ) ), dims=(1, ))
MVN_staticpsdmean = ustrip(u"cm^-3 * km^-3 * s^3", 0.5 * 1u"cm^-2 * s^-1" * (m_p/1u"eV")^2)  .* MVN_staticefluxmean .* (MVN_staticmass[:,:,:,MVN_staticswpind,1]./ MVN_staticenergy[:,:,:,MVN_staticswpind,1]).^2

# ===========================四分位数=======================================
size(MVN_staticefluxmean)
size(MVN_staticmass[:,:,:,MVN_staticswpind,1])
size(MVN_staticenergy[:,:,:,MVN_staticswpind,1])
MVN_staticefluxmid = similar(MVN_staticefluxmean)
MVN_staticefluxuquant = similar(MVN_staticefluxmean)
MVN_staticefluxlquant = similar(MVN_staticefluxmean)
for i in axes(MVN_staticeflux, 2), j in axes(MVN_staticeflux, 3), k in axes(MVN_staticeflux, 4)
    data = MVN_staticeflux[MVN_staticind, i, j, k]
    MVN_staticefluxmid[i, j, k] = quantile(data, 0.5)
    MVN_staticefluxuquant[i, j, k] = quantile(data, 0.75)
    MVN_staticefluxlquant[i, j, k] = quantile(data, 0.25)
end
i, j =2, 14
MVN_staticenergyl = zeros(32)
MVN_staticenergyu = zeros(32)
for ie in 1:32
    R = MVN_staticdenergy[1, (j-1)*4+i, ie, MVN_staticswpind, 1]./MVN_staticenergy[1, (j-1)*4+i, ie, MVN_staticswpind, 1]
    ratio = 0.5R+sqrt(1.0+0.25R^2)
    MVN_staticenergyl[ie] = MVN_staticenergy[1, (j-1)*4+i, ie, MVN_staticswpind, 1]/ratio
    MVN_staticenergyu[ie] = MVN_staticenergy[1, (j-1)*4+i, ie, MVN_staticswpind, 1]*ratio
end
MVN_staticenergygridH = vcat(MVN_staticenergyu, MVN_staticenergyl[end])
MVN_staticvelocitygridH = reverse(sqrt.(MVN_staticenergygridH.*(2.0/ustrip(u"eV", m_p*c_0^2))/1.03).*ustrip(u"km/s", c_0))
MVN_staticphigridH = [MVN_staticphi[1, (j-1)*4+1, 7, MVN_staticswpind, 1]-MVN_staticdphi[1, (j-1)*4+1, 7, MVN_staticswpind, 1]*0.5 for j in 1:16]
MVN_staticphigridH = vcat(MVN_staticphigridH, MVN_staticphigridH[begin]+360.0)
MVN_staticthetau = 90.0.-[MVN_statictheta[1, i, 7, MVN_staticswpind, 1]-MVN_staticdtheta[1, i, 7, MVN_staticswpind, 1]*0.5 for i in 1:4]
MVN_staticthetal = 90.0.-[MVN_statictheta[1, i, 7, MVN_staticswpind, 1]+MVN_staticdtheta[1, i, 7, MVN_staticswpind, 1]*0.5 for i in 1:4]
MVN_staticthetagridH = vcat(MVN_staticthetal, MVN_staticthetau[end])
# ===============================================================================
MVN_staticindmean = round(Int64, mean(MVN_staticind))
MVN_statictimemean = MVN_statictime[MVN_staticindmean]
MVN_staticjulianmean = MVN_staticjulian[MVN_staticindmean]
MVN_static2msomean = QuatRotation(MVN_staticquat_mso[MVN_staticindmean, :]) # STATIC instrument to MSO transformation matrix at mean time
MVN_Bmeanstatic = MVN_static2msomean' * MVN_Bmeanmso #magnetic field in STATIC coordinates at mean time
MVN_parastatic = MVN_Bmeanstatic ./ norm(MVN_Bmeanstatic) #unit vector of magnetic field in STATIC coordinates at mean time

# ====================太阳风方向==========================
swi, swj =2, 14
ind = argmax(MVN_staticefluxmean[1, (swj-1)*4+swi, :] ./ MVN_staticenergy[1, 1, :, MVN_staticswpind, 1].^2)
MVN_staticvswrough = CartesianFromSpherical()(
    CoordinateTransformations.Spherical(MVN_staticvelocity[1, (swj-1)*4+swi, ind, MVN_staticswpind, 1], 
    deg2rad(MVN_staticphi[1, (swj-1)*4+swi, ind, MVN_staticswpind, 1]), 
    deg2rad(MVN_statictheta[1, (swj-1)*4+swi, ind, MVN_staticswpind, 1])))
MVN_staticdirections = [CartesianFromSpherical()(CoordinateTransformations.Spherical(1.0, deg2rad(MVN_staticphi[1, ij, ind, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, ij, ind, MVN_staticswpind, 1]))) for ij in 1:64] #direction unit vectors
MVN_staticss = [] #interpolation splines
Snan = SphericalSplines.InterpolationSpline(stack(MVN_staticdirections; dims=1), fill(-30.0, 64))
for k in 1:32
    local ind = findall(MVN_staticefluxmean[1, :, k] .> 0.0)
    if !isempty(ind)
        S = SphericalSplines.InterpolationSpline(stack(MVN_staticdirections[ind]; dims=1), log10.(MVN_staticefluxmean[1, ind, k]))
        push!(MVN_staticss, S)
    else
        push!(MVN_staticss, Snan)
    end  
end
MVN_perpstatic1 = cross(MVN_parastatic, MVN_staticvswrough) #perpendicular direction1: E=B×Vsw
MVN_perpstatic1 = MVN_perpstatic1 ./ norm(MVN_perpstatic1)
MVN_perpstatic2 = cross(MVN_perpstatic1, MVN_parastatic) #perpendicular direction 2: N=E×B
MVN_static2mfa = transpose(hcat(MVN_perpstatic1, MVN_parastatic, MVN_perpstatic2))

# 预计算: 哪些能道(6~32)有有效样条（排除Snan），避免-30.0污染WENO4插值
valid_channels = [k for k in 6:32 if MVN_staticss[k] !== Snan]
valid_logv = log10.(MVN_staticvelocity[1, 1, valid_channels, MVN_staticswpind, 1])
# println("Valid WENO4 channels: $(length(valid_channels))/27 (indices $(valid_channels[1])..$(valid_channels[end]))")

function MVN_STATIC_eflux_interp(v3)
    vs = SphericalFromCartesian()(v3) 
    vtot = vs.r
    lon = vs.θ
    lat = vs.ϕ
    v_lo = MVN_staticvelocity[1, 1, valid_channels[1], MVN_staticswpind, 1]   # 最高速度
    v_hi = MVN_staticvelocity[1, 1, valid_channels[end], MVN_staticswpind, 1] # 最低速度
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

vpara = collect(range(-1000.0, 1000.0, step = 10.0))
vperp1 = collect(range(-1000.0, 1000.0, step = 10.0))
vperp2 = collect(range(-1000.0, 1000.0, step = 10.0))
MVN_static_eflux_rect = Array{Float64, 3}(undef, length(vperp1), length(vpara), length(vperp2))
@batch for k in 1:length(vperp2)
    for j in 1:length(vpara)
        for i in 1:length(vperp1)
            v3 = Array(vpara[j] .* MVN_parastatic .+ vperp1[i] .* MVN_perpstatic1 .+ vperp2[k] .* MVN_perpstatic2)
            MVN_static_eflux_rect[i, j, k] = MVN_STATIC_eflux_interp(v3)
            #MVN_static_eflux_rect[i, j, k] = MVN_STATIC_eflux_npinterp(v3)
        end
    end
end

# nan_count = count(isnan, MVN_static_eflux_rect)
# total = length(MVN_static_eflux_rect)
# finite_count = total - nan_count
# println("Grid: $total points, $finite_count finite, $nan_count NaN ($(round(100*nan_count/total, digits=1))%)")

# if finite_count > 0
#     finite_vals = MVN_static_eflux_rect[.!isnan.(MVN_static_eflux_rect)]
#     println("Finite value range: [$(minimum(finite_vals)), $(maximum(finite_vals))]")
#     @show maximum(finite_vals)
# else
#     println("ERROR: All grid points are NaN!")
# end

using GLMakie
GLMakie.activate!()

size_inches = (20, 9)
size_pt = 72 .* size_inches
fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 10, 9, 8))
title_str = Dates.format(ts, "yyyymmddTHH:MM") * "-" * Dates.format(te, "HH:MM")
fig_name = "MAVEN_STATIC_H+_VDF_3D_" * Dates.format(ts, "yyyymmddTHHMMSS") *
    "-" * Dates.format(te, "HHMMSS") * ".png"
Label(fig[1, 1, Top()], title_str, fontsize = 18, tellheight = true)
ls1 = LScene(fig[1, 1])

hm_ref = heatmap!(ls1,
    vperp1, vperp2,
    MVN_static_eflux_rect[:, 1, :],   # 用第1个 vpara 切片
    transformation = (:xz, vpara[1]),
    colormap = :gist_earth,
    colorrange = (1e4, 1e8),
    colorscale = log10)

for ipara in 1:length(vpara)
    vparai = vpara[ipara]
    hm = heatmap!(ls1, vperp1, vperp2, MVN_static_eflux_rect[:, ipara, :], transformation = (:xz, vparai), 
        colormap = :gist_earth,  colorrange = (1e4, 1e8), colorscale = log10)
    if ipara == 1
        hm_ref = hm
    end
end
# i, j = 1, 13
# ind = argmax(MVN_staticefluxmean[1, j, i, :])
# v1 = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, ind, MVN_staticswpind, 1], 
#     deg2rad(MVN_staticphi[1, (j-1)*4+i, ind, MVN_staticswpind, 1]),
#     deg2rad(MVN_statictheta[1, (j-1)*4+i, ind, MVN_staticswpind, 1])))
# # v1 = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, 22, MVN_staticswpind, 1], 
# #     deg2rad(MVN_staticphi[1, (j-1)*4+i, ind, MVN_staticswpind, 1]),
# #     deg2rad(MVN_statictheta[1, (j-1)*4+i, ind, MVN_staticswpind, 1]+22.5)))
# ps = Point3f(0.0, 0.0, 0.0)
# pe = Point3f(dot(v1, MVN_perpstatic1), dot(v1, MVN_parastatic), dot(v1, MVN_perpstatic2))
# arrows3d!(ls1, ps, pe, argmode = :endpoint, color = :red, shaftradius = 0.02, tipradius = 0.05, tailradius = 0.02)

Label(fig[1, 2, Top()], title_str, fontsize = 18, tellheight = true)
ls2 = LScene(fig[1, 2])

for iperp1 in 81:1:81
    vperp1i = vperp1[iperp1]
    # heatmap!(ls2, vpara, vperp2, TW_minpa_eflux_rect[iperp1, :, :], transformation = (:yz, vperp1i), colormap = :gist_earth,  colorrange = (1e4, 1e8), colorscale = log10)
    heatmap!(ls2, vpara, vperp2, MVN_static_eflux_rect[iperp1, :, :], transformation = (:yz, vperp1i), colormap = :gist_earth,  colorrange = (1e4, 1e8), colorscale = log10)
end

# i, j = 1, 13
# ind = argmax(MVN_staticefluxmean[1, j, i, :])
# v1 = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, ind, MVN_staticswpind, 1], 
#     deg2rad(MVN_staticphi[1, (j-1)*4+i, ind, MVN_staticswpind, 1]),
#     deg2rad(MVN_statictheta[1, (j-1)*4+i, ind, MVN_staticswpind, 1])))
# # v1 = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, 22, MVN_staticswpind, 1], 
# #     deg2rad(MVN_staticphi[1, (j-1)*4+i, ind, MVN_staticswpind, 1]),
# #     deg2rad(MVN_statictheta[1, (j-1)*4+i, ind, MVN_staticswpind, 1]+22.5)))
# ps = Point3f(0.0, 0.0, 0.0)
# pe = Point3f(dot(v1, MVN_perpstatic1), dot(v1, MVN_parastatic), dot(v1, MVN_perpstatic2))
# arrows3d!(ls2, ps, pe, argmode = :endpoint, color = :red, shaftradius = 0.02, tipradius = 0.05, tailradius = 0.02)

Colorbar(fig[1, 3], hm_ref, label = "Differential Energy Flux")

display(fig)
# 持续显示窗口，直到用户关闭
while isopen(fig.scene)
    sleep(0.1)
end