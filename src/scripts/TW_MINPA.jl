module TW_MINPA
using Glob, DelimitedFiles, DataFrames
using Dates, LinearAlgebra, Statistics, Interpolations
using ThreadsX, WENO4, Polyester, CoordinateTransformations
using Rotations
using Unitful
import PhysicalConstants.CODATA2018: m_p, m_e, e, ε_0, μ_0, c_0, G, k_B
using GLMakie

include("SphericalSplines.jl")
include("TW_load.jl")
using .SphericalSplines
using .TW_load

export compute_velocity_table, switch_mod
export TW_MINPA_eflux_interp, compute_vdf

const M_Mars = 0.64169e24u"kg"
const R_Mars = 3389.5u"km"
const au = 149597870.7u"km"

const energymod1 = [2.81, 3.548928, 4.482167, 5.660813, 7.149401, 9.029433, 11.40385,
    14.40264, 18.19001, 22.97332, 29.01447, 36.64422, 46.28032, 58.45036, 73.82067,
    93.23282, 117.7497, 148.7135, 187.8198, 237.2095, 299.587, 378.3675, 477.8644,
    603.5253, 762.2309, 962.6694, 1215.816, 1535.532, 1939.321, 2449.292, 3093.366,
    3906.809, 4934.157, 6231.661, 7870.361, 9939.98, 12553.83, 15855.03, 20024.33, 25290.0]
const energymod4 = [2.81, 3.246924, 3.751784, 4.335145, 5.009211, 5.788088, 6.688071, 7.727991,
    8.929607, 10.31806, 11.9224, 13.77621, 15.91825, 18.39336, 21.25332, 24.55798,
    28.37647, 32.7887, 37.88697, 43.77797, 50.58496, 58.45036, 67.53873, 78.04025,
    90.17464, 104.1958, 120.3971, 139.1175, 160.7487, 185.7433, 214.6243, 247.996,
    286.5566, 331.113, 382.5974, 442.087, 510.8266, 590.2544, 682.0324, 788.0808,
    910.6185, 1052.21, 1215.816, 1404.862, 1623.303, 1875.708, 2167.36, 2504.36,
    2893.76, 3343.707, 3863.617, 4464.366, 5158.525, 5960.618, 6887.428, 7958.346, 
    9195.78, 10625.62, 12277.79, 14186.84, 16392.74, 18941.63, 21886.84, 25290.0]
const θmod1 = [11.3, 33.8, 56.2, 78.7]
const ϕmod1 = [11.25, 33.75, 56.25, 78.75, 101.25, 123.75, 146.25, 168.75,
               191.25, 213.75, 236.25, 258.75, 281.25, 303.75, 326.25, 348.75]
const massmod1 = [1, 2, 4, 16, 28, 32, 44, 64]
const TW_minpa2sc = [0.0 0.0 1.0; -1.0 0.0 0.0; 0.0 -1.0 0.0]

function switch_mod(mod)
    if mod == 1
        mass = massmod1
        minpaphi = 180.0 .- ϕmod1
        minpatheta = -θmod1
        energy = energymod1
    elseif mod == 4
        mass = massmod1
        minpaphi = 180.0 .- ϕmod1
        minpatheta = -θmod1
        energy = energymod4
    end
    return mass, energy, minpaphi, minpatheta
end

function compute_velocity_table(mod)
    mass, energy, minpaphi, minpatheta = switch_mod(mod)
    v = zeros(Float64, length(mass), length(energy))
    for i in 1:length(energy)
        for j in 1:length(mass)
            v[j, i] = sqrt(energy[i] * (2.0 / ustrip(u"eV", m_p * c_0^2)) / mass[j]) *
                      ustrip(u"km/s", c_0)
        end
    end
    return v
end

function TW_MINPA_eflux_interp(v3, TW_minpass; mod=1)
    vs = SphericalFromCartesian()(v3) 
    vtot = vs.r
    lon = vs.θ
    lat = vs.ϕ

    TW_minpavelocity = compute_velocity_table(mod)
    mass, energy, minpaphi, minpatheta = switch_mod(mod)

    cond1 = rad2deg(lat) <= 0.0
    cond2 = TW_minpavelocity[1, 18] ≤ vtot ≤ TW_minpavelocity[1, 31]
    cond3 = minpaphi[15]-22.5/2.0 <= rad2deg(lon) <= minpaphi[11]+22.5/2.0
    
    if cond1 && cond2 && cond3
    #if (TW_minpavelocity[1, 18] ≤ vtot ≤ TW_minpavelocity[1, 31])
        # dir = v3 ./ norm(vtot)
        dir = v3 ./ norm(v3)
        logv = log10.(TW_minpavelocity[1, 18:31])
        logeflux = Vector{Float64}(undef, 14)
        for k in 18:31 # energy index
            S = TW_minpass[k-17]
            logeflux[k-17] = S(dir)
        end
        x = logv
        y = logeflux
        return 10.0^interpolate_weno4([log10(vtot)], x, y)[1]
    else
        return NaN64
    end
end

function compute_vdf(minpa_data, mag_data, ts, te;
                     vpara_r=(-800.0, 800.0), vperp1_r=(-800.0, 800.0),
                     vperp2_r=(-800.0, 800.0), vst=10.0,
                     external_basis=nothing,
                     label="")
    ind_mag = ThreadsX.findall(ts .<= mag_data[:epoch] .<= te)
    TW_Bmeanmso = dropdims(mean(Matrix(mag_data[:B][ind_mag, :]), dims=(1,)), dims=(1,))
    roll  = mean(filter(!isnan, mag_data[:roll][ind_mag]))  |> deg2rad
    pitch = mean(filter(!isnan, mag_data[:pitch][ind_mag])) |> deg2rad
    yaw   = mean(filter(!isnan, mag_data[:yaw][ind_mag]))   |> deg2rad

    TW_mso2scmean = RotXYZ(roll=roll, pitch=pitch, yaw=yaw)
    TW_sc2msomean = transpose(TW_mso2scmean)

    # mass, energy, minpaphi, minpatheta = switch_mod(mod)
    mod = minpa_data[:mod]
    TW_minpavelocity = compute_velocity_table(mod)
    ind_minpa = ThreadsX.findall(ts .<= minpa_data[:epoch] .<= te)

    TW_minpaefluxmean = dropdims(mean(minpa_data[:eflux][ind_minpa, :, :, :, :], 
        dims=(1,)), dims=(1,))

    # 中位数 / 四分位数
    TW_minpaefluxmid    = similar(TW_minpaefluxmean)
    TW_minpaefluxuquant = similar(TW_minpaefluxmean)
    TW_minpaefluxlquant = similar(TW_minpaefluxmean)
    for i in axes(TW_minpaefluxmean, 1), j in axes(TW_minpaefluxmean, 2),
        k in axes(TW_minpaefluxmean, 3), l in axes(TW_minpaefluxmean, 4)
        data = minpa_data[:eflux][ind_minpa, i, j, k, l]
        if length(data) > 0
            TW_minpaefluxmid[i, j, k, l]    = quantile(data, 0.50)
            TW_minpaefluxuquant[i, j, k, l] = quantile(data, 0.75)
            TW_minpaefluxlquant[i, j, k, l] = quantile(data, 0.25)
        else
            TW_minpaefluxmid[i, j, k, l]    = NaN
            TW_minpaefluxuquant[i, j, k, l] = NaN
            TW_minpaefluxlquant[i, j, k, l] = NaN
        end
    end

    # 清理可疑通道
    TW_minpaefluxmean[1, :, 4, :] .= 0.0
    TW_minpaefluxmean[1, 1:10, :, :] .= 0.0
    TW_minpaefluxmean[1, 16, :, :] .= 0.0

    # 粗略太阳风方向
    i_peak, j_peak = 1, 13
    ind_peak = argmax(TW_minpaefluxmean[1, j_peak, i_peak, :])
    TW_minpavswrough = CartesianFromSpherical()(
        CoordinateTransformations.Spherical(
            TW_minpavelocity[1, ind_peak],
            deg2rad(minpa_data[:phi][j_peak]),
            deg2rad(minpa_data[:theta][i_peak])))
    
    # MFA 基向量
    if external_basis === nothing
        TW_Bmeanminpa = TW_minpa2sc' * TW_sc2msomean' * TW_Bmeanmso
        TW_paraminpa   = TW_Bmeanminpa ./ norm(TW_Bmeanminpa)
        TW_perpminpa1  = cross(TW_paraminpa, TW_minpavswrough)
        TW_perpminpa1  = TW_perpminpa1 ./ norm(TW_perpminpa1)
        TW_perpminpa2  = cross(TW_perpminpa1, TW_paraminpa)

        angle_B_Vsw = rad2deg(atan(norm(cross(TW_Bmeanminpa, Array(TW_minpavswrough))),
                                    dot(TW_Bmeanminpa, Array(TW_minpavswrough))))
        println("  ∠(B, Vsw_rough) = $(round(angle_B_Vsw, digits=1))°")
    else
        TW_paraminpa   = external_basis.paraminpa
        TW_perpminpa1  = external_basis.perpminpa1
        TW_perpminpa2  = external_basis.perpminpa2
    end
    TW_minpa2mfa = transpose(hcat(TW_perpminpa1, TW_paraminpa, TW_perpminpa2))

    # 球面样条插值器
    TW_minpadirections = [
        CartesianFromSpherical()(CoordinateTransformations.Spherical(
            1.0, deg2rad(minpa_data[:phi][jj]), deg2rad(minpa_data[:theta][ii])))
        for jj in 1:16, ii in 1:4] |> vec

    TW_minpass = []
    Snan = SphericalSplines.InterpolationSpline(stack(TW_minpadirections; dims=1),
                                                 fill(NaN64, 64))
    for k in 18:31
        temp = [TW_minpaefluxmean[1, jj, ii, k] for jj in 1:16 for ii in 1:4]
        ind_valid = findall(temp .> 0.0)
        if !isempty(ind_valid)
            S = SphericalSplines.InterpolationSpline(
                stack(TW_minpadirections[ind_valid]; dims=1), log10.(temp[ind_valid]))
            push!(TW_minpass, S)
        else
            push!(TW_minpass, Snan)
        end
    end

    # 在选定 MFA 基向量上插值到三维速度网格
    vpara  = collect(range(vpara_r[1],  vpara_r[2],  step=vst))
    vperp1 = collect(range(vperp1_r[1], vperp1_r[2], step=vst))
    vperp2 = collect(range(vperp2_r[1], vperp2_r[2], step=vst))

    TW_minpa_efluxmean_rect = Array{Float64, 3}(undef, length(vperp1), length(vpara), length(vperp2))
    @batch for k in 1:length(vperp2)
        for j in 1:length(vpara)
            for i in 1:length(vperp1)
                v3 = Array(vpara[j] .* TW_paraminpa .+
                           vperp1[i] .* TW_perpminpa1 .+
                           vperp2[k] .* TW_perpminpa2)
                TW_minpa_efluxmean_rect[i, j, k] = TW_MINPA_eflux_interp(v3, TW_minpass, mod=mod)
            end
        end
    end

    # n_valid = count(!isnan, TW_minpa_efluxmean_rect)
    # n_total = length(TW_minpa_efluxmean_rect)
    # pct = round(n_valid / n_total * 100, digits=2)
    # println("  VDF 网格有效: $n_valid / $n_total ($pct%)")
    # if external_basis !== nothing
    #     println("  (在外部 MFA 坐标系中求值)")
    # end

    return (rect       = TW_minpa_efluxmean_rect,
            vpara      = vpara, vperp1 = vperp1, vperp2 = vperp2,
            perpminpa1 = TW_perpminpa1,
            paraminpa  = TW_paraminpa,
            perpminpa2 = TW_perpminpa2,
            efluxmean    = TW_minpaefluxmean,
            # paphi      = minpaphi,
            # patheta    = minpatheta,
            pvelocity  = TW_minpavelocity,
            ts = ts, te = te)
end                    

end