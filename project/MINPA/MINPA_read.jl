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

include("../../src/scripts/SphericalSplines.jl")
using .SphericalSplines

date = DateTime(2022, 08, 24, 00, 00, 00)
datejul = datetime2julian(date)
datestr = Dates.format.(date, "yyyymmdd")
ts = DateTime(2022, 08, 24, 08, 10, 00) 
te = DateTime(2022, 08, 24, 14, 00, 00)
# ts = DateTime(2022, 08, 24, 7, 30, 00) 
# te = DateTime(2022, 08, 24, 8, 00, 00)

# read TW data 
TW_file2c1hz = "E:/Data/Tianwen-1/MOMAG/2C.V03/2022/08/".*["TW1_MOMAG_MSO_01hz_20220824_2C_v03.dat"]
# TW_file2c1hz = glob("*$datestr*.dat", TW_datapath2c1hz)
TW_mag2c1hz = reduce(vcat, [identity.(DataFrame(readdlm(file, skipstart=19), :auto)) for file in TW_file2c1hz])
name = ["Time", "Sampling_Rate", "X_MSO", "Y_MSO", "Z_MSO", "Probe_Position_X_MSO", "Probe_Position_Y_MSO", "Probe_Position_Z_MSO", "Roll", "Pitch", "Yaw",  "Quality_Flags"]
rename!(TW_mag2c1hz, name)
TW_mag2c1hz[!, :Time] = map(x->DateTime(x[begin:end-4], DateFormat("y-m-dTH:M:S.s")), TW_mag2c1hz[!, :Time])
TW_mag2c1hz[!, :JulUT] = datetime2julian.(TW_mag2c1hz[!, :Time])
unique!(TW_mag2c1hz) #remove depulicate rows
sort!(TW_mag2c1hz) #sorting
#TW_mag2c1hz[!, [:Probe_Position_X_MSO, :Probe_Position_Y_MSO, :Probe_Position_Z_MSO]] = TW_mag2c1hz[!, [:Probe_Position_X_MSO, :Probe_Position_Y_MSO, :Probe_Position_Z_MSO]] ./R_Mars
println("TW_MAG 2C1Hz data read, number of rows: ", nrow(TW_mag2c1hz))

TW_magut1hz = TW_mag2c1hz[:, :JulUT].-datejul
TW_BMSO1hz = TW_mag2c1hz[:, [:X_MSO, :Y_MSO, :Z_MSO]]
ind0 = findall(x-> !isnan(x), TW_BMSO1hz[!, 1] ) 
ind1 = findall(x-> isnan(x), TW_BMSO1hz[!, 1]) 
for ib in 1:3 
    interp_linear = linear_interpolation(TW_magut1hz[ind0], TW_BMSO1hz[ind0, ib]; extrapolation_bc=Interpolations.Line())
    TW_BMSO1hz[ind1, ib] = interp_linear(TW_magut1hz[ind1])
end

ind = ThreadsX.findall(ts .<= TW_mag2c1hz[!, :Time] .<= te)
TW_Bmeanmso = dropdims(mean(Matrix(TW_BMSO1hz[ind, :]), dims=(1,)), dims=(1,))
roll  = mean(filter(!isnan, TW_mag2c1hz[ind, :Roll]))  |> deg2rad
pitch = mean(filter(!isnan, TW_mag2c1hz[ind, :Pitch])) |> deg2rad
yaw   = mean(filter(!isnan, TW_mag2c1hz[ind, :Yaw]))   |> deg2rad
TW_mso2scmean = RotXYZ(roll=roll, pitch=pitch, yaw=yaw)
TW_sc2msomean = transpose(TW_mso2scmean)

minpafmod1 = "E:/Data/Tianwen-1/MINPA/2B/2022/08/".*["HX1-Or_GRAS_MINPA-MOD1-DEF_SCI_N_20220824080616_20220824140203_01536_A.2B"]
# minpafmod1 = "E:/Data/Tianwen-1/MINPA/2B/2022/08/".*["HX1-Or_GRAS_MINPA-MOD4-MT_SCI_N_20220824073202_20220824080551_01536_A.2B"]
minpamod1 = reduce(vcat, [identity.(DataFrame(readdlm(file), :auto)) for file in minpafmod1])
massmod1 = [1,2,4,16,28,32,44,64]
energymod1 = [2.81, 3.548928, 4.482167, 5.660813, 7.149401, 9.029433, 11.40385, 
    14.40264, 18.19001, 22.97332, 29.01447, 36.64422, 46.28032, 58.45036, 73.82067, 
    93.23282, 117.7497, 148.7135, 187.8198, 237.2095, 299.587, 378.3675, 477.8644, 
    603.5253, 762.2305, 962.6694, 1215.816, 1535.532, 1939.321, 2449.292, 3093.366, 
    3906.809, 4934.157, 6231.661, 7870.361, 9939.98, 12553.83, 15855.03, 20024.33, 25290.0]
# energymod4 = [2.81, 3.246924, 3.751784, 4.335145, 5.009211, 5.788088, 6.688071, 7.727991,
#     8.929607, 10.31806, 11.9224, 13.77621, 15.91825, 18.39336, 21.25332, 24.55798,
#     28.37647, 32.7887, 37.88697, 43.77797, 50.58496, 58.45036, 67.53873, 78.04025,
#     90.17464, 104.1958, 120.3971, 139.1175, 160.7487, 185.7433, 214.6243, 247.996,
#     286.5566, 331.113, 382.5974, 442.087, 510.8266, 590.2544, 682.0324, 788.0808,
#     910.6185, 1052.21, 1215.816, 1404.862, 1623.303, 1875.708, 2167.36, 2504.36,
#     2893.76, 3343.707, 3863.617, 4464.366, 5158.525, 5960.618, 6887.428, 7958.346, 
#     9195.78, 10625.62, 12277.79, 14186.84, 16392.74, 18941.63, 21886.84, 25290.0]
energymod12 = [2.81, 3.267539, 3.799578, 4.418245, 5.137647, 5.974187, 6.946936, 8.078073,
    9.393388, 10.92287, 12.70139, 14.7695, 17.17435, 19.97077, 23.22251, 27.00373,
    31.40062, 36.51344, 42.45875, 49.37211, 57.41115, 66.75914, 77.62922, 90.26922,
    104.9673, 122.0587, 141.9329, 165.0432, 191.9164, 223.1653, 259.5023, 301.7558,
    350.8893, 408.023, 474.4595, 551.7135, 641.5465, 746.0065, 867.4753, 1008.722,
    1172.968, 1363.957, 1586.043, 1844.292, 2144.589, 2493.783, 2899.834, 3372]
θmod1 = [11.3, 33.8, 56.2, 78.7]
ϕmod1 = [11.25, 33.75, 56.25, 78.75, 101.25, 123.75, 146.25, 168.75, 191.25, 213.75, 236.25, 258.75, 281.25, 303.75, 326.25, 348.75]
nmod1 = length(massmod1) * length(ϕmod1) * length(θmod1) * length(energymod1)
minpamod1 = reduce(vcat, [identity.(DataFrame(readdlm(file), :auto)) for file in minpafmod1])
minpautmod1 = map(x->DateTime(x[begin:end-4], DateFormat("y-m-dTH:M:S.s")), minpamod1[!, 1])
minpajulianmod1 = datetime2julian.(minpautmod1)
minpadfimod1= reshape(Array(minpamod1[:,60:59+nmod1]), length(minpautmod1), length(massmod1), length(ϕmod1), length(θmod1), length(energymod1)) #时间 质量数 方位角 俯仰角 能道 #
ind =findall(x->x<0, minpadfimod1)
minpadfimod1[ind] .= 0.0
println("TW-MINPA data read, number of time points: ", length(minpautmod1))

TW_minpaeflux = minpadfimod1 .* reshape(energymod1, 1, 1, 1, 1, length(energymod1)) # convert to differential energy flux
TW_minpavelocity = zeros(Float64, length(massmod1), length(energymod1))
for i in 1:length(energymod1)
    for j in 1:length(massmod1)
        TW_minpavelocity[j, i] = sqrt(energymod1[i]*(2.0/ustrip(u"eV", m_p*c_0^2))/massmod1[j])*ustrip(u"km/s", c_0)
    end
end
#TW_minpaphi = 360.0 .- ϕmod1 .+ 180.0
TW_minpaphi = 180.0 .- ϕmod1
TW_minpatheta = -θmod1  # latitude

TW_minpa2sc = [0.0 0.0 1.0; -1.0 0.0 0.0; 0.0 -1.0 0.0] #TW-MINPA instrument to SC transformation matrix

ind = ThreadsX.findall(ts .<= minpautmod1 .<= te)
TW_minpaefluxmean = dropdims(mean(TW_minpaeflux[ind, :, :, :, :], dims=(1, )), dims=(1,))

TW_minpaefluxmid = similar(TW_minpaefluxmean)
TW_minpaefluxuquant = similar(TW_minpaefluxmean)
TW_minpaefluxlquant = similar(TW_minpaefluxmean)

for i in axes(TW_minpaefluxmean, 1), j in axes(TW_minpaefluxmean, 2), k in axes(TW_minpaefluxmean, 3), l in axes(TW_minpaefluxmean, 4)
    data = TW_minpaeflux[ind, i, j, k, l]
    TW_minpaefluxmid[i, j, k, l] = quantile(data, 0.5)
    TW_minpaefluxuquant[i, j, k, l] = quantile(data, 0.75)
    TW_minpaefluxlquant[i, j, k, l] = quantile(data, 0.25)
end

i, j = 1, 13
ind = argmax(TW_minpaefluxmean[1, j, i, :]./energymod1.^2)
ind = argmax(TW_minpaefluxmean[1, j, i, :])
#=
    3D spherical coordinates in the CoordinateTransformations library.

    There are many Spherical coordinate conventions and this library uses a somewhat exotic one. Given a vector v with Cartesian coordinates xyz, let v_xy = [x,y,0] be the orthogonal projection of v on the xy plane.

    r is the radius. It is given by norm(v, 2).
    θ is the azimuth. It is the angle from the x-axis to v_xy
    ϕ is the latitude. It is the angle from v_xy to v.
=#
TW_minpaefluxmean[1, :, 4, :] .= 0.0
TW_minpaefluxmean[1, 1:10, :, :] .= 0.0
TW_minpaefluxmean[1, 16, :, :] .= 0.0
TW_minpavswrough = CartesianFromSpherical()(CoordinateTransformations.Spherical(TW_minpavelocity[1, ind], deg2rad(TW_minpaphi[j]), deg2rad(TW_minpatheta[i])))

TW_minpavswrough_mso = TW_sc2msomean * TW_minpa2sc * TW_minpavswrough

TW_Bmeanminpa = TW_minpa2sc' * TW_sc2msomean' * TW_Bmeanmso  #magnetic field in MINPA coordinates at mean time
TW_paraminpa = TW_Bmeanminpa ./ norm(TW_Bmeanminpa) #unit vector of magnetic field in MINPA coordinates at mean time

# @show "angle between B and Vsw from MINPA", 
    # rad2deg(atan(norm(cross(TW_Bmeanminpa, Array(TW_minpavswrough))), dot(TW_Bmeanminpa, Array(TW_minpavswrough))))

TW_perpminpa1 = cross(TW_paraminpa, TW_minpavswrough) #perpendicular direction1: E=B×Vsw
TW_perpminpa1 = TW_perpminpa1 ./ norm(TW_perpminpa1)
TW_perpminpa2 = cross(TW_perpminpa1, TW_paraminpa) #perpendicular direction 2: N=E×B
TW_minpa2mfa = transpose(hcat(TW_perpminpa1, TW_paraminpa, TW_perpminpa2))

TW_minpadirections = [CartesianFromSpherical()(CoordinateTransformations.Spherical(1.0, deg2rad(TW_minpaphi[j]), deg2rad(TW_minpatheta[i]))) for j in 1:16 for i in 1:4] #direction unit vectors


function TW_MINPA_eflux_interp(v3, TW_minpass)
    vs = SphericalFromCartesian()(v3) 
    vtot = vs.r
    lon = vs.θ
    lat = vs.ϕ

    cond1 = rad2deg(lat) <= 0.0
    cond2 = TW_minpavelocity[1, 18] ≤ vtot ≤ TW_minpavelocity[1, 31]
    cond3 = TW_minpaphi[15]-22.5/2.0 <= rad2deg(lon) <= TW_minpaphi[11]+22.5/2.0
    
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

vpara = collect(range(-800.0, 800.0, step = 10.0))
vperp1 = collect(range(-800.0, 800.0, step = 10.0))
vperp2 = collect(range(-800.0, 800.0, step = 10.0))

TW_minpass = [] #interpolation splines
Snan = SphericalSplines.InterpolationSpline(stack(TW_minpadirections; dims=1), fill(NaN64, 64))
for k in 18:31
    temp = [TW_minpaefluxmean[1, j, i, k] for j in 1:16 for i in 1:4]
    ind = findall(temp .> 0.0)
    if !isempty(ind)
        S = SphericalSplines.InterpolationSpline(stack(TW_minpadirections[ind]; dims=1), log10.(temp[ind]))
        push!(TW_minpass, S)
    else
        push!(TW_minpass, Snan)
    end  
end

TW_minpa_efluxmean_rect = Array{Float64, 3}(undef, length(vperp1), length(vpara), length(vperp2))
@batch for k in 1:length(vperp2)
    for j in 1:length(vpara)
        for i in 1:length(vperp1)
            v3 = Array(vpara[j] .* TW_paraminpa .+ vperp1[i] .* TW_perpminpa1 .+ vperp2[k] .* TW_perpminpa2)
            TW_minpa_efluxmean_rect[i, j, k] = TW_MINPA_eflux_interp(v3, TW_minpass)
            #TW_minpa_eflux_rect[i, j, k] = TW_MINPA_eflux_npinterp(v3)
        end
    end
end

using GLMakie
GLMakie.activate!()

size_inches = (20, 9)
size_pt = 72 .* size_inches
fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 10, 9, 8))
title_str = Dates.format(ts, "yyyymmddTHH:MM") * "-" * Dates.format(te, "HH:MM")
Label(fig[1, 1, Top()], title_str, fontsize = 18, tellheight = true)
ls1 = LScene(fig[1, 1])

hm_ref = heatmap!(ls1,
    vperp1, vperp2,
    TW_minpa_efluxmean_rect[:, 1, :],   # 用第1个 vpara 切片
    transformation = (:xz, vpara[1]),
    colormap = :gist_earth,
    colorrange = (1e4, 1e8),
    colorscale = log10)

#for ipara in 80:10:120 #length(vpara)
for ipara in 1:length(vpara)
# for ipara in 81:3:90
#for ipara in 80:10:100
#for ipara in 108:3:111
    vparai = vpara[ipara]
    # heatmap!(ls1, vperp1, vperp2, TW_minpa_eflux_rect[:, ipara, :], transformation = (:xz, vparai), colormap = :gist_earth,  colorrange = (1e4, 1e8), colorscale = log10)
    hm = heatmap!(ls1, vperp1, vperp2, TW_minpa_efluxmean_rect[:, ipara, :], transformation = (:xz, vparai), 
        colormap = :gist_earth,  colorrange = (1e4, 1e8), colorscale = log10)
    if ipara == 1
        hm_ref = hm
    end
end

i, j = 1, 13
ind = argmax(TW_minpaefluxmean[1, j, i, :])
v1 = CartesianFromSpherical()(CoordinateTransformations.Spherical(TW_minpavelocity[1, ind], deg2rad(TW_minpaphi[j]), deg2rad(TW_minpatheta[i])))
v1 = CartesianFromSpherical()(CoordinateTransformations.Spherical(TW_minpavelocity[1, 22], deg2rad(TW_minpaphi[j]), deg2rad(TW_minpatheta[i]+22.5)))
ps = Point3f(0.0, 0.0, 0.0)
pe = Point3f(dot(v1, TW_perpminpa1), dot(v1, TW_paraminpa), dot(v1, TW_perpminpa2))
arrows3d!(ls1, ps, pe, argmode = :endpoint, color = :red, shaftradius = 0.02, tipradius = 0.05, tailradius = 0.02)

Label(fig[1, 2, Top()], title_str, fontsize = 18, tellheight = true)
ls2 = LScene(fig[1, 2])

for iperp1 in 81:1:81
    vperp1i = vperp1[iperp1]
    # heatmap!(ls2, vpara, vperp2, TW_minpa_eflux_rect[iperp1, :, :], transformation = (:yz, vperp1i), colormap = :gist_earth,  colorrange = (1e4, 1e8), colorscale = log10)
    heatmap!(ls2, vpara, vperp2, TW_minpa_efluxmean_rect[iperp1, :, :], transformation = (:yz, vperp1i), colormap = :gist_earth,  colorrange = (1e4, 1e8), colorscale = log10)
end

i, j = 1, 13
ind = argmax(TW_minpaefluxmean[1, j, i, :])
v1 = CartesianFromSpherical()(CoordinateTransformations.Spherical(TW_minpavelocity[1, ind], deg2rad(TW_minpaphi[j]), deg2rad(TW_minpatheta[i])))
v1 = CartesianFromSpherical()(CoordinateTransformations.Spherical(TW_minpavelocity[1, 22], deg2rad(TW_minpaphi[j]), deg2rad(TW_minpatheta[i]+22.5)))
ps = Point3f(0.0, 0.0, 0.0)
pe = Point3f(dot(v1, TW_perpminpa1), dot(v1, TW_paraminpa), dot(v1, TW_perpminpa2))
arrows3d!(ls2, ps, pe, argmode = :endpoint, color = :red, shaftradius = 0.02, tipradius = 0.05, tailradius = 0.02)

Colorbar(fig[1, 3], hm_ref, label = "Differential Energy Flux")

display(fig)
# 持续显示窗口，直到用户关闭
while isopen(fig.scene)
    sleep(0.1)
end
# save("test2.png", fig)