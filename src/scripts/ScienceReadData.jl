using Glob  # find files with a specific pattern
using DelimitedFiles, DataFrames, EzXML #read txt  and xml and then convert them to DataFrames
using CDFDatasets
using Dates
#using LaTeXStrings
#using ContinuousWavelets  #wavelet power
using Interpolations # interpolation
using PCHIPInterpolation
using StaticArrays #faster array
#using VMD # VMD decomposition
using DSP
using AbstractFFTs
#using Deconvolution
using LinearAlgebra
using ForwardDiff
using JLD2
using CairoMakie, GeometryBasics
using FortranFiles
using MAT
using Ephemerides # A pure Julia binary JPL ephemeris reader
using Printf
using SpecialFunctions
using Statistics
using Unitful
using CSV
using Rotations
using CoordinateTransformations
using CoordRefSystems
using Meshes
using AngleBetweenVectors
using LsqFit
using CurveFit
using WENO4
using Polyester
using ThreadsX
#import XLSX


import PhysicalConstants.CODATA2018: m_p, m_e, e, ε_0, μ_0, c_0, G, k_B
const M_Mars = 0.64169e24u"kg"
const R_Mars = 3389.5u"km"  #m
const au = 149597870.7u"km"


CairoMakie.activate!()
#set_theme!(theme_black())
set_theme!()


#read solar cycle data
scfile = "omni2_27day_NTS5PQs1Ow.lst"
sc = identity.(DataFrame(readdlm(scfile, skipstart=0), :auto))
name = ["year", "doy", "hour", "f107", "lyman"]
rename!(sc, name)
#Solar index F10.7 ( sfu = 10-22W m-2 Hz-1)
#Lyman Alpha Solar Index, W/m^2
sc[!, :time] = DateTime.(sc[!, :year]) .+ Day.(sc[!, :doy] .- 1)
sc[!, :julian] = datetime2julian.(sc[!, :time])

#read solar wind OMNI data
#omnifile = "OMNI2_H0_MRG1HR_4095501.txt"
omnifile = "OMNI2_H0_MRG1HR_1232574.txt"
#omnifile = "OMNI_HRO_1MIN_1194063.txt"
omni = identity.(DataFrame(readdlm(omnifile, skipstart=73), :auto))
#omnifile = "OMNI_HRO2_5MIN_2020135.txt"
#omni = identity.(DataFrame(readdlm(omnifile, skipstart=93), :auto))
name = ["date", "time", "Bt", "Tp", "Np", "Vp"]
rename!(omni, name)
omni[!, :UT] = DateTime.(omni[!, :date].*"T".*omni[!, :time], DateFormat("d-m-yTH:M:S.s"))
omni[!, :julian] = datetime2julian.(omni[!, :UT])
omni[!, :Tp] = omni[!, :Tp]./1.1604e4 #convert from K to eV
ind = findall(x->x>9000.0, omni[!, :Vp])
omni[ind, [:Bt, :Tp, :Np, :Vp]] .= NaN64

#read ACE data
ACE_mfifile = "ac_k0s_mfi_20220514000000_20220612000000_cdaweb.cdf"
ACE_mfi = cdfopen(ACE_mfifile)
ACE_mfitime = ACE_mfi["Epoch"] |> Array
ACE_mfijulian = datetime2julian.(DateTime.(ACE_mfitime))
ACE_mfiBt = ACE_mfi["Magnitude"] |> Array

ACE_swefile = "ac_k0s_swe_20220514000000_20220612000000_cdaweb.cdf"
ACE_swe = cdfopen(ACE_swefile)
ACE_swetime = ACE_swe["Epoch"] |> Array
ACE_swejulian = datetime2julian.(DateTime.(ACE_swetime))
ACE_sweNp = ACE_swe["Np"] |> Array
ACE_sweVp = ACE_swe["Vp"] |> Array
ACE_sweTp = ACE_swe["Tpr"] |> Array

#read WIND data
WIND_mfifile = "wi_h0s_mfi_20220514000030_20220611235930_cdaweb.cdf"
WIND_mfi = cdfopen(WIND_mfifile)
WIND_mfitime = WIND_mfi["Epoch"] |> Array
WIND_mfijulian = datetime2julian.(DateTime.(WIND_mfitime))
WIND_mfiBt = WIND_mfi["BF1"] |> Array
WIND_mfiind = ThreadsX.findall(x->x>0.0, WIND_mfiBt)

WIND_swefile = "wind_h1s_swe_20220514000020_20220611235809_cdaweb.cdf"
WIND_swe = cdfopen(WIND_swefile)
WIND_swetime = WIND_swe["Epoch"] |> Array
WIND_swejulian = datetime2julian.(DateTime.(WIND_swetime))
WIND_sweNa = WIND_swe["Alpha_Na_nonlin"] |> Array # He++ number density
WIND_sweNp = WIND_swe["Proton_Np_moment"] |> Array
WIND_sweVp = WIND_swe["Proton_V_moment"] |> Array
WIND_sweTp = ustrip.(u"eV", (Vector(WIND_swe["Proton_W_moment"]).*1u"km/s") .^2.0 .* m_p ./2.0)
WIND_sweind = ThreadsX.findall(x->x<9999.99, WIND_sweNp)

#Mars orbit and locations
eph = EphemerisProvider("E:\\JPLEphem\\" .* ["de440s.bsp", "mar097.bsp", "jup346.bsp", "sat452.bsp", "ura117.bsp", "nep103.bsp", "plu060.bsp"])
date0 = DateTime(2022, 06, 12, 0, 0, 0)
time0 = (datetime2julian(date0) - datetime2julian(DateTime(2000, 01, 01, 12, 00, 00)))*86400
time = time0 .+ 86400 .* (1:700)
sunpos = [ephem_vector3(eph, 0, 10, t) for t in time] |> stack .|> (x->x*1u"km")
merpos = [ephem_vector3(eph, 1, 199, t)+ephem_vector3(eph, 0, 1, t) for t in time] |> stack .|> (x->x*1u"km")
venpos = [ephem_vector3(eph, 2, 299, t)+ephem_vector3(eph, 0, 2, t) for t in time] |> stack .|> (x->x*1u"km")
earpos = [ephem_vector3(eph, 3, 399, t)+ephem_vector3(eph, 0, 3, t) for t in time] |> stack .|> (x->x*1u"km")
marpos = [ephem_vector3(eph, 4, 499, t)+ephem_vector3(eph, 0, 4, t) for t in time] |> stack .|> (x->x*1u"km")
juppos = [ephem_vector3(eph, 5, 599, t)+ephem_vector3(eph, 0, 5, t) for t in time] |> stack .|> (x->x*1u"km")
satpos = [ephem_vector3(eph, 6, 699, t)+ephem_vector3(eph, 0, 6, t) for t in time] |> stack .|> (x->x*1u"km")
urapos = [ephem_vector3(eph, 7, 799, t)+ephem_vector3(eph, 0, 7, t) for t in time] |> stack .|> (x->x*1u"km")
neppos = [ephem_vector3(eph, 8, 899, t)+ephem_vector3(eph, 0, 8, t) for t in time] |> stack .|> (x->x*1u"km")
plupos = [ephem_vector3(eph, 9, 999, t)+ephem_vector3(eph, 0, 9, t) for t in time] |> stack .|> (x->x*1u"km")

merpos = (merpos .- sunpos) ./au
venpos = (venpos .- sunpos) ./au
earpos = (earpos .- sunpos) ./au
marpos = (marpos .- sunpos) ./au
juppos = (juppos .- sunpos) ./au
satpos = (satpos .- sunpos) ./au
urapos = (urapos .- sunpos) ./au
neppos = (neppos .- sunpos) ./au
plupos = (plupos .- sunpos) ./au


r = dropdims(sum(x->x^2, marpos, dims=(1,)), dims=(1,))
marperihelion, pind = findmin(r) 
maraphelion, aind = findmax(r)
x = marpos[:, pind]
x = x ./ norm(x)
yt = marpos[:, aind]
z = cross(x, yt)
z = z ./ norm(z)
y = cross(z, x)

trm = [x y z]

merpos = trm' * merpos
venpos = trm' * venpos
earpos = trm' * earpos
marpos = trm' * marpos
juppos = trm' * juppos
satpos = trm' * satpos
urapos = trm' * urapos
neppos = trm' * neppos
plupos = trm' * plupos

function MarPos(ut)
    time = (datetime2julian(ut) - datetime2julian(DateTime(2000, 01, 01, 12, 00, 00)))*86400
    return trm' * (ephem_vector3(eph, 4, 499, time)+ephem_vector3(eph, 0, 4, time) - ephem_vector3(eph, 0, 10, time))*1u"km" / au
end

function EarthPos(ut)
    time = (datetime2julian(ut) - datetime2julian(DateTime(2000, 01, 01, 12, 00, 00)))*86400
    return trm' * (ephem_vector3(eph, 3, 399, time)+ephem_vector3(eph, 0, 3, time) - ephem_vector3(eph, 0, 10, time))*1u"km" / au
end

marp = stack(MarPos.(omni[!, :UT]))
earthp = stack(EarthPos.(omni[!, :UT]))

#Edberg, N. J. T., M. Lester, S. W. H. Cowley, and A. I. Eriksson (2008), Statistical analysis of the location of the Martian magnetic pileup boundary and bow shock and the influence of crustal magnetic fields, J. Geophys. Res., 113, A08206, doi:10.1029/2008JA013096.
#bow-shock model
function bowshock(xshock)
    xF = 0.55 # R_Mars
    ϵ = 1.05
    L = 2.10 # R_Mars
    rSD = 1.58
    temp = (ϵ^2-1.0)*(xshock-xF)^2-2ϵ*L*(xshock-xF)+L^2
    if temp>=0 
        return sqrt(temp)
    else
        return NaN64
    end
end
#magnetopause model
function magnetopause(xmp)
    rSD = 1.33
	xF = 0.86
	ϵ = 0.92
	L = 0.90
    temp = (ϵ^2-1.0)*(xmp-xF)^2-2ϵ*L*(xmp-xF)+L^2
    if temp>=0 
        return sqrt(temp)
    else
        return NaN64
    end
end
#minimum-variance-analysis
function MVA(magbc)
    nb = length(magbc[:, 1])
    bm = mean(magbc, dims=1)
    mμν = magbc' * magbc ./ nb - bm' * bm
    return eigen(mμν)
end

#calculate the norminal locations of bow shock and magnetopause
np = 5001
xshock = range(-10, 3, length=np)
ryzshock = bowshock.(xshock)
ind = findall(isfinite, ryzshock)
xshock = [xshock[ind]; reverse(xshock[ind])]
ryzshock = [ryzshock[ind]; -reverse(ryzshock[ind])]
xmp = range(-10, 3, length=np)
ryzmp = magnetopause.(xmp)
ind = findall(isfinite, ryzmp)
xmp = [xmp[ind]; reverse(xmp[ind])]
ryzmp = [ryzmp[ind]; -reverse(ryzmp[ind])]

date = DateTime(2022, 06, 12, 00, 00, 00)
datejul = datetime2julian(date)
datestr = Dates.format.(date, "yyyymmdd")
year = Dates.format(date, "yyyy")
month = Dates.format(date, "mm")

#read MVN Statistics
#data=load("Statistics_data_output.jld2")
#data=load("Statistics_data_output_20240820.jld2")
data=load("Statistics_data_output_20240822.jld2")
X = collect(get(data, "X", NaN64))
Y = collect(get(data, "Y", NaN64))
MVN_B = get(data, "B_MAVEN", NaN64) 
MVN_V = collect(get(data, "v_flow", NaN64)) 
MVN_VA = collect(get(data, "v_angle", NaN64)) 
MVN_N = get(data, "H_density", NaN64) 
ind = findall(x->x<1e-7, MVN_N)
MVN_N[ind].= NaN64 
ind = findall(x->x<1e-7, MVN_B)
MVN_B[ind].= NaN64 
ind = findall(x->x<1e-7, MVN_V)
MVN_V[ind].= NaN64 
ind = findall(x->x<1e-7, MVN_VA)
MVN_VA[ind].= NaN64 
MVN_Mₐ = NoUnits.(MVN_V.*1u"km/s" ./ (MVN_B.*1u"nT" ./ sqrt.(μ_0*m_p .* MVN_N.*1u"cm^-3")))
for (j, y0) in enumerate(Y)
    for (i, x0) in enumerate(X)
        ryz = magnetopause(x0)
        if !isnan(ryz) && y0 < ryz
            MVN_Mₐ[i, j] = NaN64
            MVN_N[i, j] = NaN64
            MVN_B[i, j] = NaN64
            MVN_V[i, j] = NaN64
            MVN_VA[i, j] = NaN64
        end
    end
end

function velocity(x, y, X, Y, xmp, ryzmp, MVN_VA)
    ryz = bowshock(x)
    if isnan(ryz) || y > ryz
        angle = π
    else 
        tmp = copy(MVN_VA)
        interp_linear = linear_interpolation((X, Y), tmp)
        angle = deg2rad(max(min(interp_linear(x, y), 180.0), 0.0)) 
        if isnan(angle)
            rt = @. abs(sqrt(x^2+y^2) - sqrt(xmp^2+ryzmp^2))
            minrt = minimum(rt)
            local ind = findall(x->isequal(x, minrt), rt)
            ind = ind[1]
            angle = atan(ForwardDiff.derivative(magnetopause, xmp[ind]))+π
        end
    end
    return angle 
end
vel(x, y) = velocity(x, y, X, Y, xmp, ryzmp, MVN_VA)

function stream(x, y)
    lx = [x]
    ly = [y]
    va = [vel(x, y)*1.0]
    xt = copy(x)
    yt = copy(y)

    Δr= 0.05
    while xt <= 10.0 && xt >= -5.0 && yt <= 4.0 && yt >= 0.0
        tp = vel(xt, yt)
        push!(va, tp)
        xt = xt - Δr*cos(tp)
        yt = yt - Δr*sin(tp)
        push!(lx, xt)
        push!(ly, yt)
    end
    xt = copy(x)
    yt = copy(y)
    while xt <= 10.0 && xt >= -5.0 && yt <= 4.0 && yt >= 0.0
        tp = vel(xt, yt)
        pushfirst!(va, tp)
        xt = xt + Δr*cos(tp)
        yt = yt + Δr*sin(tp)
        pushfirst!(lx, xt)
        pushfirst!(ly, yt)
    end
   return lx, ly, va 

end

# Define the Maxwellian distribution function for the solar wind core
#swmaxwellian(N, vthpara, vthperp, vdpara, vdperp1, vdperp2, v3) = N / (2*π)^1.5 / vthpara / vthperp^2 * exp(-0.5 * ((dot(v3, TW_perpminpa2)-vdperp2)/vthperp)^2) * exp(-0.5 * ((dot(v3, TW_perpminpa1)-vdperp1)/vthperp)^2) * exp(-0.5 * ((dot(v3, TW_paraminpa)-vdpara)/vthpara)^2) * norm(v3)^4 / (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
function swmaxwellian(N, vthpara, vthperp, vdpara, vdperp1, vdperp2, inst2mfa, v3) 
    v3mfa = inst2mfa * v3
    return N / (2*π)^1.5 / vthpara / vthperp^2 * exp(-0.5 * ((v3mfa[3]-vdperp2)/vthperp)^2) * exp(-0.5 * ((v3mfa[1]-vdperp1)/vthperp)^2) * exp(-0.5 * ((v3mfa[2]-vdpara)/vthpara)^2) * norm(v3)^4 / (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
end

#Define the Maxwellian Ring-Beam function for the pickup ions
function puringbeam(N, vthpara, vthperp, vb, vr, inst2mfa, vdperp1, vdperp2, v3)
    v3mfa = inst2mfa * v3
    vpara = v3mfa[2]
    vperp = sqrt((v3mfa[1] - vdperp1)^2 + (v3mfa[3] - vdperp2)^2)
    a = vr/sqrt(2.0)/vthperp
    return N / (2*π)^1.5 / vthpara / vthperp^2 / (exp(-a^2) + sqrt(π) * a * erfc(-a)) * exp(-0.5 * ((vpara-vb)/vthpara)^2) * exp(-0.5 * ((vperp-vr)/vthperp)^2) * norm(v3)^4 / (2.0*ustrip(u"cm^-3 * km^-3 * s^3", 1u"cm^-2 * s^-1" / (1u"km/s")^4))
end

# read TW data 
TW_datapath2c1hz = "E:/Tianwen-1/MOMAG/2C.V02"
TW_file2c1hz = glob("*$datestr*.dat", TW_datapath2c1hz)
TW_mag2c1hz = reduce(vcat, [identity.(DataFrame(readdlm(file, skipstart=19), :auto)) for file in TW_file2c1hz])
name = ["Time", "Sampling_Rate", "X_MSO", "Y_MSO", "Z_MSO", "Probe_Position_X_MSO", "Probe_Position_Y_MSO", "Probe_Position_Z_MSO", "Roll", "Pitch", "Yaw",  "Quality_Flags"]
rename!(TW_mag2c1hz, name)
TW_mag2c1hz[!, :Time] = map(x->DateTime(x[begin:end-4], DateFormat("y-m-dTH:M:S.s")), TW_mag2c1hz[!, :Time])
TW_mag2c1hz[!, :JulUT] = datetime2julian.(TW_mag2c1hz[!, :Time])
unique!(TW_mag2c1hz) #remove depulicate rows
sort!(TW_mag2c1hz) #sorting
#TW_mag2c1hz[!, [:Probe_Position_X_MSO, :Probe_Position_Y_MSO, :Probe_Position_Z_MSO]] = TW_mag2c1hz[!, [:Probe_Position_X_MSO, :Probe_Position_Y_MSO, :Probe_Position_Z_MSO]] ./R_Mars

TW_magut1hz = TW_mag2c1hz[:, :JulUT].-datejul
TW_BMSO1hz = TW_mag2c1hz[:, [:X_MSO, :Y_MSO, :Z_MSO]]
ind0 = findall(x-> !isnan(x), TW_BMSO1hz[!, 1] ) 
ind1 = findall(x-> isnan(x), TW_BMSO1hz[!, 1]) 
for ib in 1:3 
    interp_linear = linear_interpolation(TW_magut1hz[ind0], TW_BMSO1hz[ind0, ib]; extrapolation_bc=Interpolations.Line())
    TW_BMSO1hz[ind1, ib] = interp_linear(TW_magut1hz[ind1])
end

ns1hz = length(TW_magut1hz)
dt1hz = 1.0
mother = "MORLET"

TW_wave1hz, TW_period1hz, TW_scale1hz, TW_coi1hz = wavelet(reshape(TW_BMSO1hz[:, 1], ns1hz), dt1hz; pad=1, mother=mother)
xpower1hz = abs.(TW_wave1hz).^2
TW_wave1hz, TW_period1hz, TW_scale1hz, TW_coi1hz = wavelet(reshape(TW_BMSO1hz[:, 2], ns1hz), dt1hz; pad=1, mother=mother)
ypower1hz = abs.(TW_wave1hz).^2
TW_wave1hz, TW_period1hz, TW_scale1hz, TW_coi1hz = wavelet(reshape(TW_BMSO1hz[:, 3], ns1hz), dt1hz; pad=1, mother=mother)
zpower1hz = abs.(TW_wave1hz).^2
TW_Bpower1hz =xpower1hz+ypower1hz+zpower1hz

ns1hz = length(TW_magut1hz)
fs = 1
f0 = 0.01
f1 = 0.08
f2 = 0.12
TW_B0MSO1hz = stack([filtfilt(digitalfilter(Lowpass(f0), Butterworth(4); fs=fs), TW_BMSO1hz[:, i]) for i in 1:3])
TW_fcp = dropdims(28.0/1836.0 .* sqrt.(sum(x->x^2, TW_B0MSO1hz, dims=2)); dims=2)
TW_Tcp = 1.0 ./ TW_fcp


ts = DateTime(2022, 06, 12, 9, 51, 00) 
te = DateTime(2022, 06, 12, 9, 55, 00)

#=
ts = DateTime(2022, 06, 12, 9, 59, 00) 
te = DateTime(2022, 06, 12, 10, 01, 00)

ts = DateTime(2022, 06, 12, 10, 10, 00) 
te = DateTime(2022, 06, 12, 10, 15, 00)
=#

ind = ThreadsX.findall(ts .<= TW_mag2c1hz[!, :Time] .<= te)
TW_Bmeanmso = dropdims(mean(Matrix(TW_BMSO1hz[ind, :]), dims=(1,)), dims=(1,))
roll = mean(TW_mag2c1hz[ind, :Roll])|> deg2rad
pitch = mean(TW_mag2c1hz[ind, :Pitch])|> deg2rad
yaw = mean(TW_mag2c1hz[ind, :Yaw])|> deg2rad
TW_mso2scmean = RotXYZ(roll=roll, pitch=pitch, yaw=yaw)
TW_sc2msomean = transpose(TW_mso2scmean)

#=
rot_matrix = zeros(3, 3)
r,p,y = roll,pitch,yaw
cp = cos(p)
sp = sin(p)
cy = cos(y)
sy = sin(y)
cr = cos(r)
sr = sin(r)
rot_matrix[1,:] = [cp*cy, sy*cr+sr*sp*cy, sr*sy-cr*sp*cy]
rot_matrix[2,:] =  [-cp*sy, cr*cy-sr*sp*sy, sr*cy+cr*sp*sy]
rot_matrix[3,:] = [sp, -sr*cp, cr*cp]
=#

#read TW-MINPA data
minpafmod1 = "F:/Work/Data/Tianwen-1/MINPA/2B/2022/".*["HX1-Or_GRAS_MINPA-MOD1-DEF_SCI_N_20220612034833_20220612105300_01288_A.2B"]
massmod1 = [1,2,4,16,28,32,44,64]
energymod1 = [2.81, 3.548928, 4.482167, 5.660813, 7.149401, 9.029433, 11.40385, 14.40264, 18.19001, 22.97332, 29.01447, 36.64422, 46.28032, 58.45036, 73.82067, 93.23282, 117.7497, 148.7135, 187.8198, 237.2095, 299.587, 378.3675, 477.8644, 603.5253, 762.2305, 962.6694, 1215.816, 1535.532, 1939.321, 2449.292, 3093.366, 3906.809, 4934.157, 6231.661, 7870.361, 9939.98, 12553.83, 15855.03, 20024.33, 25290.0]
θmod1 = [11.3, 33.8, 56.2, 78.7]
ϕmod1 = [11.25, 33.75, 56.25, 78.75, 101.25, 123.75, 146.25, 168.75, 191.25, 213.75, 236.25, 258.75, 281.25, 303.75, 326.25, 348.75]

minpamod1 = reduce(vcat, [identity.(DataFrame(readdlm(file), :auto)) for file in minpafmod1])
minpautmod1 = map(x->DateTime(x[begin:end-4], DateFormat("y-m-dTH:M:S.s")), minpamod1[!, 1])
minpajulianmod1 = datetime2julian.(minpautmod1)
minpadfimod1= reshape(Array(minpamod1[:,60:20539]), length(minpautmod1), length(massmod1), length(ϕmod1), length(θmod1), length(energymod1)) #时间 质量数 方位角 俯仰角 能道 #
ind =findall(x->x<0, minpadfimod1)
minpadfimod1[ind] .= 0.0

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

@show "angle between B and Vsw from MINPA", rad2deg(angle(TW_Bmeanminpa, Array(TW_minpavswrough)))

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
    if (rad2deg(lat) <= 0.0) && (TW_minpavelocity[1, 18] ≤ vtot ≤ TW_minpavelocity[1, 31]) && TW_minpaphi[15]-22.5/2.0 <= lon <= TW_minpaphi[11]+22.5/2.0 
    #if (TW_minpavelocity[1, 18] ≤ vtot ≤ TW_minpavelocity[1, 31])
        dir = v3 ./ norm(vtot)
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

minimum(TW_minpa_efluxmean_rect)
#=
TW_minpass = [] #interpolation splines
Snan = SphericalSplines.InterpolationSpline(stack(TW_minpadirections; dims=1), fill(NaN64, 64))
for k in 18:31
    temp = [TW_minpaefluxmid[1, j, i, k] for j in 1:16 for i in 1:4]
    ind = findall(temp .> 0.0)
    if !isempty(ind)
        S = SphericalSplines.InterpolationSpline(stack(TW_minpadirections[ind]; dims=1), log10.(temp[ind]))
        push!(TW_minpass, S)
    else
        push!(TW_minpass, Snan)
    end  
end
TW_minpa_efluxmid_rect = Array{Float64, 3}(undef, length(vperp1), length(vpara), length(vperp2))
@batch for k in 1:length(vperp2)
    for j in 1:length(vpara)
        for i in 1:length(vperp1)
            v3 = Array(vpara[j] .* TW_paraminpa .+ vperp1[i] .* TW_perpminpa1 .+ vperp2[k] .* TW_perpminpa2)
            TW_minpa_efluxmid_rect[i, j, k] = TW_MINPA_eflux_interp(v3, TW_minpass)
            #TW_minpa_eflux_rect[i, j, k] = TW_MINPA_eflux_npinterp(v3)
        end
    end
end

TW_minpass = [] #interpolation splines
Snan = SphericalSplines.InterpolationSpline(stack(TW_minpadirections; dims=1), fill(NaN64, 64))
for k in 18:31
    temp = [TW_minpaefluxlquant[1, j, i, k] for j in 1:16 for i in 1:4]
    ind = findall(temp .> 0.0)
    if !isempty(ind)
        S = SphericalSplines.InterpolationSpline(stack(TW_minpadirections[ind]; dims=1), log10.(temp[ind]))
        push!(TW_minpass, S)
    else
        push!(TW_minpass, Snan)
    end  
end
TW_minpa_efluxlquant_rect = Array{Float64, 3}(undef, length(vperp1), length(vpara), length(vperp2))
@batch for k in 1:length(vperp2)
    for j in 1:length(vpara)
        for i in 1:length(vperp1)
            v3 = Array(vpara[j] .* TW_paraminpa .+ vperp1[i] .* TW_perpminpa1 .+ vperp2[k] .* TW_perpminpa2)
            TW_minpa_efluxlquant_rect[i, j, k] = TW_MINPA_eflux_interp(v3, TW_minpass)
            #TW_minpa_eflux_rect[i, j, k] = TW_MINPA_eflux_npinterp(v3)
        end
    end
end

TW_minpass = [] #interpolation splines
Snan = SphericalSplines.InterpolationSpline(stack(TW_minpadirections; dims=1), fill(NaN64, 64))
for k in 18:31
    temp = [TW_minpaefluxuquant[1, j, i, k] for j in 1:16 for i in 1:4]
    ind = findall(temp .> 0.0)
    if !isempty(ind)
        S = SphericalSplines.InterpolationSpline(stack(TW_minpadirections[ind]; dims=1), log10.(temp[ind]))
        push!(TW_minpass, S)
    else
        push!(TW_minpass, Snan)
    end  
end
=#


function TW_MINPA_eflux_slice(p0, d)
    v = Vector{Float64}()
    v3 = Vector{Vector{Float64}}()
    eflux = Vector{Float64}()
    for k in 18:31 # energy index
        roots = line_sphere_intersection(p0, d, [0.0, 0.0, 0.0], TW_minpavelocity[1, k])
        if !isnothing(roots)
            for root in roots
                push!(v, dot(root, d))
                push!(v3, root)
                sph = SphericalFromCartesian()(root)
                lat = rad2deg(sph.ϕ) 
                lon = rad2deg(sph.θ)
                if rad2deg(lat) <= 0.0 && TW_minpaphi[15]-22.5/2.0 <= lon <= TW_minpaphi[11]+22.5/2.0 
                    S = TW_minpass[k-17]
                    new_direction = Vector(root ./ norm(root))
                    eflux_val = 10.0^S(new_direction)
                    push!(eflux, eflux_val)                   
                else
                    push!(eflux, NaN64)
                end    
            end
        end
    end
    ind = sortperm(v)
    v = v[ind]
    v3 = v3[ind]
    eflux = eflux[ind]
    
    return v, v3, eflux
end


#=
#read TW MINPA moments data
#revised 4pi view
file = matopen("MINPA_FOV4PI_data_20220612.mat")
time = identity.(read(file, "time"))
nH = read(file, "nH")
TH = read(file, "TH")
VH = read(file, "VH")
nO = read(file, "nO")
TO = read(file, "TO")
VO = read(file, "VO")
minpautmod1 = map(x->DateTime(x[begin:end-7], DateFormat("y-m-dTH:M:S.s")), time)
minpajulianmod1 = datetime2julian.(minpautmod1)
=#


#read mvn data
#mvn mag
MVN_datapath2c1hz = "E:/MAVEN/MAG/l2/$year/$month/"
MVN_file2c1hz = glob("*ss1s_*$datestr*.sts", MVN_datapath2c1hz)
MVN_file2c1hz=MVN_file2c1hz[1]
skip = 1 
for line in eachline(MVN_file2c1hz)
    if line[1:6] != "  $year"
        global skip = skip+1
    else
        break
    end
end
MVN_mag2c1hz = identity.(DataFrame(readdlm(MVN_file2c1hz, skipstart=skip), :auto))
name = ["Year", "DOY", "Hour", "Minute", "Second", "Millisecond", "DDOY", "X_MSO", "Y_MSO", "Z_MSO", "Range_MSO", "Probe_Position_X_MSO", "Probe_Position_Y_MSO", "Probe_Position_Z_MSO", "X_PLD", "Y_PLD", "Z_PLD", "Range_PLD"]

rename!(MVN_mag2c1hz, name)
MVN_mag2c1hz[!, :Time] = @. DateTime(MVN_mag2c1hz[:, :Year]) + Day(MVN_mag2c1hz[:, :DOY]-1) + Hour(MVN_mag2c1hz[:, :Hour]) + Minute(MVN_mag2c1hz[:, :Minute]) + Second(MVN_mag2c1hz[:, :Second]) + Millisecond(MVN_mag2c1hz[:, :Millisecond])
MVN_mag2c1hz[!, :JulUT] = datetime2julian.(MVN_mag2c1hz[!, :Time])
MVN_mag2c1hz[!, [:Probe_Position_X_MSO, :Probe_Position_Y_MSO, :Probe_Position_Z_MSO]] = MVN_mag2c1hz[!, [:Probe_Position_X_MSO, :Probe_Position_Y_MSO, :Probe_Position_Z_MSO]].*1u"km" ./R_Mars

MVN_magut1hz = MVN_mag2c1hz[:, :JulUT].-datejul
MVN_BMSO1hz = MVN_mag2c1hz[:, [:X_MSO, :Y_MSO, :Z_MSO]]

ns1hz = length(MVN_magut1hz)
dt1hz = 1.0
mother = "MORLET"
MVN_wave1hz, MVN_period1hz, MVN_scale1hz, MVN_coi1hz = wavelet(reshape(MVN_BMSO1hz[:, 1], ns1hz), dt1hz; pad=1, mother=mother)
xpower1hz = abs.(MVN_wave1hz).^2
MVN_wave1hz, MVN_period1hz, MVN_scale1hz, MVN_coi1hz = wavelet(reshape(MVN_BMSO1hz[:, 2], ns1hz), dt1hz; pad=1, mother=mother)
ypower1hz = abs.(MVN_wave1hz).^2
MVN_wave1hz, MVN_period1hz, MVN_scale1hz, MVN_coi1hz = wavelet(reshape(MVN_BMSO1hz[:, 3], ns1hz), dt1hz; pad=1, mother=mother)
zpower1hz = abs.(MVN_wave1hz).^2
MVN_Bpower1hz =xpower1hz+ypower1hz+zpower1hz


ts = DateTime(2022, 06, 12, 10, 05, 00) 
te = DateTime(2022, 06, 12, 10, 15, 00)
#ts = DateTime(2022, 06, 12, 10, 10, 00) 
#te = DateTime(2022, 06, 12, 10, 15, 00)

#ts = DateTime(2022, 06, 12, 10, 05, 00) 
#te = DateTime(2022, 06, 12, 10, 20, 00)
ind = ThreadsX.findall(ts .<= MVN_mag2c1hz[!, :Time] .<= te)
MVN_Bmeanmso = dropdims(mean(Matrix(MVN_BMSO1hz[ind, :]), dims=(1,)), dims=(1,))

ns1hz = length(MVN_magut1hz)
fs = 1
f0 = 0.01
f1 = 0.08
f2 = 0.12
MVN_B0MSO1hz = stack([filtfilt(digitalfilter(Lowpass(f0), Butterworth(4); fs=fs), MVN_BMSO1hz[:, i]) for i in 1:3])
MVN_fcp = dropdims(28.0/1836.0 .* sqrt.(sum(x->x^2, MVN_B0MSO1hz, dims=2)); dims=2)
MVN_Tcp = 1.0 ./ MVN_fcp

#mvn kp
dstr = Dates.format.(date .+ Day(1) * range(-30, 40), "yyyymmdd")
for (id, ds) in enumerate(dstr)
    MVN_datapathkp = "E:/MAVEN/KP/insitu/$year/$(ds[5:6])/"
    MVN_filekp = glob("*$ds*.tab", MVN_datapathkp)
    MVN_filekp=MVN_filekp[1]
    local skip = 0 
    for line in eachline(MVN_filekp)
        if line[1] == '#'
            skip = skip+1
        else
            break
        end
    end
    
    if id == 1
        global MVN_kp = identity.(DataFrame(readdlm(MVN_filekp, skipstart=skip), :auto))
    else
        global MVN_kp = vcat(MVN_kp, identity.(DataFrame(readdlm(MVN_filekp, skipstart=skip), :auto)))
    end
    @show id, skip, size(MVN_kp)
end
MVN_kp[!, :Time] = DateTime.(MVN_kp[!, 1], DateFormat("y-m-dTH:M:S.s"))
MVN_kp[!, :JulUT] = datetime2julian.(MVN_kp[!, :Time])

#mvn swia
MVN_swiapath = "E:/MAVEN/SWIA/l2/$year/$month/"
MVN_swiafile = glob("*coarsesvy3d_*$datestr*_v02_r01.cdf", MVN_swiapath)
MVN_swia = cdfopen(MVN_swiafile)
MVN_swiatime = MVN_swia["epoch"] |> Array
MVN_swiajulian=datetime2julian.(DateTime.(MVN_swiatime))
MVN_swiaeflux = MVN_swia["diff_en_fluxes"] |> Array
MVN_swiaenergy = MVN_swia["energy_coarse"] |> Array
MVN_swiadenergy = Array(MVN_swia["de_over_e_coarse"])[1] .* MVN_swiaenergy
MVN_swiatheta = MVN_swia["theta_coarse"] |> Array
MVN_swiaphi = MVN_swia["phi_coarse"] |> Array


MVN_swiavelocity = sqrt.(MVN_swiaenergy.*(2.0/ustrip(u"eV", m_p*c_0^2))).*ustrip(u"km/s", c_0)

MVN_swiaind = ThreadsX.findall(ts .<= MVN_swiatime .<= te)
MVN_swiaefluxmean = dropdims(mean(MVN_swiaeflux[:, :, :, MVN_swiaind], dims = (4, ) ), dims=(4, ))
MVN_swiaindmean = round(Int64, mean(MVN_swiaind))
MVN_swiatimemean = MVN_swiatime[MVN_swiaindmean]
MVN_swiajulianmean = MVN_swiajulian[MVN_swiaindmean]

MVN_swiaefluxmid = similar(MVN_swiaefluxmean)
MVN_swiaefluxuquant = similar(MVN_swiaefluxmean)
MVN_swiaefluxlquant = similar(MVN_swiaefluxmean)
for i in axes(MVN_swiaeflux, 1), j in axes(MVN_swiaeflux, 2), k in axes(MVN_swiaeflux, 3)
    swiaeflux = MVN_swiaeflux[i, j, k, MVN_swiaind]
    MVN_swiaefluxmid[i, j, k] = median(swiaeflux)
    MVN_swiaefluxuquant[i, j, k] = quantile(swiaeflux, 0.75)
    MVN_swiaefluxlquant[i, j, k] = quantile(swiaeflux, 0.25)
end

MVN_swia2sc = [0.0 -1.0 0.0; 0.0 0.0 -1.0; 1.0 0.0 0.0] #SWIA instrument to SC transformation matrix
MVN_sc2swia = MVN_swia2sc' #SC to SWIA instrument transformation matrix
ind = argmin(abs.(MVN_kp[!, :Time] .- MVN_swiatimemean))
MVN_sc2msomean = transpose(reshape(Array(MVN_kp[ind, 227:235]), (3, 3))) #SC to MSO transformation matrix at mean time
MVN_swia2msomean = MVN_sc2msomean * MVN_swia2sc #SWIA instrument to MSO transformation matrix at mean time
#QuatRotation(0.5682896622, -0.7372041626, -0.2934748941, 0.2178287607) at meantime unixtime 1655028598  
MVN_Bmeanswia = MVN_swia2msomean' * MVN_Bmeanmso #magnetic field in SWIA coordinates at mean time
MVN_paraswia = MVN_Bmeanswia ./ norm(MVN_Bmeanswia) #unit vector of magnetic field in SWIA coordinates at mean time
i, j =4, 15
ind = argmax(MVN_swiaefluxmean[j, i, :] ./ MVN_swiaenergy[:, 1].^2)
#=
    3D spherical coordinates in the CoordinateTransformations library.

    There are many Spherical coordinate conventions and this library uses a somewhat exotic one. Given a vector v with Cartesian coordinates xyz, let v_xy = [x,y,0] be the orthogonal projection of v on the xy plane.

    r is the radius. It is given by norm(v, 2).
    θ is the azimuth. It is the angle from the x-axis to v_xy
    ϕ is the latitude. It is the angle from v_xy to v.
=#
MVN_swiavswrough = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_swiavelocity[ind, 1], deg2rad(MVN_swiaphi[j, 1]), deg2rad(MVN_swiatheta[i, ind, 1])))
MVN_swiavswrough = Array(MVN_swiavswrough)
@show "angle between B and Vsw from SWIA", rad2deg(angle(MVN_Bmeanswia, MVN_swiavswrough))
MVN_swiadirections = [CartesianFromSpherical()(CoordinateTransformations.Spherical(1.0, deg2rad(MVN_swiaphi[j, 1]), deg2rad(MVN_swiatheta[i, ind, 1]))) for j in 1:16 for i in 1:4] #direction unit vectors
MVN_swiass = [] #interpolation splines
Snan = SphericalSplines.InterpolationSpline(stack(MVN_swiadirections; dims=1), fill(-30.0, 64))
for k in 1:48
    swiaefluxmean = reshape(MVN_swiaefluxmean[:, :, k]', :)
    ind = findall(swiaefluxmean .> 0.0)
    if !isempty(ind)
        S = SphericalSplines.InterpolationSpline(stack(MVN_swiadirections[ind]; dims=1), log10.(swiaefluxmean[ind]))
        push!(MVN_swiass, S)
    else
        push!(MVN_swiass, Snan)
    end
end

#Because of the limited coverage of SWIA， we use the solar wind velocity vector determined by STATIC instrument
#MVN_swiavswrough = 0.5*(MVN_swia2msomean' * MVN_static2msomean * MVN_staticvswrough + MVN_swiavswrough)

MVN_perpswia1 = cross(MVN_paraswia, MVN_swiavswrough) #perpendicular direction1: E=B×Vsw
MVN_perpswia1 = MVN_perpswia1 ./ norm(MVN_perpswia1)
MVN_perpswia2 = cross(MVN_perpswia1, MVN_paraswia) #perpendicular direction 2: N=E×B
MVN_swia2mfa = transpose(hcat(MVN_perpswia1, MVN_paraswia, MVN_perpswia2))

"""
    line_sphere_intersection(p0, d, center, radius)

Compute intersection points between an infinite line and a sphere.

# Arguments
- `p0::Vec`: A point on the line (3D)
- `d::Vec`: Direction vector of the line (3D, non-zero)
- `center::Vec`: Sphere center (3D)
- `radius::Quantity`: Sphere radius (with units)

# Returns
- `nothing` if no intersection
- `(Point,)` if tangent (one point)
- `(Point, Point)` if two intersections
"""
function line_sphere_intersection(p0, d, center, radius)
    oc = p0 - center
    a = dot(d, d)
    b = 2 * dot(d, oc)
    c = dot(oc, oc) - radius^2

    Δ = b^2 - 4 * a * c

    if Δ < 0
        return nothing
    elseif isapprox(Δ, 0; atol=1e-12)
        t = -b / (2 * a)
        return (p0 + t * d,)
    else
        sqrtΔ = sqrt(Δ)
        t1 = (-b - sqrtΔ) / (2 * a)
        t2 = (-b + sqrtΔ) / (2 * a)
        return (p0 + t1 * d, p0 + t2 * d)
    end
end

function MVN_SWIA_eflux_slice(p0, d)
    v = Vector{Float64}()
    v3 = Vector{Vector{Float64}}()
    eflux = Vector{Float64}()
    for k in 13:48 # energy index
        roots = line_sphere_intersection(p0, d, [0.0, 0.0, 0.0], MVN_swiavelocity[k, 1])
        if !isnothing(roots)
            for root in roots
                push!(v, dot(root, d))
                push!(v3, root)
                sph = SphericalFromCartesian()(root)
                lat = rad2deg(sph.ϕ) 
                if abs(lat) <= 45.0
                    S = MVN_swiass[k]
                    new_direction = Vector(root ./ norm(root))
                    eflux_val = 10.0^S(new_direction)
                    push!(eflux, eflux_val)                   
                else
                    push!(eflux, NaN64)
                end    
            end
        end
    end
    ind = sortperm(v)
    v = v[ind]
    v3 = v3[ind]
    eflux = eflux[ind]
    
    return v, v3, eflux
end

function MVN_SWIA_eflux_npslice(p0, d)
    v = Vector{Float64}()
    v3 = Vector{Vector{Float64}}()
    eflux = Vector{Float64}()
    for k in 13:48 # energy index
        roots = line_sphere_intersection(p0, d, [0.0, 0.0, 0.0], MVN_swiavelocity[k, 1])
        if !isnothing(roots)
            for root in roots
                push!(v, dot(root, d))
                push!(v3, root)
                sph = SphericalFromCartesian()(root)
                lat = rad2deg(sph.ϕ) 
                if abs(lat) <= 45.0
                    swiaefluxmean = reshape(MVN_swiaefluxmean[:, :, k]', :)
                    ind = argmin([angle(Array(x), Array(root)) for x in Array(MVN_swiadirections)])  #find the closest direction
                    push!(eflux, swiaefluxmean[ind])                  
                else
                    push!(eflux, NaN64)
                end    
            end
        end
    end
    ind = sortperm(v)
    v = v[ind]
    v3 = v3[ind]
    eflux = eflux[ind]
    
    return v, v3, eflux
end
#x->E, y->B, z->N coordinates

MVN_swiavtot = MVN_swiavelocity[13:48, 1]
agls = range(0, 180, 91)
MVN_swiavarphi = deg2rad.(range(0, 360, 181))
MVN_swiavarphi = MVN_swiavarphi[begin:end-1]
MVN_swiaefluxyz = Matrix{Float64}(undef, 180, 36)
for (ia, agl) in enumerate(agls[begin:end-1])
    d = MVN_paraswia .* cosd(agl) .+ MVN_perpswia2 .* sind(agl)
    #v, v3, eflx = MVN_SWIA_eflux_slice([0.0, 0.0, 0.0], d)
    v, v3, eflx = MVN_SWIA_eflux_slice([0.0, 0.0, 0.0], d)
    MVN_swiaefluxyz[ia, :] = reverse(eflx[begin+36:end])
    MVN_swiaefluxyz[90+ia, :] = eflx[begin:begin+35]
end

MVN_swiaefluxxy = Matrix{Float64}(undef, 180, 36)
for (ia, agl) in enumerate(agls[begin:end-1])
    d = MVN_perpswia1 .* cosd(agl) .+ MVN_paraswia .* sind(agl)
    #v, v3, eflx = MVN_SWIA_eflux_slice([0.0, 0.0, 0.0], d)
    v, v3, eflx = MVN_SWIA_eflux_slice([0.0, 0.0, 0.0], d)
    MVN_swiaefluxxy[ia, :] = reverse(eflx[begin+36:end])
    MVN_swiaefluxxy[90+ia, :] = eflx[begin:begin+35]
end

MVN_swiaefluxxz = Matrix{Float64}(undef, 180, 36)
for (ia, agl) in enumerate(agls[begin:end-1])
    d = MVN_perpswia1 .* cosd(agl) .+ MVN_perpswia2 .* sind(agl)
    #v, v3, eflx = MVN_SWIA_eflux_slice([0.0, 0.0, 0.0], d)
    v, v3, eflx = MVN_SWIA_eflux_slice([0.0, 0.0, 0.0], d)
    MVN_swiaefluxxz[ia, :] = reverse(eflx[begin+36:end])
    MVN_swiaefluxxz[90+ia, :] = eflx[begin:begin+35]
end

function MVN_SWIA_eflux_interp(v3)
    vs = SphericalFromCartesian()(v3) 
    vtot = vs.r
    lon = vs.θ
    lat = vs.ϕ
    if (abs(rad2deg(lat)) <= 45.0) & (MVN_swiavelocity[13, 1] ≥ vtot ≥ MVN_swiavelocity[48, 1])
        dir = v3 ./ norm(vtot)
        logv = log10.(MVN_swiavelocity[13:48, 1])
        logeflux = Vector{Float64}(undef, 36)
        for k in 13:48 # energy index
            S = MVN_swiass[k]
            logeflux[k-12] = S(dir)
        end
        x = reverse(logv[begin:2:end])
        y = reverse(logeflux[begin:2:end])
        return 10.0^interpolate_weno4([log10(vtot)], x, y)[1]
    else
        return NaN64
    end

    #=
    if (abs(rad2deg(lat)) <= 45.0) & (MVN_swiavelocity[13, 1] ≥ vtot ≥ MVN_swiavelocity[48, 1])
        dir = v3 ./ norm(vtot)
        logv = log10.(MVN_swiavelocity[13:48, 1])
        logeflux = Vector{Float64}(undef, 36)
        for k in 13:48 # energy index
            S = MVN_swiass[k]
            logeflux[k-12] = S(dir)
        end
        x = reverse(logv)
        y = reverse(logeflux)
        return 10.0^interpolate_weno4([log10(vtot)], x, y)[1]
    else
        return NaN64
    end
    =#
end

function MVN_SWIA_eflux_npinterp(v3)
    vs = SphericalFromCartesian()(v3) 
    vtot = vs.r
    lon = vs.θ
    lat = vs.ϕ

    vtot = norm(v3)
    swialon = MVN_swiaphi[:, 1] .- 22.5/2.0  #-180 to 180
    ind = findall(x -> x > 180.0, swialon)
    swialon[ind] = swialon[ind] .- 360.0
    swialat = MVN_swiatheta[:, 13, 1] .- 22.5/2.0 #-45 to 45
    append!(swialat, 45.0) # add the last latitude bin edge
    swiav = MVN_swiavelocity[:, 1] # from high to low
    if (abs(rad2deg(lat)) <= 45.0) & (MVN_swiavelocity[13, 1] ≥ vtot ≥ MVN_swiavelocity[48, 1])
        indlon = findfirst(x -> x ≥ rad2deg(lon), swialon) - 1
        if indlon == 0
            indlon = 16
        end
        indlat = findfirst(x -> x ≥ rad2deg(lat), swialat) -1
        indv = findfirst(x -> x ≤ vtot, swiav) -1

        #@show lon, lat, vtot, indlon, indlat
        return MVN_swiaefluxmean[indlon, indlat, indv]
    else
        return NaN64
    end
end

vpara = collect(range(-800.0, 800.0, step = 10.0))
vperp1 = collect(range(-800.0, 800.0, step = 10.0))
vperp2 = collect(range(-800.0, 800.0, step = 10.0))
MVN_swia_eflux_rect = Array{Float64, 3}(undef, length(vperp1), length(vpara), length(vperp2))
@batch for k in 1:length(vperp2)
    for j in 1:length(vpara)
        for i in 1:length(vperp1)
            v3 = Array(vpara[j] .* MVN_paraswia .+ vperp1[i] .* MVN_perpswia1 .+ vperp2[k] .* MVN_perpswia2)
            MVN_swia_eflux_rect[i, j, k] = MVN_SWIA_eflux_interp(v3)
            #MVN_swia_eflux_rect[i, j, k] = MVN_SWIA_eflux_npinterp(v3)
        end
    end
end



function density_from_swiaflux(eflux)
    #total density
    Ntot = 0.0
    for k in 1:48
        for i in 1:4
            for j in 1:16
                Ntot += eflux[j, i, k] / MVN_swiavelocity[k, 1] / MVN_swiaenergy[k, 1] * MVN_swiadenergy[k, 1] * cosd(MVN_swiatheta[i, k, 1]) * deg2rad(MVN_swiatheta[2, k, 1]-MVN_swiatheta[1, k, 1]) * deg2rad(22.5)
            end
        end
    end
    Ntot *= ustrip(u"cm^-3", 1u"cm^-2*s^-1" / 1u"km/s")

    #solar wind core
    N = 0.0
    for i in 3:4
        for j in 13:16 
            for k in 26:48 
                N += eflux[j, i, k] / MVN_swiavelocity[k, 1] / MVN_swiaenergy[k, 1] * MVN_swiadenergy[k, 1] * cosd(MVN_swiatheta[i, k, 1]) * deg2rad(MVN_swiatheta[2, k, 1]-MVN_swiatheta[1, k, 1]) * deg2rad(22.5)
            end
        end
    end

    for i in 1:2
        for j in 13:16 #1:16
            for k in 29:48
                N += eflux[j, i, k] / MVN_swiavelocity[k, 1] / MVN_swiaenergy[k, 1] * MVN_swiadenergy[k, 1] * cosd(MVN_swiatheta[i, k, 1]) * deg2rad(MVN_swiatheta[2, k, 1]-MVN_swiatheta[1, k, 1]) * deg2rad(22.5)
            end
        end
    end
    N *= ustrip(u"cm^-3", 1u"cm^-2*s^-1" / 1u"km/s")
    Nsw = N

    #pickup ions
    N = 0.0
    for i in 1:1
        for j in 11:16 #1:16
            for k in 1:28 #7:32
                N += eflux[j, i, k] / MVN_swiavelocity[k, 1] / MVN_swiaenergy[k, 1] * MVN_swiadenergy[k, 1] * cosd(MVN_swiatheta[i, k, 1]) * deg2rad(MVN_swiatheta[2, k, 1]-MVN_swiatheta[1, k, 1]) * deg2rad(22.5)
            end
        end
    end

    for i in 2:2
        for j in 9:16 #1:16
            for k in 1:28 #7:32
                N += eflux[j, i, k] / MVN_swiavelocity[k, 1] / MVN_swiaenergy[k, 1] * MVN_swiadenergy[k, 1] * cosd(MVN_swiatheta[i, k, 1]) * deg2rad(MVN_swiatheta[2, k, 1]-MVN_swiatheta[1, k, 1]) * deg2rad(22.5)
            end
        end
    end

    for i in 2:2
        for j in 1:1 #1:16
            for k in 1:28 #7:32
                N += eflux[j, i, k] / MVN_swiavelocity[k, 1] / MVN_swiaenergy[k, 1] * MVN_swiadenergy[k, 1] * cosd(MVN_swiatheta[i, k, 1]) * deg2rad(MVN_swiatheta[2, k, 1]-MVN_swiatheta[1, k, 1]) * deg2rad(22.5)
            end
        end
    end

    for i in 3:4
        for j in 13:16 #1:16
            for k in 1:25 #7:32
                N += eflux[j, i, k] / MVN_swiavelocity[k, 1] / MVN_swiaenergy[k, 1] * MVN_swiadenergy[k, 1] * cosd(MVN_swiatheta[i, k, 1]) * deg2rad(MVN_swiatheta[2, k, 1]-MVN_swiatheta[1, k, 1]) * deg2rad(22.5)
            end
        end
    end

    for i in 3:4
        for j in 1:12 #1:16
            for k in 1:29 #7:32
                N += eflux[j, i, k] / MVN_swiavelocity[k, 1] / MVN_swiaenergy[k, 1] * MVN_swiadenergy[k, 1] * cosd(MVN_swiatheta[i, k, 1]) * deg2rad(MVN_swiatheta[2, k, 1]-MVN_swiatheta[1, k, 1]) * deg2rad(22.5)
            end
        end
    end
    N *= ustrip(u"cm^-3", 1u"cm^-2*s^-1" / 1u"km/s")
    Npu = N

    return Ntot, Nsw, Npu
end




#mvn static
MVN_staticpath = "E:/MAVEN/STATIC/l2/$year/$month/"
MVN_staticfile = glob("*32e4d16a8m_*$datestr*_v02_r07.cdf", MVN_staticpath)
MVN_static = cdfopen(MVN_staticfile)
MVN_staticswp = MVN_static["swp_ind"] |> Array 
MVN_staticswp .+= 1  # +1 because Julia is 1-based indexing while the original data is 0-based indexing
#Index that identifies the energy and deflector sweep look up tables (LUT) for the sensor. SWP_IND is an index that selects the following support data arrays: ENERGY, DENERGY, THETA, DTHETA, PHI, DPHI, DOMEGA, GF and MASS_ARR. (NUM_DISTS elements), SWP_IND Le NSWP
MVN_statictime = MVN_static["epoch"] |> Array
MVN_staticjulian=datetime2julian.(DateTime.(MVN_statictime))
MVN_staticeflux = MVN_static["eflux"] |> Array
MVN_staticmass = MVN_static["mass_arr"] |> Array
MVN_staticenergy = MVN_static["energy"] |> Array
MVN_statictheta = MVN_static["theta"] |> Array #-45 to 45 deg, corresponding to latitude in the STATIC/APP coordinate system
MVN_staticphi = MVN_static["phi"] |> Array #-180 to 180 deg, corresponding to longitude in the STATIC/APP coordinate system
MVN_staticdenergy = MVN_static["denergy"] |> Array
MVN_staticdtheta = MVN_static["dtheta"] |> Array
MVN_staticdphi = MVN_static["dphi"] |> Array
MVN_staticdomega = MVN_static["domega"] |> Array
MVN_staticvelocity = sqrt.(MVN_staticenergy.*(2.0/ustrip(u"eV", m_p*c_0^2))./MVN_staticmass).*ustrip(u"km/s", c_0)
MVN_staticpos_sc_mso = MVN_static["pos_sc_mso"] |> Array
MVN_staticpos_sc_mso /= ustrip(u"km", R_Mars)

# Doc about quaternion
# https://danceswithcode.net/engineeringnotes/quaternions/quaternions.html
MVN_staticquat_mso = MVN_static["quat_mso"] |> Array #Quaternion elements to rotate from STATIC coordinates (same as APP coordinates) to MSO coordinates (NUM_DISTS, 3)
MVN_staticquat_sc = MVN_static["quat_sc"] |> Array #Quaternion elements to rotate from STATIC coordinates (same as APP coordinates) to SC coordinates (NUM_DISTS, 3)
MVN_staticmsotrans = [QuatRotation(quat) for quat in eachcol(MVN_staticquat_mso)] # active (right handed) rotations where applicable.

MVN_staticind = ThreadsX.findall(ts .<= MVN_statictime .<= te)
@show minimum(MVN_staticswp[MVN_staticind]), maximum(MVN_staticswp[MVN_staticind])
MVN_staticswpind = MVN_staticswp[MVN_staticind[1]]
MVN_staticefluxmean = dropdims(mean(MVN_staticeflux[:, :, :, MVN_staticind], dims = (4, ) ), dims=(4, ))
MVN_staticpsdmean = ustrip(u"cm^-3 * km^-3 * s^3", 0.5 * 1u"cm^-2 * s^-1" * (m_p/1u"eV")^2)  .* MVN_staticefluxmean .* (MVN_staticmass[:,:,:,MVN_staticswpind,1]./ MVN_staticenergy[:,:,:,MVN_staticswpind,1]).^2
size(MVN_staticefluxmean)
size(MVN_staticmass[:,:,:,MVN_staticswpind,1])
size(MVN_staticenergy[:,:,:,MVN_staticswpind,1])
MVN_staticefluxmid = similar(MVN_staticefluxmean)
MVN_staticefluxuquant = similar(MVN_staticefluxmean)
MVN_staticefluxlquant = similar(MVN_staticefluxmean)
for i in axes(MVN_staticeflux, 1), j in axes(MVN_staticeflux, 2), k in axes(MVN_staticeflux, 3)
    data = MVN_staticeflux[i, j, k, MVN_staticind]
    MVN_staticefluxmid[i, j, k] = quantile(data, 0.5)
    MVN_staticefluxuquant[i, j, k] = quantile(data, 0.75)
    MVN_staticefluxlquant[i, j, k] = quantile(data, 0.25)
end
# lower and upper bounds of energy channels during the selected time interval
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
MVN_staticindmean = round(Int64, mean(MVN_staticind))
MVN_statictimemean = MVN_statictime[MVN_staticindmean]
MVN_staticjulianmean = MVN_staticjulian[MVN_staticindmean]
MVN_static2msomean = QuatRotation(MVN_staticquat_mso[:, MVN_staticindmean]) # STATIC instrument to MSO transformation matrix at mean time
MVN_Bmeanstatic = MVN_static2msomean' * MVN_Bmeanmso #magnetic field in STATIC coordinates at mean time
MVN_parastatic = MVN_Bmeanstatic ./ norm(MVN_Bmeanstatic) #unit vector of magnetic field in STATIC coordinates at mean time
i, j =2, 14
ind = argmax(MVN_staticefluxmean[1, (j-1)*4+i, :])
ind = argmax(MVN_staticefluxmean[1, (j-1)*4+i, :] ./ MVN_staticenergy[1, 1, :, MVN_staticswpind, 1].^2)
#=
    3D spherical coordinates in the CoordinateTransformations library.

    There are many Spherical coordinate conventions and this library uses a somewhat exotic one. Given a vector v with Cartesian coordinates xyz, let v_xy = [x,y,0] be the orthogonal projection of v on the xy plane.

    r is the radius. It is given by norm(v, 2).
    θ is the azimuth. It is the angle from the x-axis to v_xy
    ϕ is the latitude. It is the angle from v_xy to v.
=#
MVN_staticvswrough = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, ind, MVN_staticswpind, 1], deg2rad(MVN_staticphi[1, (j-1)*4+i, ind, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, (j-1)*4+i, ind, MVN_staticswpind, 1])))
@show "angle between Vsw from SWIA and STATIC", rad2deg(angle(MVN_swia2msomean * MVN_swiavswrough, Array(MVN_static2msomean * MVN_staticvswrough)) )

#According to the symmetry of elfux profiles, the SWIA velocity appears to be more reasoanable than STATIC velocity
#MVN_staticvswrough = MVN_static2msomean' * MVN_swia2msomean * MVN_swiavswrough

#MVN_staticvswrough = MVN_staticvswrough .+ 100 .* MVN_parastatic #adjust the solar wind velocity magnitude to be more reasonable

MVN_staticdirections = [CartesianFromSpherical()(CoordinateTransformations.Spherical(1.0, deg2rad(MVN_staticphi[1, ij, ind, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, ij, ind, MVN_staticswpind, 1]))) for ij in 1:64] #direction unit vectors
MVN_staticss = [] #interpolation splines
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
@show "angle between B and Vsw from STATIC", rad2deg(angle(MVN_Bmeanstatic, MVN_staticvswrough))
MVN_perpstatic1 = cross(MVN_parastatic, MVN_staticvswrough) #perpendicular direction1: E=B×Vsw
MVN_perpstatic1 = MVN_perpstatic1 ./ norm(MVN_perpstatic1)
MVN_perpstatic2 = cross(MVN_perpstatic1, MVN_parastatic) #perpendicular direction 2: N=E×B
MVN_static2mfa = transpose(hcat(MVN_perpstatic1, MVN_parastatic, MVN_perpstatic2))

#MVN_staticphi[1, 1:4:64, :, 24, 1]
#MVN_statictheta[1, 1:4, :, 24, 1]

function MVN_STATIC_eflux_slice(p0, d)
    v = Vector{Float64}()
    v3 = Vector{Vector{Float64}}()
    eflux = Vector{Float64}()
    for k in 6:32 # energy index
        roots = line_sphere_intersection(p0, d, [0.0, 0.0, 0.0], MVN_staticvelocity[1, 1, k, MVN_staticswpind, 1])
        if !isnothing(roots)
            for root in roots
                push!(v, dot(root, d))
                push!(v3, root)
                sph = SphericalFromCartesian()(root)
                lat = rad2deg(sph.ϕ) 
                if abs(lat) <= 45.0
                    S = MVN_staticss[k]
                    new_direction = Vector(root ./ norm(root))
                    eflux_val = 10.0^S(new_direction)
                    push!(eflux, eflux_val)                  
                else
                    push!(eflux, NaN64)
                end    
            end
        end
    end
    ind = sortperm(v)
    v = v[ind]
    v3 = v3[ind]
    eflux = eflux[ind]
    
    return v, v3, eflux
end

function MVN_STATIC_eflux_npslice(p0, d)
    v = Vector{Float64}()
    v3 = Vector{Vector{Float64}}()
    eflux = Vector{Float64}()
    for k in 6:32 # energy index
        roots = line_sphere_intersection(p0, d, [0.0, 0.0, 0.0], MVN_staticvelocity[1, 1, k, MVN_staticswpind, 1])
        if !isnothing(roots)
            for root in roots
                push!(v, dot(root, d))
                push!(v3, root)
                sph = SphericalFromCartesian()(root)
                lat = rad2deg(sph.ϕ) 
                if abs(lat) <= 45.0
                    ind = argmin([angle(Array(x), Array(root)) for x in Array(MVN_staticdirections)])  #find the closest direction
                    push!(eflux, MVN_staticefluxmean[1, ind, k])                   
                else
                    push!(eflux, NaN64)
                end    
            end
        end
    end
    ind = sortperm(v)
    v = v[ind]
    v3 = v3[ind]
    eflux = eflux[ind]
    
    return v, v3, eflux
end

#=
i , j = 2, 14
MVN_staticvtot = MVN_staticvelocity[1, (j-1)*4+i, 6:32, MVN_staticswpind, 1]
agls = range(0, 180, 91)
MVN_staticvarphi = deg2rad.(range(0, 360, 181))
MVN_staticvarphi = MVN_staticvarphi[begin:end-1]
#x->E, y->B, z->N coordinates
MVN_staticefluxyz = Matrix{Float64}(undef, 180, 26)
for (ia, agl) in enumerate(agls[begin:end-1])
    d = MVN_parastatic .* cosd(agl) .+ MVN_perpstatic2 .* sind(agl)
    #v, v3, eflx = MVN_STATIC_eflux_slice([0.0, 0.0, 0.0], d)
    v, v3, eflx = MVN_STATIC_eflux_slice([0.0, 0.0, 0.0], d)
    MVN_staticefluxyz[ia, :] = reverse(eflx[begin+26:end])
    MVN_staticefluxyz[90+ia, :] = eflx[begin:begin+25]
end

MVN_staticefluxxy = Matrix{Float64}(undef, 180, 26)
for (ia, agl) in enumerate(agls[begin:end-1])
    d = MVN_perpstatic1 .* cosd(agl) .+ MVN_parastatic .* sind(agl)
    #v, v3, eflx = MVN_STATIC_eflux_slice([0.0, 0.0, 0.0], d)
    v, v3, eflx = MVN_STATIC_eflux_slice([0.0, 0.0, 0.0], d)
    MVN_staticefluxxy[ia, :] = reverse(eflx[begin+26:end])
    MVN_staticefluxxy[90+ia, :] = eflx[begin:begin+25]
end

MVN_staticefluxxz = Matrix{Float64}(undef, 180, 26)
for (ia, agl) in enumerate(agls[begin:end-1])
    d = MVN_perpstatic1 .* cosd(agl) .+ MVN_perpstatic2 .* sind(agl)
    #v, v3, eflx = MVN_STATIC_eflux_slice([0.0, 0.0, 0.0], d)
    v, v3, eflx = MVN_STATIC_eflux_slice([0.0, 0.0, 0.0], d)
    MVN_staticefluxxz[ia, :] = reverse(eflx[begin+26:end])
    MVN_staticefluxxz[90+ia, :] = eflx[begin:begin+25]
end
=#

function MVN_STATIC_eflux_interp(v3)
    vs = SphericalFromCartesian()(v3) 
    vtot = vs.r
    lon = vs.θ
    lat = vs.ϕ
    if (abs(rad2deg(lat)) <= 45.0) & (MVN_staticvelocity[1, 1, 6, MVN_staticswpind, 1] ≥ vtot ≥ MVN_staticvelocity[1, 1, 32, MVN_staticswpind, 1])
        dir = v3 ./ norm(vtot)
        logv = log10.(MVN_staticvelocity[1, 1, 6:32, MVN_staticswpind, 1])
        logeflux = Vector{Float64}(undef, 27)
        for k in 6:32 # energy index
            S = MVN_staticss[k]
            logeflux[k-5] = S(dir)
        end
        x = reverse(logv)
        y = reverse(logeflux)
        return 10.0^interpolate_weno4([log10(vtot)], x, y)[1]
    else
        return NaN64
    end
end

function MVN_STATIC_eflux_npinterp(v3)
    vs = SphericalFromCartesian()(v3) 
    vtot = vs.r
    lon = vs.θ
    lat = vs.ϕ
    vtot = norm(v3)
    staticlon = MVN_staticphigridH  #-180 to 180
    staticlat = 90.0 .- MVN_staticthetagridH #45 to -45
    staticv = reverse(MVN_staticvelocitygridH) # from high to low
    if (abs(rad2deg(lat)) <= 45.866) & (MVN_staticvelocity[1, 1, 7, MVN_staticswpind, 1] ≥ vtot ≥ MVN_staticvelocity[1, 1, 32, MVN_staticswpind, 1])
        indlon = findfirst(x -> x ≥ rad2deg(lon), staticlon) - 1
        if indlon == 0
            indlon = 16
        end
        indlat = findfirst(x -> x ≤ rad2deg(lat), staticlat) -1
        indv = findfirst(x -> x ≤ vtot, staticv) -1

        #@show lon, lat, vtot, indlon, indlat
        return MVN_staticefluxmean[1, (indlon-1)*4+indlat, indv]
    else
        return NaN64
    end
end

#=
i=2
j=13
MVN_staticefluxmean[1, (j-1)*4+i, 9] = NaN64
i=1
j=13
MVN_staticefluxmean[1, (j-1)*4+i, 9] = NaN64
i=3
j=13
MVN_staticefluxmean[1, (j-1)*4+i, 9] = NaN64
i=3
j=14
MVN_staticefluxmean[1, (j-1)*4+i, 9] = NaN64

i=2
for j in 11:16
    MVN_staticefluxmean[1, (j-1)*4+i, 10] = NaN64
end
=#
vpara = collect(range(-800.0, 800.0, step = 10.0))
vperp1 = collect(range(-800.0, 800.0, step = 10.0))
vperp2 = collect(range(-800.0, 800.0, step = 10.0))
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

#density 
function density_from_staticeflux(eflux)
    #total density
    Ntot = sum(eflux[1, :, :] ./ MVN_staticvelocity[1, :, :, MVN_staticswpind, 1] ./ MVN_staticenergy[1, :, :, MVN_staticswpind, 1] .* MVN_staticdenergy[1, :, :, MVN_staticswpind, 1] .* MVN_staticdomega[1, :, :, MVN_staticswpind, 1]) * ustrip(u"cm^-3", 1u"cm^-2*s^-1" / 1u"km/s")

    #solar wind core
    N = 0.0
    for i in 1:3
        for j in 13:16 #1:16
            for k in 12:32 #7:32
                N += eflux[1, (j-1)*4+i, k] / MVN_staticvelocity[1, (j-1)*4+i, k, MVN_staticswpind, 1] / MVN_staticenergy[1, (j-1)*4+i, k, MVN_staticswpind, 1] * MVN_staticdenergy[1, (j-1)*4+i, k, MVN_staticswpind, 1] * MVN_staticdomega[1, (j-1)*4+i, k, MVN_staticswpind, 1]
            end
        end
    end

    for i in 4:4
        for j in 12:13 #1:16
            for k in 12:32 #7:32
                N += eflux[1, (j-1)*4+i, k] / MVN_staticvelocity[1, (j-1)*4+i, k, MVN_staticswpind, 1] / MVN_staticenergy[1, (j-1)*4+i, k, MVN_staticswpind, 1] * MVN_staticdenergy[1, (j-1)*4+i, k, MVN_staticswpind, 1] * MVN_staticdomega[1, (j-1)*4+i, k, MVN_staticswpind, 1]
            end
        end
    end
    N *= ustrip(u"cm^-3", 1u"cm^-2*s^-1" / 1u"km/s")
    Nsw = N

    #pickup ions
    N = 0.0
    for i in 1:3
        for j in 13:16 #1:16
            for k in 1:11 #7:32
                N += eflux[1, (j-1)*4+i, k] / MVN_staticvelocity[1, (j-1)*4+i, k, MVN_staticswpind, 1] / MVN_staticenergy[1, (j-1)*4+i, k, MVN_staticswpind, 1] * MVN_staticdenergy[1, (j-1)*4+i, k, MVN_staticswpind, 1] * MVN_staticdomega[1, (j-1)*4+i, k, MVN_staticswpind, 1]
            end
        end
    end

    for i in 4:4
        for j in 14:16 #1:16
            for k in 1:13 #7:32
                N += eflux[1, (j-1)*4+i, k] / MVN_staticvelocity[1, (j-1)*4+i, k, MVN_staticswpind, 1] / MVN_staticenergy[1, (j-1)*4+i, k, MVN_staticswpind, 1] * MVN_staticdenergy[1, (j-1)*4+i, k, MVN_staticswpind, 1] * MVN_staticdomega[1, (j-1)*4+i, k, MVN_staticswpind, 1]
            end
        end
    end

    for i in 1:4
        for j in 1:3 #1:16
            for k in 1:13 #7:32
                N += eflux[1, (j-1)*4+i, k] / MVN_staticvelocity[1, (j-1)*4+i, k, MVN_staticswpind, 1] / MVN_staticenergy[1, (j-1)*4+i, k, MVN_staticswpind, 1] * MVN_staticdenergy[1, (j-1)*4+i, k, MVN_staticswpind, 1] * MVN_staticdomega[1, (j-1)*4+i, k, MVN_staticswpind, 1]
            end
        end
    end

    for i in 1:4
        for j in 7:12 #1:16
            for k in 1:14 #7:32
                N += eflux[1, (j-1)*4+i, k] / MVN_staticvelocity[1, (j-1)*4+i, k, MVN_staticswpind, 1] / MVN_staticenergy[1, (j-1)*4+i, k, MVN_staticswpind, 1] * MVN_staticdenergy[1, (j-1)*4+i, k, MVN_staticswpind, 1] * MVN_staticdomega[1, (j-1)*4+i, k, MVN_staticswpind, 1]
            end
        end
    end
    N *= ustrip(u"cm^-3", 1u"cm^-2*s^-1" / 1u"km/s")
    Npu = N

    return Ntot, Nsw, Npu
end

Ntot, Nsw, Npu = density_from_staticeflux(MVN_staticefluxmean)
Ntot, Nsw, Npu = density_from_staticeflux(MVN_staticefluxmid)
Ntot, Nsw, Npu = density_from_staticeflux(MVN_staticefluxlquant)
Ntot, Nsw, Npu = density_from_staticeflux(MVN_staticefluxuquant)


# solar wind core fitting
v_data = []
f_data = []
vbulk = MVN_staticvswrough
mag = MVN_Bmeanstatic
#solar wind core
for i in 1:2
    for j in 13:15 #1:16
        for k in 12:19 #7:32
            push!(v_data, Array(CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, k, MVN_staticswpind, 1], deg2rad(MVN_staticphi[1, (j-1)*4+i, k, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, (j-1)*4+i, k, MVN_staticswpind, 1])))))
            push!(f_data, log(MVN_staticpsdmean[1, (j-1)*4+i, k]))
        end
    end
end

for i in 4:4
    for j in 12:13 #1:16
        for k in 12:19 #7:32
            push!(v_data, Array(CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, k, MVN_staticswpind, 1], deg2rad(MVN_staticphi[1, (j-1)*4+i, k, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, (j-1)*4+i, k, MVN_staticswpind, 1])))))
            push!(f_data, log(MVN_staticpsdmean[1, (j-1)*4+i, k]))
        end
    end
end
initial_guess = [log(5.0), ustrip(u"km/s", sqrt(10.9 * 1u"eV" / m_p)), ustrip(u"km/s", sqrt(27.7 * 1u"eV" / m_p)), MVN_staticvswrough...] 
lower_bounds = [log(3.0), ustrip(u"km/s", sqrt(9.0 * 1u"eV" / m_p)), ustrip(u"km/s", sqrt(22.0 * 1u"eV" / m_p)), -300, -300.0, -300.0]  
upper_bounds = [log(25.0), ustrip(u"km/s", sqrt(15.0 * 1u"eV" / m_p)), ustrip(u"km/s", sqrt(32.0 * 1u"eV" / m_p)), 300.0, 300, 300.0] 

#=
eperp1, epara, eperp2 = field_aligned_coord(mag; vbulk = vbulk)
initial_guess = [log(5.0), ustrip(u"km/s", sqrt(100.0 * 1u"eV" / m_p)), ustrip(u"km/s", sqrt(100.0 * 1u"eV" / m_p)), dot(MVN_staticvswrough, epara), dot(MVN_staticvswrough, eperp1), dot(MVN_staticvswrough, eperp2)] 
lower_bounds = [log(0.1), ustrip(u"km/s", sqrt(1.0 * 1u"eV" / m_p)), ustrip(u"km/s", sqrt(1.0 * 1u"eV" / m_p)), 100, -100.0, 0.0]  
upper_bounds = [log(20.0), ustrip(u"km/s", sqrt(300.0 * 1u"eV" / m_p)), ustrip(u"km/s", sqrt(300.0 * 1u"eV" / m_p)), 300.0, 100, 200.0] 
=#


#=
MVN_staticNsw, MVN_staticTparasw, MVN_staticTperpsw, MVN_staticvdsw, result = least_squares_fit(v_data, f_data, mag, initial_guess,lower_bounds, upper_bounds; vbulk = vbulk) #initial guess for the fitting parameters: density, temperature, and bulk velocity

MVN_staticNsw, MVN_staticTparasw, MVN_staticTperpsw, MVN_staticvdsw, result = least_squares_fit_bb(v_data, f_data, mag, initial_guess,lower_bounds, upper_bounds; vbulk = vbulk, method = :pure_bbo)

eperp1, epara, eperp2 = field_aligned_coord(mag; vbulk = vbulk)
MVN_staticparams_fit = BiMaxwellianParams(log(MVN_staticNsw), ustrip(u"km/s", sqrt(MVN_staticTparasw * 1u"eV" / m_p)), ustrip(u"km/s", sqrt(MVN_staticTperpsw * 1u"eV" / m_p)), dot(MVN_staticvdsw, epara), dot(MVN_staticvdsw, eperp1), dot(MVN_staticvdsw, eperp2))
exp(log_bi_maxwellian_density([dot(MVN_staticvswrough, eperp1), dot(MVN_staticvswrough, epara), dot(MVN_staticvswrough, eperp2)], MVN_staticparams_fit))

norm(MVN_staticvswrough), norm(MVN_staticvdsw)
norm(MVN_swiavswrough)

i, j, k = 1, 14, 14
MVN_staticpsdmean[1, (j-1)*4+i, k]
vec = CartesianFromSpherical()(CoordinateTransformations.Spherical(MVN_staticvelocity[1, (j-1)*4+i, k, MVN_staticswpind, 1], deg2rad(MVN_staticphi[1, (j-1)*4+i, k, MVN_staticswpind, 1]), deg2rad(MVN_statictheta[1, (j-1)*4+i, k, MVN_staticswpind, 1])))
norm(vec)

rad2deg(angle(vec, MVN_Bmeanstatic))

=#

data = load("sta_H+fitting_20220612.jld2")
MVN_staticHfitting = Array{Any, 4}(undef, 4, 4, 16, 2)
fill!(MVN_staticHfitting, 0.0)
for ϕ in 1:16
    for θ in 1:4
        fitting = data["fit_parameters"]
        if haskey(fitting, (ϕ, θ))
            MVN_staticHfitting[1, θ, ϕ, 1] = fitting[(ϕ, θ)].pickup_part.N * 1u"m^-3" #density parameter in the Maxwellian distribution
            MVN_staticHfitting[2, θ, ϕ, 1] = fitting[(ϕ, θ)].pickup_part.T * 1u"K" #temperature
            MVN_staticHfitting[3, θ, ϕ, 1] = fitting[(ϕ, θ)].pickup_part.v0 *1u"m/s" #drift velocity
            MVN_staticHfitting[4, θ, ϕ, 1] = fitting[(ϕ, θ)].pickup_part.density *1u"m^-3" #total number of particles

            MVN_staticHfitting[1, θ, ϕ, 2] = fitting[(ϕ, θ)].sw_part.N * 1u"m^-3" #density parameter in the Maxwellian distribution
            MVN_staticHfitting[2, θ, ϕ, 2] = fitting[(ϕ, θ)].sw_part.T * 1u"K" #temperature
            MVN_staticHfitting[3, θ, ϕ, 2] = fitting[(ϕ, θ)].sw_part.v0 *1u"m/s" #drift velocity
            MVN_staticHfitting[4, θ, ϕ, 2] = fitting[(ϕ, θ)].sw_part.density *1u"m^-3" #total number of particles
        end
    end
end





data["fit_parameters"][(11, 1)].pickup_part.N * 1u"m^-3"

MVN_staticHfitting[1, 1, 16, 2] 



#mvn NGIMS
ngipath = "E:/MAVEN/NGIMS/l3/2022/" .* ["04", "05","06", "07"]
ngifile = reduce(vcat, [glob("mvn_ngi_l3_res*.csv", path) for path in ngipath])
ngi = reduce(vcat, [identity.(DataFrame(CSV.File(file))) for file in ngifile])
ngi_julian = datetime2julian.(ngi[:, :t_utc])

#solar wind time IntervalsBetween
swtime = identity.(DataFrame(readdlm("solarwind.txt", skipstart=0), :auto))
swtime[!, 1] = DateTime.(swtime[!, 1], DateFormat("y-m-dTH:M:S"))
swtime[!, 2] = DateTime.(swtime[!, 2], DateFormat("y-m-dTH:M:S"))




