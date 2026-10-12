using Glob
using CommonDataFormat
using Dates
using Statistics
using CairoMakie
using LinearAlgebra
include("../../src/scripts/MAG_plot.jl")
include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_plot.jl")
using .MAG_plot
using .MAVEN_load
using .MAVEN_plot
const Rm = 3389.5

date = Date(2022, 8, 24)
datestr = Dates.format.(date, "yyyymmdd")
year_num = Dates.format(date, "yyyy")
month_num = Dates.format(date, "mm")

data_dir = "E:\\Data"
output_dir = "E:\\Output"
MAVEN_path = joinpath(data_dir, "MAVEN")
MAG_path = joinpath(MAVEN_path, "MAG", "l2", year_num, month_num)
KP_path = joinpath(MAVEN_path, "KP", "l3", year_num, month_num)
SWIA_path = joinpath(MAVEN_path, "SWIA", "l2", year_num, month_num)
STATIC_path = joinpath(MAVEN_path, "STATIC", "l2", year_num, month_num)
MAVEN_output_path = joinpath(output_dir, "MAVEN")
STATIC_output_path = joinpath(MAVEN_output_path, "STATIC")
KP_output_path = joinpath(MAVEN_output_path, "KP")
MAG_output_path = joinpath(MAVEN_output_path, "MAG")
Overview_path = joinpath(MAG_output_path, "Overview", year_num, month_num)
mkpath(Overview_path)
MAG_file = glob("mvn_mag_l2*$datestr*.sts", MAG_path)
KP_file = glob("mvn_kp*$datestr*.jld2", KP_path)
SWIA_file = glob("mvn_swi_l2_coarsesvy3d*$datestr*.cdf", SWIA_path)
STA_file = glob("mvn_sta_l2*32e4d16a8m*$datestr*.cdf", STATIC_path)
println("MAG file: ", MAG_file)
println("KP file: ", KP_file)
println("SWIA file: ", SWIA_file)
println("STATIC file: ", STA_file)
println("Output path: ", Overview_path)

mag_data = load_mag_l2(MAG_file[1])
mag_data[:position] = mag_data[:position]./Rm
println("Data Loaded.")

shock_time = DateTime(date, Time(7, 27, 20))
wave_time = DateTime(date, Time(7, 30, 18))
shockJUL = datetime2julian.(shock_time)
waveJUL = datetime2julian.(wave_time)
ind_shock = findmin(abs.(mag_data[:epoch].-shock_time))[2]
shock_pos = mag_data[:position][ind_shock, :]

function bowshock(xshock)
    xF = 0.55
    ϵ  = 1.05
    L  = 2.10
    temp = (ϵ^2 - 1.0) * (xshock - xF)^2 -
           2ϵ * L * (xshock - xF) +
           L^2
    return temp >= 0 ? sqrt(temp) : Inf64
end

using Optim

function nearest_bowshock_3D(x0, y0, z0;
                             xmin = -5.0,
                             xmax = 5.0)
    # MAVEN 到 x 轴的距离
    r0 = sqrt(y0^2 + z0^2)
    # distance^2
    function dist2(x)
        rs = bowshock(x)
        if !isfinite(rs)
            return Inf
        end
        return (x - x0)^2 + (rs - r0)^2
    end
    result = optimize(dist2, xmin, xmax)

    xs = Optim.minimizer(result)
    rs = bowshock(xs)

    # 最近点的 y,z
    if r0 > 0
        ys = rs * y0 / r0
        zs = rs * z0 / r0
    else
        # MAVEN 正好在 x 轴上
        ys = 0.0
        zs = 0.0
    end
    distance = sqrt(Optim.minimum(result))
    return (
        xshock = xs,
        yshock = ys,
        zshock = zs,
        rshock = rs,
        distance = distance
    )
end

function bowshock_normal_3D(x, y, z)
    xF = 0.55
    ϵ  = 1.05
    L  = 2.10
    dx = x - xF
    # dT/dx
    dTdx = 2 * (ϵ^2 - 1.0) * dx - 2ϵ * L
    # gradient of F
    nx = -dTdx
    ny = 2y
    nz = 2z
    # normalize
    norm_n = sqrt(nx^2 + ny^2 + nz^2)
    return [
        nx / norm_n,
        ny / norm_n,
        nz / norm_n
    ]
end

shock_p = nearest_bowshock_3D(shock_pos[1], shock_pos[2], shock_pos[3])
shock_n_mol = bowshock_normal_3D(shock_p.xshock, shock_p.yshock, shock_p.zshock)
println("model results:=====================")
# println("Shock position: ", shock_p.xshock, ", ", shock_p.yshock, ", ", shock_p.zshock)
println("Shock normal: ", shock_n_mol)

shock_time_dat = shock_time - Dates.Minute(1) .+ Dates.Minute(1) .* range(0, 2)
ind_shock_dat = findall(minimum(shock_time_dat).<= mag_data[:epoch] .<= maximum(shock_time_dat))

magShock = mag_data[:B][ind_shock_dat, :]
bShock = Array(magShock[:, 1:3])
bmShock = MAG_plot.MVA(bShock)
shock_n_dat = bmShock.vectors[:, 1]
println("Data for shock intersection:=====================")
println("Shock normal: ", shock_n_dat)

dot(shock_n_mol, shock_n_dat) / (norm(shock_n_mol) * norm(shock_n_dat))