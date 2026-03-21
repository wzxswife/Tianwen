"""
画图示例：弓激波法向与太阳风磁场夹角分析
"""

using Dates
using DataFrames
using Statistics
using LinearAlgebra

include("bowshock_angle.jl")
include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/TW_load.jl")
using .TW_load
using .MAVEN_load
using .BowshockAngle

# 配置参数
const Data_PATH = joinpath(@__DIR__, "..", "..", "Data", "MAVEN")
const DATE_STR = "20220824"
const THRESHOLD = 0.5  # 弓激波附近距离阈值 (Rm)

# 常量定义
const Rm = 3390.0  # 火星半径 (km)
# 弓激波模型参数 
const BS_XF = 0.55
const BS_EPSILON = 1.05
const BS_L = 2.10

# 加载数据
filename = "mvn_mag_l2_2022236ss_20220824_v01_r01.sts"
filepath = joinpath(Data_PATH, filename)
mag2c32hz = load_mag_l2(filepath)

shock_time = DateTime("2022-08-24T07:27:20")
timeShock = shock_time - Dates.Minute(1) .+ Dates.Minute(1) .* range(0, 2)
timeSolar = shock_time .- Dates.Minute(1) .* range(0,5)
indShock = findall((minimum(timeShock) .<= mag2c32hz[:epoch] .<= maximum(timeShock))
    .& .!isnan.(mag2c32hz[:position][:, 1])
    .& .!isnan.(mag2c32hz[:position][:, 2])
    .& .!isnan.(mag2c32hz[:position][:, 3]))
indSolar = findall((minimum(timeSolar) .<= mag2c32hz[:epoch] .<= maximum(timeSolar))
    .& .!isnan.(mag2c32hz[:position][:, 1])
    .& .!isnan.(mag2c32hz[:position][:, 2])
    .& .!isnan.(mag2c32hz[:position][:, 3]))

function MVA(magbc)
    nb = length(magbc[:, 1])
    bm = mean(magbc, dims=1)
    muv = magbc' * magbc ./ nb - bm' * bm
    return eigen(muv)
end

function angle_between(v1, v2)
    cosθ = dot(v1, v2) / (norm(v1) * norm(v2))
    return acosd(clamp(cosθ, -1.0, 1.0))  # 限制范围防止数值误差
end

magShock = mag2c32hz[:B][indShock, :]
magSolar = mag2c32hz[:B][indSolar, :]
bShock = Array(magShock[:, 1:3])
bSolar = Array(magSolar[:, 1:3])
bmShock = MVA(bShock)
BmeanSolar = dropdims(mean(bSolar, dims=1), dims=1)
kShock = bmShock.vectors[:, 1]
ShockAngle = angle_between(BmeanSolar, kShock)
println("激波法向: n = ($(kShock[1]), $(kShock[2]), $(kShock[3]))")
println("弓激波法向量与太阳风磁场的夹角: $(ShockAngle)°")

_, time_ind = findmin(d -> abs(d - shock_time), mag2c32hz[:epoch])
x0, y0, z0 = mag2c32hz[:position][time_ind, :] ./ Rm
result = BowshockAngle.nearest_point_on_bowshock(x0, y0, z0)
nShock = bowshock_normal_analytic(result.x, result.y, result.z)
println("该点法向量: n = ($(nShock[1]), $(nShock[2]), $(nShock[3]))")
angleBowShock = angle_between(BmeanSolar, nShock)
println("法向量与太阳风磁场的夹角: $(angleBowShock)°")

"bow-shock model"
function bowshock(xshock)
    xF = 0.55 # Rm
    ϵ = 1.05
    L = 2.10 # rm
    rSD = 1.58
    temp = (ϵ^2-1.0)*(xshock-xF)^2-2ϵ*L*(xshock-xF)+L^2
    if temp>=0 
        return sqrt(temp)
    else
        return Inf64
    end
end
"magnetopause model"
function magnetopause(xmp)
    rSD = 1.33
	xF = 0.86
	ϵ = 0.92
	L = 0.90
    temp = (ϵ^2-1.0)*(xmp-xF)^2-2ϵ*L*(xmp-xF)+L^2
    if temp>=0 
        return sqrt(temp)
    else
        return Inf64
    end
end

# 计算弓激波和磁层边界位置
np = 5001
xshock = range(-10, 3, length=np)
ryzshock = bowshock.(xshock)
ind1 = findall(isfinite, ryzshock)
xshock = [xshock[ind1]; reverse(xshock[ind1])]
ryzshock = [ryzshock[ind1]; -reverse(ryzshock[ind1])]
xmp = range(-10, 3, length=np)
ryzmp = magnetopause.(xmp)
ind1 = findall(isfinite, ryzmp)
xmp = [xmp[ind1]; reverse(xmp[ind1])]
ryzmp = [ryzmp[ind1]; -reverse(ryzmp[ind1])]

OUTPUT_PATH = @__DIR__
#overview
using CairoMakie, GeometryBasics
CairoMakie.activate!()
println("Drawing...")
fig = Figure(resolution=(800, 1200))
title = "Bowshock Analysis - MAVEN 2022-08-24"

ax1 = Axis(fig[1,1], aspect=DataAspect(), xlabel=L"$x$ ($R_{\mathrm{M}}$)", 
    ylabel=L"$\sqrt{y^2+z^2}$ ($R_{\mathrm{M}}$)", xticks=range(-10, 10), yticks=range(-10, 10))
xlims!(ax1, -5, 5) 
ylims!(ax1, 0, 5) 
a1 = arc!(ax1, Point2f(0), 1, 0, π; color=:black, linewidth=2)
phi = range(0, 0.5π; length=180)
x = [cos.(phi); 0.0]
y = [sin.(phi); 0.0]
pn = Polygon(Point2f[(xi, yi) for (xi, yi) in zip(x, y)])
poly!(ax1, pn)
l1 = lines!(ax1, xshock, ryzshock; linewidth=2, overdraw = true)
l2 = lines!(ax1, xmp, ryzmp; linewidth=2, overdraw = true)
l3 = lines!(ax1, mag2c32hz[:position][:, 1]./ Rm, 
    sqrt.(mag2c32hz[:position][:, 2].^2+mag2c32hz[:position][:, 3].^2)./ Rm; 
    color=:black, overdraw = true)
scatter!(ax1, mag2c32hz[:position][time_ind, 1]./ Rm, 
    sqrt.(mag2c32hz[:position][time_ind, 2].^2+mag2c32hz[:position][time_ind, 3].^2)./ Rm, 
    markersize=20)
# 太阳风箭头、激波法向箭头
arrow_length = 1.5
B_unit = BmeanSolar ./ norm(BmeanSolar)
end_x = mag2c32hz[:position][time_ind, 1]./ Rm
end_r = sqrt.(mag2c32hz[:position][time_ind, 2].^2+mag2c32hz[:position][time_ind, 3].^2)./ Rm
start_r = end_r - sqrt.(B_unit[2].^2+B_unit[3].^2) * arrow_length
start_x = end_x + B_unit[1] * arrow_length 
dbx = -B_unit[1] * arrow_length 
dbr = sqrt.(B_unit[2].^2+B_unit[3].^2) * arrow_length
dkx = -nShock[1] * arrow_length
dkr = sqrt.(nShock[2].^2+nShock[3].^2) * arrow_length
arrows2d!(ax1, [start_x], [start_r], 
    [dbx],
    [dbr];  
    color=:purple, shaftwidth=3, tipwidth=8, tiplength=10
)
arrows2d!(ax1, [end_x], [end_r], 
    [dkx],
    [dkr];  
    color=:red, shaftwidth=3, tipwidth=8, tiplength=10
)
text!(ax1, start_x-0.2*dbr, start_r-0.2*dbx, text=L"sw", 
    fontsize=18, color=:purple, align=(:center, :center))

ax2 = Axis(fig[2,1], aspect=DataAspect(), 
    xlabel=L"$x$ ($R_{\mathrm{M}}$)", ylabel=L"$y$ ($R_{\mathrm{M}}$)", 
    xticks=range(-10, 10), yticks=range(-10, 10))
xlims!(ax2, -6, 6) 
ylims!(ax2, -3, 3) 
a2 = arc!(ax2, Point2f(0), 1, 0, 2π; color=:black, linewidth=2)
phi = range(-0.5π, 0.5π; length=180)
x = [cos.(phi); 0.0]
y = [sin.(phi); 0.0]
pn = Polygon(Point2f[(xi, yi) for (xi, yi) in zip(x, y)])
poly!(ax2, pn)
l1 = lines!(ax2, xshock, ryzshock; linewidth=2, overdraw = true)
l2 = lines!(ax2, xmp, ryzmp; linewidth=2, overdraw = true)
l3 = lines!(ax2, mag2c32hz[:position][:, 1]./ Rm, 
    mag2c32hz[:position][:, 2]./ Rm; 
    color=:black, overdraw = true)
scatter!(ax2, mag2c32hz[:position][time_ind, 1]./ Rm, 
    mag2c32hz[:position][time_ind, 2]./ Rm, 
    markersize=20)
# 太阳风箭头
end_x = mag2c32hz[:position][time_ind, 1]./ Rm
end_y = mag2c32hz[:position][time_ind, 2]./ Rm
start_x = end_x + B_unit[1] * arrow_length 
start_y = end_y + B_unit[2] * arrow_length
dby = -B_unit[2] * arrow_length
dky = -nShock[2] * arrow_length
arrows2d!(ax2, 
    [start_x], [start_y],                           
    [dbx],  
    [dby];  
    color=:purple, shaftwidth=3, tipwidth=8, tiplength=10
)
arrows2d!(ax2, 
    [end_x], [end_y],                           
    [dkx],  
    [dky];  
    color=:red, shaftwidth=3, tipwidth=8, tiplength=10
)

ax3 = Axis(fig[3,1], aspect=DataAspect(), 
    xlabel=L"$x$ ($R_{\mathrm{M}}$)", ylabel=L"$z$ ($R_{\mathrm{M}}$)", 
    xticks=range(-10, 10, 11), yticks=range(-10, 10, 11))
xlims!(ax3, -6, 6) 
ylims!(ax3, -3, 3) 
a3 = arc!(ax3, Point2f(0), 1, 0, 2π; color=:black, linewidth=2)
phi = range(-0.5π, 0.5π; length=180)
x = [cos.(phi); 0.0]
y = [sin.(phi); 0.0]
pn = Polygon(Point2f[(xi, yi) for (xi, yi) in zip(x, y)])
poly!(ax3, pn)
l1 = lines!(ax3, xshock, ryzshock; linewidth=2, overdraw = true)
l2 = lines!(ax3, xmp, ryzmp; linewidth=2, overdraw = true)
l3 = lines!(ax3, mag2c32hz[:position][:, 1]./ Rm, mag2c32hz[:position][:, 3]./ Rm, 
    color=:black, overdraw = true)
scatter!(ax3, mag2c32hz[:position][time_ind, 1]./ Rm, mag2c32hz[:position][time_ind, 3]./ Rm, 
    markersize=20)
# 太阳风箭头
end_x = mag2c32hz[:position][time_ind, 1]./ Rm
end_z = mag2c32hz[:position][time_ind, 3]./ Rm
start_x = end_x + B_unit[1] * arrow_length 
start_z = end_z + B_unit[3] * arrow_length 
dbz = -B_unit[3] * arrow_length
dkz = -nShock[3] * arrow_length
arrows2d!(ax3, 
    [start_x], [start_z],                           
    [dbx],  
    [dbz];
    color=:purple, shaftwidth=3, tipwidth=8, tiplength=10
) 
arrows2d!(ax3, 
    [end_x], [end_z],                           
    [dkx],  
    [dkz];  
    color=:red, shaftwidth=3, tipwidth=8, tiplength=10
)

save(joinpath(OUTPUT_PATH, "bowshock_overview.png"), fig)