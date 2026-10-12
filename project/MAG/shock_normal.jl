"""
计算给定时刻MAVEN航天器位置处的弓激波法向量
基于 Edberg et al. (2008) 弓激波模型
"""

using Dates
using LinearAlgebra
using Glob

include("../../src/scripts/MAVEN_load.jl")
using .MAVEN_load

# ============================================================
# 弓激波模型参数 (Edberg et al., 2008)
# ============================================================
const Rm = 3390.0  # 火星半径 (km)
const BS_XF = 0.55       # 焦点x坐标 (Rm)
const BS_EPSILON = 1.05  # 离心率
const BS_L = 2.10        # 半通径 (Rm)

"""
    bowshock_r(x) -> r

弓激波在给定x坐标处的半径 (Rm)。
超出模型范围返回 Inf64。
"""
function bowshock_r(x)
    temp = (BS_EPSILON^2 - 1.0) * (x - BS_XF)^2 -
           2 * BS_EPSILON * BS_L * (x - BS_XF) + BS_L^2
    return temp >= 0 ? sqrt(temp) : Inf64
end

"""
    bowshock_dr_dx(x) -> dr/dx

弓激波半径对x的解析导数。
"""
function bowshock_dr_dx(x)
    r_bs = bowshock_r(x)
    if !isfinite(r_bs) || r_bs < 1e-10
        return Inf64
    end
    numerator = (BS_EPSILON^2 - 1.0) * (x - BS_XF) - BS_EPSILON * BS_L
    return numerator / r_bs
end

"""
    bowshock_normal(x, y, z) -> [nx, ny, nz]

解析计算弓激波表面法向量（单位向量）。

数学原理:
    弓激波隐函数 F(x,y,z) = sqrt(y²+z²) - r(x) = 0
    梯度 ∇F = [-dr/dx, y/r_perp, z/r_perp]
    归一化后即为单位法向量，指向太阳风侧（外侧）。

参数:
    x, y, z: MSO坐标 (Rm)

返回值:
    单位法向量 [nx, ny, nz]
"""
function bowshock_normal(x, y, z)
    r_perp = sqrt(y^2 + z^2)
    r_bs = bowshock_r(x)
    if !isfinite(r_bs)
        return [0.0, 0.0, 0.0]
    end

    drdx = bowshock_dr_dx(x)
    if !isfinite(drdx)
        return [0.0, 0.0, 0.0]
    end

    # 在x轴上 (r_perp ≈ 0)，法向沿x方向
    if r_perp < 1e-10
        nx = drdx > 0 ? 1.0 : -1.0
        return [nx, 0.0, 0.0]
    end

    # 法向量分量（未归一化）
    nx = -drdx
    ny = y / r_perp
    nz = z / r_perp

    # 归一化
    n = [nx, ny, nz] / norm([nx, ny, nz])

    # 确保指向外侧（太阳风侧）
    if n[1] > 0
        n = -n
    end

    return n
end

"""
    nearest_point_on_bowshock(x0, y0, z0) -> (x, y, z, distance)

找到弓激波表面上距离给定点最近的点。
利用旋转对称性简化为一维搜索。
"""
function nearest_point_on_bowshock(x0, y0, z0)
    r0 = sqrt(y0^2 + z0^2)

    if r0 < 1e-10
        best_x = _find_nearest_x_on_axis(x0)
        return (x=best_x, y=0.0, z=0.0, distance=abs(x0 - best_x))
    end

    function f(x)
        r = bowshock_r(x)
        isfinite(r) || return Inf
        return (x - x0) + (r - r0) * bowshock_dr_dx(x)
    end
    function df(x)
        r = bowshock_r(x)
        (!isfinite(r) || r < 1e-10) && return 1.0
        drdx = bowshock_dr_dx(x)
        d2rdx2 = (BS_EPSILON^2 - 1.0 - drdx^2) / r
        return 1.0 + drdx^2 + (r - r0) * d2rdx2
    end

    x = isfinite(bowshock_r(x0)) ? x0 : _find_valid_x(x0)

    for _ in 1:100
        fx = f(x)
        abs(fx) < 1e-10 && break
        dfx = df(x)
        abs(dfx) < 1e-15 && break
        x = clamp(x - fx / dfx, -20.0, 5.0)
    end

    r_bs = bowshock_r(x)
    if !isfinite(r_bs)
        x = _golden_search(x0, r0)
        r_bs = bowshock_r(x)
    end

    y_nearest = r_bs * y0 / r0
    z_nearest = r_bs * z0 / r0
    dist = sqrt((x - x0)^2 + (y_nearest - y0)^2 + (z_nearest - z0)^2)

    return (x=x, y=y_nearest, z=z_nearest, distance=dist)
end

function _find_valid_x(x0)
    for x in range(-15.0, 3.0, length=200)
        isfinite(bowshock_r(x)) && return x
    end
    return 0.0
end

function _find_nearest_x_on_axis(x0)
    phi = (1.0 + sqrt(5.0)) / 2.0
    resphi = 2.0 - phi
    a, b = -15.0, 3.0
    x1, x2 = a + resphi * (b - a), b - resphi * (b - a)
    f1 = isfinite(bowshock_r(x1)) ? abs(x1 - x0) : Inf
    f2 = isfinite(bowshock_r(x2)) ? abs(x2 - x0) : Inf
    for _ in 1:100
        if f1 < f2
            b = x2; x2 = x1; f2 = f1
            x1 = a + resphi * (b - a)
            f1 = isfinite(bowshock_r(x1)) ? abs(x1 - x0) : Inf
        else
            a = x1; x1 = x2; f1 = f2
            x2 = b - resphi * (b - a)
            f2 = isfinite(bowshock_r(x2)) ? abs(x2 - x0) : Inf
        end
        (b - a) < 1e-10 && break
    end
    return (a + b) / 2.0
end

function _golden_search(x0, r0)
    phi = (1.0 + sqrt(5.0)) / 2.0
    resphi = 2.0 - phi
    a, b = -15.0, 3.0
    function obj(x)
        r = bowshock_r(x)
        return isfinite(r) ? (x - x0)^2 + (r - r0)^2 : Inf
    end
    x1, x2 = a + resphi * (b - a), b - resphi * (b - a)
    f1, f2 = obj(x1), obj(x2)
    for _ in 1:100
        if f1 < f2
            b = x2; x2 = x1; f2 = f1
            x1 = a + resphi * (b - a); f1 = obj(x1)
        else
            a = x1; x1 = x2; f1 = f2
            x2 = b - resphi * (b - a); f2 = obj(x2)
        end
        (b - a) < 1e-10 && break
    end
    return (a + b) / 2.0
end

# ============================================================
# 主流程：计算给定时刻的激波法向
# ============================================================

# --- 参数设置 ---
date = Date(2022, 8, 24)
datestr = Dates.format(date, "yyyymmdd")
year_num = Dates.format(date, "yyyy")
month_num = Dates.format(date, "mm")

# 关注时刻（与 sta_VDF.jl 一致）
shock_time = DateTime("2022-08-24T07:26:00")

# 数据路径
data_dir = "E:\\Data"
MAVEN_path = joinpath(data_dir, "MAVEN")
MAG_path = joinpath(MAVEN_path, "MAG", "l2", year_num, month_num)

# 加载磁场数据
MAG_file = glob("mvn_mag_l2*$datestr*.sts", MAG_path)
println("MAG file: ", MAG_file)
mag_data = MAVEN_load.load_mag_l2(MAG_file[1])

# --- 查找最近时刻的位置 ---
_, time_ind = findmin(d -> abs(d - shock_time), mag_data[:epoch])
position_km = mag_data[:position][time_ind, :]  # km, MSO
position_Rm = position_km ./ Rm                  # Rm, MSO

x0, y0, z0 = position_Rm
println("\n========================================")
println("弓激波法向计算 (Edberg et al., 2008 模型)")
println("========================================")
println("时刻: $(mag_data[:epoch][time_ind])")
println("位置 (MSO): x=$(round(x0, digits=3)), y=$(round(y0, digits=3)), z=$(round(z0, digits=3)) Rm")

# --- 到弓激波表面的距离 ---
r_perp = sqrt(y0^2 + z0^2)
r_bs = bowshock_r(x0)
dist = isfinite(r_bs) ? r_perp - r_bs : Inf
println("到弓激波距离: $(round(dist, digits=3)) Rm ($(dist > 0 ? "太阳风侧" : "磁鞘侧"))")

# --- 找弓激波上最近点 ---
nearest = nearest_point_on_bowshock(x0, y0, z0)
println("弓激波上最近点: ($(round(nearest.x, digits=3)), $(round(nearest.y, digits=3)), $(round(nearest.z, digits=3))) Rm")

# --- 计算法向量 ---
normal = bowshock_normal(nearest.x, nearest.y, nearest.z)
println("\n弓激波法向量 (MSO):")
println("  n = ($(round(normal[1], digits=4)), $(round(normal[2], digits=4)), $(round(normal[3], digits=4)))")

# --- 与磁场方向对比 ---
b0 = mag_data[:B][time_ind, :]
println("\n局地磁场 (MSO, nT):")
println("  B = ($(round(b0[1], digits=2)), $(round(b0[2], digits=2)), $(round(b0[3], digits=2)))")
println("  |B| = $(round(norm(b0), digits=2)) nT")

cos_theta = dot(normal, b0) / (norm(b0) + 1e-10)
theta = acosd(clamp(abs(cos_theta), -1.0, 1.0))
println("\n法向量与磁场夹角: $(round(theta, digits=1))°")
