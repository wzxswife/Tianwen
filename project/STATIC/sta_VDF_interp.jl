using Glob
using CommonDataFormat
using Dates
using Statistics, LinearAlgebra
using Rotations
using CairoMakie
include("../../src/scripts/MAVEN_STATIC.jl")
include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_plot.jl")
include("../../src/scripts/SphericalSplines.jl")
using .MAVEN_STATIC
using .MAVEN_load
using .MAVEN_plot
using .SphericalSplines

date = Date(2022, 8, 24)
datestr = Dates.format.(date, "yyyymmdd")
year_num = Dates.format(date, "yyyy")
month_num = Dates.format(date, "mm")

ts = DateTime.(date + Dates.Hour(7) + Dates.Minute(26))
te = ts + Dates.Minute(1) 
# + Dates.Second(30)

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

mag_data = MAVEN_load.load_mag_l2(MAG_file[1])
sta_data = MAVEN_load.load_cdf(STA_file[1])
println("Data loaded.")

df_unit = ["(cm",superscript("-2")," s",superscript("-1")," sr",superscript("-1"),")"]

spice = Dict(
    :mass_range => [0, 2],
    :m_int => 1,
)
spice_name = "H+"
X_vdf_label = rich("V", subscript("x,MSO"), "(km/s)")
Y_vdf_label = rich("V", subscript("y,MSO"), "(km/s)")
Z_vdf_label = rich("V", subscript("z,MSO"), "(km/s)")
# ------------------STA part-------------------------
V_30ev  = MAVEN_STATIC.ion_energy2v(30,  spice[:m_int]) / 1000
V_494ev = MAVEN_STATIC.ion_energy2v(494, spice[:m_int]) / 1000

time_b  = mag_data[:epoch]
position_mso = mag_data[:position]

# ═══════════════════════════════════════════════════════
#  找到 [ts, te] 内所有 STA 时间索引
# ═══════════════════════════════════════════════════════
sta_epoch = sta_data[:epoch]
time_inds = findall(ts .≤ sta_epoch .≤ te)
if isempty(time_inds)
    error("No STA data found in [ts, te]")
end
@info "Found $(length(time_inds)) STA time frames in [$(Dates.format(ts,"HH:MM:SS")), $(Dates.format(te,"HH:MM:SS"))]"

n_frames = length(time_inds)
n_energy = size(sta_data[:energy], 3)
n_bins   = size(sta_data[:energy], 2)

# ═══════════════════════════════════════════════════════════════
#  步骤 1: 定义公共 MSO 方向集（用第一帧的 SC 方向经旋转得到）
# ═══════════════════════════════════════════════════════════════
first_ind = time_inds[1]
slip_first = MAVEN_STATIC.static_slip(sta_data, first_ind)
slip_first = MAVEN_STATIC.STA_count2df(slip_first)
Q_first = QuatRotation(slip_first[:quat_mso])

# 每能道: 公共 MSO 方向 + 累加器
common = Vector{Dict{Symbol,Any}}(undef, n_energy)
for ie in 1:n_energy
    theta_ie = slip_first[:theta][1, :, ie]
    phi_ie   = slip_first[:phi][1, :, ie]

    ct = cosd.(theta_ie)
    st = sind.(theta_ie)
    cp = cosd.(phi_ie)
    sp = sind.(phi_ie)
    dirs_SC = hcat(ct.*cp, ct.*sp, st)

    # 第一帧的 SC 方向 → MSO → 作为"公共 MSO 方向"
    dirs_MSO = similar(dirs_SC)
    for i in 1:size(dirs_SC, 1)
        dirs_MSO[i, :] = Q_first * dirs_SC[i, :]
    end

    common[ie] = Dict(
        :dirs_MSO  => dirs_MSO,   # 公共 MSO 方向 (N×3)
        :n_total   => size(dirs_SC, 1),
        :df_sum    => zeros(Float64, size(dirs_SC, 1)),
        :cnt       => zeros(Int, size(dirs_SC, 1)),
    )
end

# ═══════════════════════════════════════════════════════════════
#  步骤 2: 逐帧构建球面样条 → 插值到公共 MSO 方向 → 累加
# ═══════════════════════════════════════════════════════════════
vsc_sum     = zeros(3)
vel_mso_sum = zeros(3)
n_frames_ok = 0

for time_ind in time_inds
    sta_slip_data = MAVEN_STATIC.static_slip(sta_data, time_ind)
    sta_slip_data = MAVEN_STATIC.STA_count2df(sta_slip_data)
    Q_f = QuatRotation(sta_slip_data[:quat_mso])

    # MAG → vsc
    _, time_b_ind = findmin(d -> abs(d - sta_data[:epoch][time_ind]), time_b)
    vsc = (position_mso[time_b_ind+1, :] .- position_mso[time_b_ind-1, :]) ./
        (datetime2unix(time_b[time_b_ind+1]) .- datetime2unix(time_b[time_b_ind-1]))

    # 体速度 → MSO
    vel, _, _ = MAVEN_STATIC.sta_v_4d(sta_slip_data; energy_range=[30, 1e8], spice...)
    vel_mso = Q_f * (vel .+ vsc)

    # V_data[:df] 已做质量维求和: (nbins, nenergy)
    V_data = MAVEN_STATIC.static_slip_2_V(sta_slip_data; spice...)

    # 逐能道构建球面样条 (输入: SC 方向, 输出: log10(df))
    for ie in 1:n_energy
        theta_ie = sta_slip_data[:theta][1, :, ie]
        phi_ie   = sta_slip_data[:phi][1, :, ie]
        df_ie    = V_data[:df][:, ie] .* 1e-3

        # 筛选有效数据
        valid = .!isnan.(theta_ie) .& .!isnan.(phi_ie) .& .!isnan.(df_ie) .& (df_ie .> 0)
        nv = count(valid)
        nv ≥ 4 || continue

        ct = cosd.(theta_ie[valid])
        st = sind.(theta_ie[valid])
        cp = cosd.(phi_ie[valid])
        sp = sind.(phi_ie[valid])
        dirs_SC_this = hcat(ct.*cp, ct.*sp, st)

        # 球面样条: 在 SC 球面上插值 df
        spline = InterpolationSpline(dirs_SC_this, log10.(df_ie[valid]))
        # 注意: spline 的输入是 SC 方向, 输出是 log10(df)

        # 对每个公共 MSO 方向, 求该帧在该方向的 df
        #   d_mso (公共) → 该帧的 SC 方向: d_sc = Q_f' * d_mso
        #   用 spline 在 d_sc 处插值
        c = common[ie]
        for i in 1:c[:n_total]
            d_mso = c[:dirs_MSO][i, :]
            d_sc  = normalize(Q_f' * d_mso)   # MSO → SC (该帧姿态下)
            log_df_val = spline(Vector(d_sc))
            if !isnan(log_df_val)
                c[:df_sum][i] += 10.0^log_df_val  # 还原线性 df
                c[:cnt][i]    += 1
            end
        end
    end

    vsc_sum     .+= vsc
    vel_mso_sum .+= vel_mso
    global n_frames_ok += 1
end

@info "Accumulated $n_frames_ok frames"

# ═══════════════════════════════════════════════════════════════
#  步骤 3: 平均 → 构建最终球面样条 (MSO 球面)
# ═══════════════════════════════════════════════════════════════
df_avg = zeros(Float64, n_bins, n_energy)
for ie in 1:n_energy
    c = common[ie]
    for i in 1:c[:n_total]
        if c[:cnt][i] > 0
            df_avg[i, ie] = c[:df_sum][i] / c[:cnt][i]
        end
    end
end

vsc_avg     = vsc_sum ./ n_frames_ok
vel_mso_avg = vel_mso_sum ./ n_frames_ok

# 构建最终样条 (与原始代码相似, 但基于正确的平均 df)
splines = Vector{Union{Nothing, SplineSolution}}(undef, n_energy)
for ie in 1:n_energy; splines[ie] = nothing; end

n_spline_ok = 0
for ie in 1:n_energy
    dirs_MSO = common[ie][:dirs_MSO]
    df_ie    = df_avg[:, ie]

    valid = .!isnan.(df_ie) .& (df_ie .> 0)
    nv = count(valid)
    nv ≥ 4 || continue

    # 样条建在 MSO 球面上 (输入 MSO 方向, 输出 log10(df))
    splines[ie] = InterpolationSpline(dirs_MSO[valid, :], log10.(df_ie[valid]))
    global n_spline_ok += 1
end
@info "Built $n_spline_ok/$n_energy MSO splines (corrected avg)"

# 能道速度 (SC 系, 固定不变)
energy_vals = zeros(n_energy)
for ie in 1:n_energy
    energy_vals[ie] = MAVEN_STATIC.ion_energy2v(slip_first[:energy][1, 1, ie], spice[:m_int]) / 1000
end
v_ch = copy(energy_vals)
reverse!(v_ch)
reverse!(splines)

# ═══════════════════════════════════════════
#  步骤 4: MSO 球面样条 → 2D 热图
# ═══════════════════════════════════════════

# 在稠密格点上用 MSO 球面样条插值求 log10(df)
# 改自原 build_vdf_grid: spline 在 MSO 球面上, 直接评估 MSO 方向
# speed 仍需转 SC 来确定能道 (vsc 影响能量)
function build_vdf_grid(Q_ref, vsc_avg, splines, v_ch, n_grid, v_centers; plane="xy")
    grid_log = fill(NaN, n_grid, n_grid)
    n_energy = length(v_ch)
    v_min, v_max = extrema(v_ch)

    for ix in 1:n_grid, iy in 1:n_grid
        if plane == "xy"
            v_mso = [v_centers[ix], v_centers[iy], 0.0]
        else  # "xz"
            v_mso = [v_centers[ix], 0.0, v_centers[iy]]
        end

        # SC 系速度 → 确定能道 (vsc + 姿态投影)
        v_sc = Q_ref' * (v_mso .- vsc_avg)
        speed = norm(v_sc)
        v_min ≤ speed ≤ v_max || continue

        dir_mso = normalize(v_mso)
        any(isnan, dir_mso) && continue

        i_high = searchsortedfirst(v_ch, speed)
        i_low  = i_high - 1
        (i_low ≥ 1 && i_high ≤ n_energy) || continue

        S_lo = splines[i_low]
        S_hi = splines[i_high]
        (S_lo ≡ nothing || S_hi ≡ nothing) && continue

        # MSO 球面样条: 输入 MSO 方向
        log_df_lo = S_lo(Vector(dir_mso))
        log_df_hi = S_hi(Vector(dir_mso))
        (isnan(log_df_lo) || isnan(log_df_hi)) && continue

        # log-log 能道插值
        log_v   = log10(speed)
        log_v_lo = log10(v_ch[i_low])
        log_v_hi = log10(v_ch[i_high])
        t = (log_v - log_v_lo) / (log_v_hi - log_v_lo)
        # Makie: heatmap!(ax, x, y, z) → z[i, j] → (x[i], y[j])
        # ix → Vx (x-axis), iy → Vy (y-axis) → grid_log[ix, iy]
        grid_log[ix, iy] = log_df_lo + t * (log_df_hi - log_df_lo)
    end
    return grid_log
end

# ═══════════════════════════════════════════
#  绘图
# ═══════════════════════════════════════════
fig_name = "MAVEN_STATIC_H+_VDF_interp_avg_" * Dates.format(ts, "yyyymmddTHHMM")*
    "-"*Dates.format(te, "HHMM")*".pdf"
fig_VDF = Figure(size = (1000, 500))
ticks = [-600, 0, 600]
plot_range = (-800, 800)
title_avg = "MAVEN STATIC H⁺ f(v) avg $(datestr) $(Dates.format(ts,"HH:MM"))-$(Dates.format(te,"HH:MM")) MSO"
colorbar_label = rich("STATIC ", spice_name, " VDF (avg) ", df_unit...)

n_grid = 80
v_edges = range(plot_range[1], plot_range[2], length=n_grid + 1)
v_centers = (v_edges[1:end-1] + v_edges[2:end]) / 2

ax2 = Axis(fig_VDF[1, 1], aspect=1, yticks=ticks, xticks=ticks,
    limits=(plot_range, plot_range), alignmode=Inside())
hidexdecorations!(ax2, grid=false, ticks=false)
ax3 = Axis(fig_VDF[1, 2], aspect=1, yticks=ticks, xticks=ticks,
    limits=(plot_range, plot_range), alignmode=Inside())
linkxaxes!(ax2, ax3)

# 先计算两个截面, 再根据数据确定色标范围
grid_xy = build_vdf_grid(Q_first, vsc_avg, splines, v_ch, n_grid, v_centers; plane="xy")
grid_xz = build_vdf_grid(Q_first, vsc_avg, splines, v_ch, n_grid, v_centers; plane="xz")

# 自适应色标: 取两个网格 2%~98% 分位数
all_valid = Float64[v for v in vec(grid_xy) if !isnan(v)]
append!(all_valid, (v for v in vec(grid_xz) if !isnan(v)))
if length(all_valid) > 10
    lo, hi = quantile(all_valid, [0.02, 0.98])
    cr_log = (lo, hi)
else
    cr_log = (-14.0, -9.0)  # fallback 默认 PSD 量级
end

# Vx-Vy 截面 (z=0)
heatmap!(ax2, v_centers, v_centers, grid_xy;
    colormap=:viridis, colorrange=cr_log, nan_color=:transparent)

lines!(ax2, Circle(Point2f(0, 0), V_30ev),  color=:white)
lines!(ax2, Circle(Point2f(0, 0), V_494ev), color=:white)
hlines!(ax2, 0, linestyle=:dash, color=:white)
vlines!(ax2, 0, linestyle=:dash, color=:white)

arrows2d!(ax2, [Point2f(0, 0)], [Vec2f(vel_mso_avg[1], vel_mso_avg[2])];
    color=:cyan, shaftwidth=3, tiplength=15, tipwidth=8)

# Vx-Vz 截面 (y=0)
heatmap!(ax3, v_centers, v_centers, grid_xz;
    colormap=:viridis, colorrange=cr_log, nan_color=:transparent)

lines!(ax3, Circle(Point2f(0, 0), V_30ev),  color=:white)
lines!(ax3, Circle(Point2f(0, 0), V_494ev), color=:white)
hlines!(ax3, 0, linestyle=:dash, color=:white)
vlines!(ax3, 0, linestyle=:dash, color=:white)

arrows2d!(ax3, [Point2f(0, 0)], [Vec2f(vel_mso_avg[1], vel_mso_avg[3])];
    color=:cyan, shaftwidth=3, tiplength=15, tipwidth=8)

ax2.ylabel = Y_vdf_label
ax3.ylabel = Z_vdf_label
hidexdecorations!(ax3, label=false, grid=false, ticks=false, ticklabels=false)

Colorbar(fig_VDF[1, 3], label=colorbar_label,
    limits=cr_log, colormap=:viridis)
Label(fig_VDF[0, :], title_avg, fontsize=20, tellwidth=false)

save(joinpath(VDF_path, fig_name), fig_VDF)
println("Saved figure: ", joinpath(VDF_path, fig_name))
