using Glob
using CommonDataFormat
using Dates
using Statistics
using CairoMakie
include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_plot.jl")
using .MAVEN_load
using .MAVEN_plot

# ========== 配置：日期范围 ==========
start_date = Date(2022, 8, 18)
end_date   = Date(2022, 8, 28)
dates = start_date:Dates.Day(1):end_date
# ===================================

# ========== 弓激波筛选配置 ==========
sw_bowshock_margin = 0.1   # 卫星需在弓激波外多少 Rm
# ===================================

data_dir = "E:\\Data"
MAVEN_path = joinpath(data_dir, "MAVEN")

# ---------- 弓激波模型 ----------
function bowshock(xshock)
    xF = 0.55
    ϵ  = 1.05
    L  = 2.10
    temp = (ϵ^2-1.0)*(xshock-xF)^2 - 2ϵ*L*(xshock-xF) + L^2
    return temp >= 0 ? sqrt(temp) : Inf64
end

# ---------- 加载数据 ----------
all_epoch = DateTime[]
all_MSO_x = Float32[]
all_MSO_y = Float32[]
all_MSO_z = Float32[]

for date in dates
    datestr = Dates.format(date, "yyyymmdd")
    year_num = Dates.format(date, "yyyy")
    month_num = Dates.format(date, "mm")
    KP_path = joinpath(MAVEN_path, "KP", "l3", year_num, month_num)
    KP_file = glob("mvn_kp*$datestr*.jld2", KP_path)
    if isempty(KP_file)
        @warn "No KP file for $datestr, skipping"
        continue
    end
    println("Loading: ", KP_file[1])
    data = MAVEN_load.load_kp_l3(KP_file[1])
    append!(all_epoch, data[:time])
    append!(all_MSO_x, data[:vars][:var_190])
    append!(all_MSO_y, data[:vars][:var_191])
    append!(all_MSO_z, data[:vars][:var_192])
end

if isempty(all_epoch)
    error("No data loaded!")
end

# ---------- 转换单位 ----------
Rm = 3393.5  # 火星半径 km
mso_x = all_MSO_x ./ Rm
mso_y = all_MSO_y ./ Rm
mso_z = all_MSO_z ./ Rm
mso_r = sqrt.(mso_y.^2 .+ mso_z.^2)  # 距 X 轴距离

# ---------- 弓激波筛选 ----------
rho_bs_raw = bowshock.(mso_x)                     # 弓激波半径（无余量）
is_upstream = isinf.(rho_bs_raw)                   # Inf = 日出方向上游，无条件为太阳风
rho_bs = rho_bs_raw .+ sw_bowshock_margin          # 含余量的弓激波半径

# 非上游区域：X-R、X-Y、X-Z 三个方向都需在弓激波外
mask_r  = mso_r .> rho_bs
mask_y  = abs.(mso_y) .> rho_bs
mask_z  = abs.(mso_z) .> rho_bs
mask_outside = mask_r .& mask_y .& mask_z

sw_mask = is_upstream .| mask_outside
n_sw = count(sw_mask)
n_all = length(sw_mask)
println("Solar wind: $n_sw / $n_all ($(round(100*n_sw/n_all, digits=1))%)")
println("  上游(Inf): $(count(is_upstream)) | 弓激波外: $(count(mask_outside)) | X-R: $(count(mask_r)) | X-Y: $(count(mask_y)) | X-Z: $(count(mask_z))")

# ---------- 弓激波理论曲线（用于绘图参考）----------
x_bs = collect(range(-3.0, 3.0, length=500))
r_bs = bowshock.(x_bs)
# 只画有限值部分（过滤 Inf）
valid_bs = isfinite.(r_bs)
x_bs_plot = x_bs[valid_bs]
r_bs_plot = r_bs[valid_bs]

# ---------- 绘图 ----------
CairoMakie.activate!()
fig = Figure(size = (1500, 500), fontsize = 16)

# 三类数据：磁鞘(灰)、上游无条件太阳风(蓝)、弓激波外太阳风(橙)
not_sw = .!sw_mask
upstream_sw = sw_mask .& is_upstream
outside_sw = sw_mask .& .!is_upstream

# X-R 图
ax1 = Axis(fig[1,1], xlabel = "X (Rm)", ylabel = "R = √(Y²+Z²) (Rm)",
    title = "Solar wind  (X-R)  —  3D outside bow shock")
scatter!(ax1, mso_x[not_sw],      mso_r[not_sw];      color = (:gray, 0.3), markersize=2, label="Magnetosheath")
scatter!(ax1, mso_x[upstream_sw], mso_r[upstream_sw]; color = (:dodgerblue, 0.6), markersize=3, label="Upstream (x > nose)")
scatter!(ax1, mso_x[outside_sw],  mso_r[outside_sw];  color = (:orange, 0.6), markersize=3, label="Outside bow shock")
lines!(ax1, x_bs_plot, r_bs_plot; color = :red, linewidth=2, label="Bow shock")
lines!(ax1, x_bs_plot, r_bs_plot .+ sw_bowshock_margin; color = :red, linewidth=1, linestyle=:dash, label="+ $(sw_bowshock_margin)Rm margin")
axislegend(ax1, position = :lt)

# X-Y 图
ax2 = Axis(fig[1,2], xlabel = "X (Rm)", ylabel = "Y (Rm)",
    title = "Satellite Position (X-Y)")
scatter!(ax2, mso_x[not_sw],      mso_y[not_sw];      color = (:gray, 0.3), markersize=2)
scatter!(ax2, mso_x[upstream_sw], mso_y[upstream_sw]; color = (:dodgerblue, 0.6), markersize=3)
scatter!(ax2, mso_x[outside_sw],  mso_y[outside_sw];  color = (:orange, 0.6), markersize=3)
lines!(ax2, x_bs_plot,  r_bs_plot; color = :red, linewidth=2, label="Bow shock")
lines!(ax2, x_bs_plot, -r_bs_plot; color = :red, linewidth=2)
lines!(ax2, x_bs_plot,  r_bs_plot .+ sw_bowshock_margin; color = :red, linewidth=1, linestyle=:dash)
lines!(ax2, x_bs_plot, -(r_bs_plot .+ sw_bowshock_margin); color = :red, linewidth=1, linestyle=:dash)
axislegend(ax2, position = :lt)

# X-Z 图
ax3 = Axis(fig[1,3], xlabel = "X (Rm)", ylabel = "Z (Rm)",
    title = "Satellite Position (X-Z)")
scatter!(ax3, mso_x[not_sw],      mso_z[not_sw];      color = (:gray, 0.3), markersize=2)
scatter!(ax3, mso_x[upstream_sw], mso_z[upstream_sw]; color = (:dodgerblue, 0.6), markersize=3)
scatter!(ax3, mso_x[outside_sw],  mso_z[outside_sw];  color = (:orange, 0.6), markersize=3)
lines!(ax3, x_bs_plot,  r_bs_plot; color = :red, linewidth=2, label="Bow shock")
lines!(ax3, x_bs_plot, -r_bs_plot; color = :red, linewidth=2)
lines!(ax3, x_bs_plot,  r_bs_plot .+ sw_bowshock_margin; color = :red, linewidth=1, linestyle=:dash)
lines!(ax3, x_bs_plot, -(r_bs_plot .+ sw_bowshock_margin); color = :red, linewidth=1, linestyle=:dash)
axislegend(ax3, position = :lt)

# 所有子图 x 范围一致
x_lo = minimum(mso_x)
x_hi = max(3.0, maximum(mso_x))
xlims!(ax1, x_lo, x_hi)
xlims!(ax2, x_lo, x_hi)
xlims!(ax3, x_lo, x_hi)

output_path = joinpath("E:\\Output", "MAVEN", "KP", "Overview",
    Dates.format(start_date, "yyyy"), Dates.format(start_date, "mm"))
mkpath(output_path)
name = "Position_check_" * Dates.format(start_date, "yyyymmdd") * 
       "_to_" * Dates.format(end_date, "yyyymmdd") * ".pdf"
save(joinpath(output_path, name), fig, pt_per_unit = 1)
println("Saved: ", joinpath(output_path, name))
