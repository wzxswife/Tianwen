using CommonDataFormat
using Dates
using Statistics
using CairoMakie

omni_file = joinpath(pwd(), "Data", "MAVEN", "mvn_swi_l2_coarsearc3d_20210824_v02_r00.cdf")
#omni_file = joinpath(pwd(), "../Data/MAVEN/mvn_swi_l2_finearc3d_20210824_v02_r00.cdf")
data = CDFDataset(omni_file)

target_time_str = "2021-08-24T07:30:00"
start_time_str = "2021-08-24T06:00:00"
end_time_str = "2021-08-24T09:00:00"
# 将目标时间字符串解析为 DateTime 对象
target_dt = DateTime(target_time_str, dateformat"yyyy-mm-ddTHH:MM:SS")
start_dt = DateTime(start_time_str, dateformat"yyyy-mm-ddTHH:MM:SS")
end_dt = DateTime(end_time_str, dateformat"yyyy-mm-ddTHH:MM:SS")

flux_4d = data["diff_en_fluxes"] # 维度: [Phi, Theta, Energy, Time]
energies = data["energy_coarse"][:] # 维度: [Energy]
epochs = data["epoch"][:]           # 维度: [Time] (需确保转换为 DateTime)
times_dates = DateTime.(epochs)  # 将 CDF 时间转换为 DateTime 对象数组
target_val, target_idx = findmin(d -> abs(d - target_dt), times_dates)
start_val, start_idx = findmin(d -> abs(d - start_dt), times_dates)
end_val, end_idx = findmin(d -> abs(d - end_dt), times_dates)

# 2. 降维 (变为 Time x Energy)
# 先替换填充值
# replace!(flux_4d, -1.0e31 => NaN)
# 对角度维度求平均，得到 [Energy, Time]
flux_2d = dropdims(mean(flux_4d, dims=(1,2)), dims=(1,2))
# 转置为 [Time, Energy] 以适配 Makie
flux_final = flux_2d' 

CairoMakie.activate!()
time_range = times_dates[start_idx: end_idx]
time_numeric = [datetime2unix(dt) for dt in time_range]
flux_at_target = flux_final[start_idx: end_idx, :]./1000
fig = Figure()
ax = Axis(fig[1, 1])
hm = heatmap!(ax, time_numeric, energies, flux_at_target)
Colorbar(fig[1, 2], hm; label = "Differential Energy Flux (cm² s sr keV)⁻¹", width = 15)
xticks_pos = range(extrema(time_numeric)..., length=6)
xticks_labels = [Dates.format(Dates.unix2datetime(t), "HH:MM") 
                 for t in xticks_pos]
save("MAVEN_i_flux_heatmap.png", fig)