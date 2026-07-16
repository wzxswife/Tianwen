using Glob
using CommonDataFormat
using Dates
using Statistics, LinearAlgebra
using Rotations
using CairoMakie
include("../../src/scripts/MAVEN_STATIC.jl")
include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_plot.jl")
using .MAVEN_STATIC
using .MAVEN_load
using .MAVEN_plot

date = Date(2022, 8, 24)
datestr = Dates.format.(date, "yyyymmdd")
year_num = Dates.format(date, "yyyy")
month_num = Dates.format(date, "mm")

ts = DateTime.(date + Dates.Hour(7) + Dates.Minute(26))
time_stemp = ts .+ Dates.Minute(1) .* range(0, 3)
# time_stemp = ts .+ Dates.Second(30) .* range(0, 3)

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

VDFS_data = Dict()
df_unit = ["(cm",superscript("-2")," s",superscript("-1")," sr",superscript("-1"),")"]

spice = Dict(
    :mass_range => [0, 2],
    :m_int => 1,
)
ion_vel = 160.0
spice_name = "H+"
X_vdf_label = rich("V", subscript("x,MSO"), "(km/s)")
Y_vdf_label = rich("V", subscript("y,MSO"), "(km/s)")
Z_vdf_label = rich("V", subscript("z,MSO"), "(km/s)")
colorrange = (1e-13, 1e-10)
# ------------------STA part-------------------------
V_30ev = MAVEN_STATIC.ion_energy2v(30, spice[:m_int]) / 1000
V_494ev = MAVEN_STATIC.ion_energy2v(494, spice[:m_int]) / 1000

time_b = mag_data[:epoch]
B_ss = mag_data[:B]
position_mso = mag_data[:position]

fig_name = "MAVEN_STATIC_H+_VDF_2D_" * Dates.format(time_stemp[1], "yyyymmddTHHMM")*
    "-"*Dates.format(time_stemp[end], "HHMM")*".pdf"
fig_VDF = Figure(resolution = (1000, 500))
ticks = [-600, 0, 600]
plot_range = (-800, 800)
colorbar_label = rich("STATIC ", spice_name, " VDF ", df_unit...)

for II in eachindex(time_stemp)
    _, time_ind = findmin(d -> abs(d - time_stemp[II]), sta_data[:epoch])
    sta_time = sta_data[:epoch][time_ind]
    sta_slip_data = MAVEN_STATIC.static_slip(sta_data, time_ind)
    sta_slip_data = MAVEN_STATIC.STA_count2df(sta_slip_data)
    rotation_Q = QuatRotation(sta_slip_data[:quat_mso])
    _, time_b_ind = findmin(d -> abs(d - sta_time), time_b)
    vsc = (position_mso[time_b_ind.+1, :] .- position_mso[time_b_ind.-1, :]) ./ 
        (datetime2unix(time_b[time_b_ind.+1]) .- datetime2unix(time_b[time_b_ind.-1]))
    b0 = B_ss[time_b_ind, :]
    # b0 = MAVEN_STATIC.rotate_vector_with_quat_reverse(b0,rotation_Q)
    b0 = rotation_Q' * b0
    vel, _, den = MAVEN_STATIC.sta_v_4d(sta_slip_data; energy_range=[30, 1e8], spice...)
    vel = vel .+ vsc
    V_data = MAVEN_STATIC.static_slip_2_V(sta_slip_data; spice...)
    v = V_data[:v]
    v0 = sqrt.(sum(v[:, :, :] .^ 2, dims=3))[:, :, 1]
    V_data[:v0] = v0
    V_data[:x0] = rotation_Q' * sta_slip_data[:pos_sc_mso]

    x0 = V_data[:x0]
    v0 = reshape(V_data[:v0], 32 * 64)
    v = reshape(V_data[:v], 32 * 64, 3)
    dF = reshape(V_data[:df], 32 * 64) .* 1e-3 # IS化 # unit sec^3/(cm^-3 km^-3) to sec^3 m^-6
    
    # phi_bin_idx = 13
    # energy_range = [1000, 10000] 
    # # ---- phi[13] 以外的 df 置 0 ----
    # nbins = 64
    # nenergy = 32
    # phi_bins = phi_bin_idx:16:nbins      # 所有 phi = phi[13] 的 bin 行索引
    # phi_mask = falses(nbins, nenergy)
    # phi_mask[phi_bins, :] .= true
    # phi_mask_flat = vec(phi_mask)
    # # ---- 能量范围筛选 ----
    # energy_arr = sta_slip_data[:energy][1, :, :]  # (nbins, nenergy)
    # energy_mask = (energy_arr .>= energy_range[1]) .& (energy_arr .<= energy_range[2])
    # energy_mask_flat = vec(energy_mask)
    # # 合并: 只保留 phi[13] 且 在能量范围内的数据
    # keep_mask = phi_mask_flat .& energy_mask_flat
    # # keep_mask = phi_mask_flat
    # dF[.!keep_mask] .= 0.0
    # n_keep = sum(keep_mask)
    # println("$(Dates.format(sta_time, "HH:MM:SS")): $n_keep bins (phi=$(phi_bin_idx), energy $(energy_range[1])-$(energy_range[2]) eV)")

    bins_sc = reshape(sta_slip_data[:bins_sc],64,1) .* ones(64,32)
    mask_bins_sc = reshape(bins_sc, 32 * 64)

    keys_range = [
            :energy,:denergy,:theta,:dtheta,:phi,:dphi,
        ]
    range_data =Dict(
        :shape => size(V_data[:df]),
        :mode => "sta",
        :m_int => spice[:m_int],
    )
    title = Dates.format(sta_time, "HH:MM:SS")
    ax2 = Axis(fig_VDF[1, II], title=title, aspect=1, yticks=ticks, xticks=ticks, alignmode=Inside())
    hidexdecorations!(ax2, grid=false, ticks=false)
    ax3 = Axis(fig_VDF[2, II], aspect=1, yticks=ticks, xticks=ticks, alignmode=Inside())
    linkxaxes!(ax2, ax3)
    ax2, rot2, ind_mask = MAVEN_plot.VDF_2d_slip(
            ax2,v,dF; 
            normal_vectors=[rotation_Q'*v for v in [[1,0,0], [0,1,0]]], 
            colorrange=colorrange, vbluk=vel, plot_range=plot_range, 
            vsc=vsc,magf=b0,return_rot_matrix=true, angle_range=[-30, 30])
    for i in keys_range
        local val = reshape(sta_slip_data[i][1,:,:],32*64)
        range_data[i] = val[ind_mask.*mask_bins_sc .==1]
    end
    ax2 = MAVEN_plot.VDF_2d_mask(ax2,range_data;rot = rot2, vsc=vsc)
    ax3, rot3, ind_mask = MAVEN_plot.VDF_2d_slip(
        ax3,v,dF; 
        normal_vectors=[rotation_Q'*v for v in [[1,0,0], [0,0,1]]], 
        colorrange=colorrange, vbluk=vel, plot_range=plot_range, 
        vsc=vsc,magf=b0, return_rot_matrix=true, angle_range=[-30, 30], 
        xlabel=X_vdf_label)
    for i in keys_range
        local val = reshape(sta_slip_data[i][1,:,:],32*64)
        range_data[i] = val[ind_mask.*mask_bins_sc .==1]
    end
    ax3 = MAVEN_plot.VDF_2d_mask(ax3,range_data;rot = rot3, vsc=vsc)

    for i in 1:2
        rott = [rot2, rot3][i]
        ax = [ax2, ax3][i]
        xx = -x0
        xx = rott * xx * 100
        # lines!(ax, [0, xx[1]], [0, xx[2]], color=:red, linewidth=3, linestyle=:dash) # 火星方向
        vel_local = rott * vel
        # scatter!(ax, vel_local[1], vel_local[2], color=:white, markersize=15,marker='X') # 速度矢量
        # bb = rott * b0 * 100
        # lines!(ax, [0, bb[1]], [0, bb[2]], color=:white, linewidth=3, linestyle=:dash)
        lines!(ax, Circle(Point2f(0, 0), V_30ev), color=:white) # 30eV的能量圈
        lines!(ax, Circle(Point2f(0, 0), V_494ev), color=:white) 
        # arrows2d!(ax, [Point2f(0, 0)],[Vec2f(vel_local[1],vel_local[2])]; 
        #     color=:blue, shaftwidth=3, tiplength = 15, tipwidth  = 8)   #离子漂移速度
        # b_local = rott * b0
        # arrows2d!(ax, [Point2f(0, 0)],[Vec2f(b_local[1], b_local[2]) * 200 / norm(b_local)]; 
        #     color=:cyan, shaftwidth=3, tiplength = 15, tipwidth  = 8)   # 磁场方向
    end
    if II == 1
        # ax1.ylabel = Y_label
        ax2.ylabel = Y_vdf_label
        ax3.ylabel = Z_vdf_label
    else
        # hideydecorations!(ax1, grid=false, ticks=false)
        hideydecorations!(ax2, grid=false, ticks=false)
        hideydecorations!(ax3, grid=false, ticks=false)
    end
    hidexdecorations!(ax3, label=false, grid=false, ticks=false, ticklabels=false)
end
Colorbar(fig_VDF[1:2, 5], label=colorbar_label, 
    limits=colorrange, colormap=:viridis, scale=log10)

save(joinpath(VDF_path, fig_name), fig_VDF)
println("Saved figure: ", joinpath(VDF_path, fig_name))