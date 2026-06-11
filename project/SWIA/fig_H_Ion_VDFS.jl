include("../code/head.jl")
CairoMakie.activate!()

hour_range = [11, 11]
min_range = [31, 34]

time_stemp = [
    DateTime(yyyy, mm, dd, 11, 32, 34 - 2),
    # DateTime(yyyy, mm, dd, 11, 32, 38 - 2),
    DateTime(yyyy, mm, dd, 11, 32, 42 - 2),
    # DateTime(yyyy, mm, dd, 11, 32, 46 - 2),
    DateTime(yyyy, mm, dd, 11, 32, 50 - 2),
    # DateTime(yyyy, mm, dd, 11, 32, 54 - 2),
    DateTime(yyyy, mm, dd, 11, 32, 58 - 2),
    # DateTime(yyyy, mm, dd, 11, 33, 02 - 2),
];

sta_data = data_dict["STATIC_d1"];
swi_data = data_dict["SWIA_coarse_svy_3d"];
time_b = data_dict["MAG_ss1s_l3"][:epoch]
B_ss = data_dict["MAG_ss1s_l3"][:B]
position_mso = data_dict["MAG_ss1s_l3"][:position]

# VDFS
fig_VDF = Figure(; size=(9*inch, 9*inch))

# Label(fig_VDF[0, :], "$spice_name", valign=:top, halign=:center, tellwidth=false)
VDFS_data = Dict()
df_unit = ["(cm",superscript("-2")," s",superscript("-1")," sr",superscript("-1"),")"]
colorbar_label = Dict(
    "O2+" =>  rich(
        O2_label," ",
        PSD_label),
    "O+" =>   rich(
        O_label," ",
        PSD_label),
    "H+" =>   rich(
        H_label," ",
        PSD_label),
        )
spice = Dict(
    :mass_range => [0, 2],
    :m_int => 1,
)
ion_vel = 160.0
spice_name = "H+"
colorrange = (1e-13, 1e-10)
# ------------------STA part-------------------------
V_30ev = MAVEN_STATIC.ion_energy2v(30, spice[:m_int]) / 1000
V_494ev = MAVEN_STATIC.ion_energy2v(494, spice[:m_int]) / 1000
for II in eachindex(time_stemp)
    local plot_range = (-300, 300)
    ticks = [-200, 0, 200]
    println("Doing $(time_stemp[II]) $spice_name")
    local time_ind = find_time(data_dict["STATIC_d1"][:epoch],time_stemp[II])
    local sta_time = data_dict["STATIC_d1"][:epoch][time_ind]
    local sta_slip_data = MAVEN_STATIC.static_slip(sta_data, time_ind)
    local sta_slip_data = MAVEN_STATIC.STA_count2df(sta_slip_data)
    # sta_slip_data = MAVEN_STATIC.static_rotation(sta_slip_data; frame="MSO")

    local rotation_Q = QuaternionF64(sta_slip_data[:quat_mso]...);

    local time_b_ind = find_time(sta_time, time_b)
    local vsc = (position_mso[time_b_ind+1, :] .- position_mso[time_b_ind-1, :]) ./ (datetime2unix(time_b[time_b_ind+1]) - datetime2unix(time_b[time_b_ind-1]))
    vsc = MAVEN_STATIC.rotate_vector_with_quat_reverse(vsc,rotation_Q)

    local b0 = B_ss[time_b_ind, :]
    b0 = MAVEN_STATIC.rotate_vector_with_quat_reverse(b0,rotation_Q)

    local vel, _, den = MAVEN_STATIC.sta_v_4d(sta_slip_data; energy_range=[30, 1e8], spice...)
    vel = vel .+ vsc

    local V_data = MAVEN_STATIC.static_slip_2_V(sta_slip_data; spice...)
    v = V_data[:v]
    v0 = sqrt.(sum(v[:, :, :] .^ 2, dims=3))[:, :, 1]
    V_data[:v0] = v0
    V_data[:x0] = MAVEN_STATIC.rotate_vector_with_quat_reverse(sta_slip_data[:pos_sc_mso],rotation_Q)

    x0 = V_data[:x0]

    v0 = reshape(V_data[:v0], 32 * 64)
    v = reshape(V_data[:v], 32 * 64, 3)
    dF = reshape(V_data[:df], 32 * 64) .* 1e-3 # IS化 # unit sec^3/(cm^-3 km^-3) to sec^3 m^-6
    bins_sc = reshape(sta_slip_data[:bins_sc],64,1) .* ones(64,32)
    mask_bins_sc = reshape(bins_sc, 32 * 64)

    keys_range = [
            :energy,:denergy,:theta,:dtheta,:phi,:dphi,
        ]
    local range_data =Dict(
        :shape => size(V_data[:df]),
        :mode => "sta",
        :m_int => spice[:m_int],
    )

    title = Dates.format(sta_time, "HH:MM:SS")
    # ax1 = Axis(fig_VDF[1+JJ2, II], aspect=2, limits=((plot_range[1], plot_range[2]), (0, plot_range[2])), xlabel=X_label, ylabel=Y_label, yticks=ticks, xticks=ticks, alignmode=Inside())
    # hidexdecorations!(ax1, grid=false, ticks=false, label=false)
    ax2 = Axis(fig_VDF[1, II], title=title, aspect=1, yticks=ticks, xticks=ticks, alignmode=Inside())
    hidexdecorations!(ax2, grid=false, ticks=false)
    ax3 = Axis(fig_VDF[2, II], aspect=1, yticks=ticks, xticks=ticks, alignmode=Inside())
    linkxaxes!(ax2, ax3)
    # hidethetadecorations!(ax1);hiderdecorations!(ax1)

    # ax2, rot2 = MAVEN_plot.STA_2d_slip(ax2, V_data; frame="xy", colorrange=colorrange, vbluk=vel, plot_range=plot_range, vsc=vsc, return_rot_matrix=true, angle_range=[-30, 30], energy_range=[0, 1e4])
    # ax3, rot3 = MAVEN_plot.STA_2d_slip(ax3, V_data; frame="xz", colorrange=colorrange, vbluk=vel, plot_range=plot_range, vsc=vsc, return_rot_matrix=true, angle_range=[-30, 30], xlabel=X_vdf_label, energy_range=[0, 1e4])
    ax2, rot2, ind_mask = MAVEN_plot.VDF_2d_slip(
        ax2,v,dF; 
        normal_vectors=MAVEN_STATIC.rotate_vector_with_quat_reverse.([[1,0,0],[0,1,0]],rotation_Q), 
        colorrange=colorrange, vbluk=vel, plot_range=plot_range, vsc=vsc,magf=b0,return_rot_matrix=true, angle_range=[-30, 30])
    for i in keys_range
        local val = reshape(sta_slip_data[i][1,:,:],32*64)
        range_data[i] = val[ind_mask.*mask_bins_sc .==1]
    end
    ax2 = MAVEN_plot.VDF_2d_mask(ax2,range_data;rot = rot2, vsc=vsc)
    ax3, rot3, ind_mask = MAVEN_plot.VDF_2d_slip(
        ax3,v,dF; 
        normal_vectors=MAVEN_STATIC.rotate_vector_with_quat_reverse.([[1,0,0],[0,0,1]],rotation_Q), 
        colorrange=colorrange, vbluk=vel, plot_range=plot_range, vsc=vsc,magf=b0, return_rot_matrix=true, angle_range=[-30, 30],xlabel=X_vdf_label)
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
        lines!(ax, [0, xx[1]], [0, xx[2]], color=:red, linewidth=3, linestyle=:dash) # 火星方向
        vel_local = rott * vel
        # scatter!(ax, vel_local[1], vel_local[2], color=:white, markersize=15,marker='X') # 速度矢量
        # bb = rott * b0 * 100
        # lines!(ax, [0, bb[1]], [0, bb[2]], color=:white, linewidth=3, linestyle=:dash)
        lines!(ax, Circle(Point2f(0, 0), V_30ev), color=:white) # 30eV的能量圈
        lines!(ax, Circle(Point2f(0, 0), V_494ev), color=:white) 
        arrows!(ax, [Point2f(0, 0)],[Vec2f(vel_local[1],vel_local[2])]; color=:blue,arrowsize=15,arrowtail=15,linewidth=3,align = :endpoint) # 速度矢量arrowtail

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
    # hidexdecorations!(ax3, grid=false, ticks=false)
end
Colorbar(fig_VDF[1:2, 5], limits=colorrange, label=colorbar_label[spice_name], colormap=:viridis, scale=log10)

# ------------------SWI part-------------------------
swi_data = data_dict["SWIA_coarse_svy_3d"]
swi_data = MAVEN_SWIA.eflux2df!(swi_data)
for II in eachindex(time_stemp)
    plot_range = (-800, 800)
    ticks = [-600, 0, 600]
    println("Doing $(time_stemp[II]) $spice_name")
    "SWIA_coarse_svy_3d","SWIA_quat"
    time_ind = find_time(swi_data[:epoch],time_stemp[II])
    swi_time = swi_data[:epoch][time_ind]

    local rotation_Q = QuaternionF64(data_dict["SWIA_quat"][:quat][find_time(data_dict["SWIA_quat"][:epoch],swi_time)]...)

    time_b_ind = find_time(swi_time, time_b)
    local vsc = (position_mso[time_b_ind+1, :] .- position_mso[time_b_ind-1, :]) ./ (datetime2unix(time_b[time_b_ind+1]) - datetime2unix(time_b[time_b_ind-1]))
    vsc = MAVEN_STATIC.rotate_vector_with_quat_reverse(vsc,rotation_Q)
    local p0 = position_mso[time_b_ind, :]
    p0 = MAVEN_STATIC.rotate_vector_with_quat_reverse(p0,rotation_Q)

    b0 = B_ss[time_b_ind, :]
    b0 = MAVEN_STATIC.rotate_vector_with_quat_reverse(b0,rotation_Q)

    vel, _, den = MAVEN_SWIA.v_3d(swi_data,time_ind;energy_range=[494,1e8])
    vel =vel[1,:] .+ vsc

    local energy = reshape(swi_data[:energy],1,1,48)
    local theta  = reshape(swi_data[:theta][time_ind,:,:],1,4,48)
    local phi    = reshape(swi_data[:phi],16,1,1)
    v0 = MAVEN_SWIA.ion_energy2v.(energy,1)./1000

    vv = vec(MAVEN_SWIA.sphere2xyz_for_SWIA.(v0,theta,phi))
    # vv = MAVEN_SWIA.rotate_vector_with_quat.(vv,quat)

    v = hcat([x[1] for x in vv],[x[2] for x in vv],[x[3] for x in vv])

    dF = vec(swi_data[:df][time_ind,:,:,:]) .* 1e-3 # IS化

    keys_range = [
            :energy,:denergy,:theta,:dtheta,:phi,:dphi,
    ]
    shape_matrix = ones(16,4,48)
    local range_data =Dict(
        :shape => ones(size(swi_data[:df][time_ind,:,:,:])),
        :mode => "swi",
        :m_int => spice[:m_int],
        :energy => energy.*shape_matrix,
        :denergy => reshape(swi_data[:denergy],1,1,48).*shape_matrix,
        :theta => theta.*shape_matrix,
        :dtheta => reshape(swi_data[:dtheta][time_ind,:,:],1,4,48).*shape_matrix,
        :phi => reshape(phi,16,1,1).*shape_matrix,
        :dphi => swi_data[:dphi].*shape_matrix,
    )
    for i in keys_range
        range_data[i] = reshape(range_data[i],16*4*48)
    end

    title = Dates.format(swi_time, "HH:MM:SS")
    ax2 = Axis(fig_VDF[1+2, II], title=title, aspect=1, yticks=ticks, xticks=ticks, alignmode=Inside())
    hidexdecorations!(ax2, grid=false, ticks=false)
    ax3 = Axis(fig_VDF[2+2, II], aspect=1, yticks=ticks, xticks=ticks, alignmode=Inside())
    linkxaxes!(ax2, ax3)
    ax2, rot2, ind_mask = MAVEN_plot.VDF_2d_slip(
        ax2,v,dF; 
        normal_vectors=MAVEN_SWIA.rotate_vector_with_quat_reverse.([[1,0,0],[0,1,0]],rotation_Q), 
        colorrange=colorrange, vbluk=vel, plot_range=plot_range, vsc=vsc,magf=b0,return_rot_matrix=true, angle_range=[-30, 30])
    local masked_range_data = Dict(
        :shape => ones(size(swi_data[:df][time_ind,:,:,:])),
        :mode => "swi",
        :m_int => spice[:m_int],
    )
    for i in keys_range
        masked_range_data[i] = range_data[i][ind_mask]
    end
    ax2 = MAVEN_plot.VDF_2d_mask(ax2,masked_range_data;rot = rot2, vsc=vsc)

    ax3, rot3, ind_mask = MAVEN_plot.VDF_2d_slip(
        ax3,v,dF; 
        normal_vectors=MAVEN_SWIA.rotate_vector_with_quat_reverse.([[1,0,0],[0,0,1]],rotation_Q), 
        colorrange=colorrange, vbluk=vel, plot_range=plot_range, vsc=vsc,magf=b0, return_rot_matrix=true, angle_range=[-30, 30],xlabel=X_vdf_label)
    local masked_range_data = Dict(
        :shape => ones(size(swi_data[:df][time_ind,:,:,:])),
        :mode => "swi",
        :m_int => spice[:m_int],
    )
    for i in keys_range
        masked_range_data[i] = range_data[i][ind_mask]
    end
    ax3 = MAVEN_plot.VDF_2d_mask(ax3,masked_range_data;rot = rot3, vsc=vsc)

    for i in 1:2
        rott = [rot2, rot3][i]
        ax = [ax2, ax3][i]
        xx = -p0
        xx = rott * xx * 100
        lines!(ax, [0, xx[1]], [0, xx[2]], color=:red, linewidth=3, linestyle=:dash) # 火星方向
        # bb = rott * b0 * 100
        vel_local = rott * vel
        # scatter!(ax, vel_local[1], vel_local[2], color=:white, markersize=15, marker='X') # 速度矢量
        # lines!(ax, [0, bb[1]], [0, bb[2]], color=:white, linewidth=3, linestyle=:dash)
        lines!(ax, Circle(Point2f(0, 0), V_30ev), color=:white) # 30eV的能量圈
        lines!(ax, Circle(Point2f(0, 0), V_494ev), color=:white) 
        arrows!(ax, [Point2f(0, 0)],[Vec2f(vel_local[1],vel_local[2])]; color=:blue,arrowsize=15,arrowtail=15,linewidth=3,align = :endpoint) # 速度矢量arrowtail

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
end

Colorbar(fig_VDF[1+2:2+2, 5], limits=colorrange, label=rich("i",superscript("+")," ",PSD_label), colormap=:viridis, scale=log10)

for np in 1:4*4
    text_color = :white
    text_char = Char('a' + np - 1)
    text_char = "($text_char)"
    JJ = mod(np-1 , 4) + 1
    II = div(np-1 , 4) + 1
    text!(fig_VDF[II, JJ], 0, 1, text=text_char, font=:bold, align=(:left, :top),offset=(4, -2), space=:relative,color = text_color)
end
# save("fig_code/fig_output/Fig_H_VDF.png", fig_VDF)
save("publish_data/trace_change_R3/Fig6.eps", fig_VDF)


# E = data_dict["SWIA_coarse_svy_3d"][:energy]
# dE = data_dict["SWIA_coarse_svy_3d"][:denergy]

# E1 = E 
# E2 = E .+ dE
# E1 = E .- dE./2
# E2 = E .+ dE./2

# yy = zeros(48)#1:48
# fig = Figure()
# ax = Axis(fig[1,1],xscale=log10)
# for i in 1:48
#     poly!(ax, [E1[1,i], E2[1,i], E2[1,i], E1[1,i]], [yy[i]-0.05, yy[i]-0.05, yy[i]+0.05, yy[i]+0.05], color = :gray, alpha=0.5)
# end
# scatter!(ax, E1[1,:],yy, color = :red,marker = '|',markersize=10)
# scatter!(ax, E2[1,:],yy, color = :blue,marker = '|',markersize=10)
# # scatter!(ax, E[1,:],1:48, color = :black)
# display(fig)