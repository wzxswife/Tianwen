using CommonDataFormat
using Dates, CSV, DataFrames
using Statistics
using CairoMakie
include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_plot.jl")
include("../../src/scripts/MAVEN_SWIA.jl")
include("../../src/scripts/MAVEN_STATIC.jl")
using .MAVEN_load
using .MAVEN_plot
using .MAVEN_SWIA
using .MAVEN_STATIC

yyyy = 2021
mm = 08
dd = 24
hh = 07
mi = 00
ss = 30

time_stemp = [
    DateTime(yyyy, mm, dd, hh, mi, ss),
    DateTime(yyyy, mm, dd, hh, mi+1, ss),
    DateTime(yyyy, mm, dd, hh, mi+2, ss),  
    DateTime(yyyy, mm, dd, hh, mi+3, ss),
];

date_str = Dates.format.(time_stemp[1], "yyyymmdd")
date_str_dye = "$(year(time_stemp[1]))$(lpad(dayofyear(time_stemp[1]), 3, '0'))"

save_name = "MAVEN_SWIA_Ion_VDF_2D_$(replace(string(time_stemp[1]), ":" => "-")).png"
save_path = joinpath(pwd(), "results", "IonBeam", save_name)
data_path = joinpath(pwd(), "Data", "MAVEN", "mvn_swi_l2_coarsearc3d_$(date_str)_v02_r00.cdf")
quat_path = joinpath(pwd(), "Data", "MAVEN", "mvn_spice_swia_qu_$(date_str).csv")
mag_path = joinpath(pwd(), "Data", "MAVEN", "mvn_mag_l3_$(date_str_dye)ss1s_$(date_str)_v01_r01.f77_unformatted")

data = load_cdf(data_path)
quat_data = load_quat(quat_path)
quat_data[:data_load_flag] = true
swi_data = get_3dc!(data, quat_data = quat_data)
mag_data = load_mag_l3(mag_path)
time_b = mag_data[:epoch]
B_ss = mag_data[:B]
position_mso = mag_data[:position]

fig_VDF = Figure(resolution = (1000, 500))

df_unit = ["(cm",superscript("-2")," s",superscript("-1")," sr",superscript("-1"),")"]

spice = Dict(
    :mass_range => [0, 2],
    :m_int => 1,
)

colorrange = (1e0, 1e8)
spice_name = "H+"
X_vdf_label = "V_x (km/s)"
Y_vdf_label = "V_y (km/s)"
Z_vdf_label = "V_z (km/s)"
V_30ev = MAVEN_STATIC.ion_energy2v(30, spice[:m_int]) / 1000
V_494ev = MAVEN_STATIC.ion_energy2v(494, spice[:m_int]) / 1000

for II in eachindex(time_stemp)
    plot_range = (-800, 800)
    ticks = [-600, 0, 600]
    println("Doing $(time_stemp[II]) $spice_name")
    _, time_ind = findmin(d -> abs(d - time_stemp[II]), swi_data[:epoch])
    swi_time = swi_data[:epoch][time_ind]

    _, rotation_ind = findmin(d -> abs(d - swi_time), quat_data[:epoch])
    local rotation_Q = quat_data[:quat][rotation_ind]

    _, time_b_ind = findmin(d -> abs(d - swi_time), time_b)
    local vsc = (position_mso[time_b_ind+1, :] .- position_mso[time_b_ind-1, :]) ./ (datetime2unix(time_b[time_b_ind+1]) - datetime2unix(time_b[time_b_ind-1]))
    #vsc = MAVEN_STATIC.rotate_vector_with_Matrix(vsc,rotation_Q)
    local p0 = position_mso[time_b_ind, :]
    #p0 = MAVEN_STATIC.rotate_vector_with_Matrix(p0,rotation_Q)
    
    b0 = B_ss[time_b_ind, :]
    #b0 = MAVEN_STATIC.rotate_vector_with_Matrix(b0,rotation_Q)

    vel, _, den = MAVEN_SWIA.v_3d(swi_data,time_ind;energy_range=[494,1e8])
    vel =vel[1,:] .+ vsc

    local energy = reshape(swi_data[:energy],1,1,48)
    local theta  = reshape(swi_data[:theta][time_ind,:,:],1,4,48)
    local phi    = reshape(swi_data[:phi],16,1,1)
    v0 = MAVEN_SWIA.ion_energy2v.(energy,1)./1000
    vv = vec(MAVEN_SWIA.sphere2xyz_for_SWIA.(v0,theta,phi))
    v = hcat([x[1] for x in vv],[x[2] for x in vv],[x[3] for x in vv])

    dF = vec(swi_data[:diff_en_fluxes][time_ind,:,:,:]) .* 1e-3 # IS化

    keys_range = [
        :energy,:denergy,:theta,:dtheta,:phi,:dphi,
    ]
    shape_matrix = ones(16,4,48)
    local range_data =Dict(
        :shape => ones(size(swi_data[:diff_en_fluxes][time_ind,:,:,:])),
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
    ax2 = Axis(fig_VDF[1, II], title=title, aspect=1, yticks=ticks, xticks=ticks, alignmode=Inside())
    hidexdecorations!(ax2, grid=false, ticks=false)
    ax3 = Axis(fig_VDF[2, II], aspect=1, yticks=ticks, xticks=ticks, alignmode=Inside())
    linkxaxes!(ax2, ax3)
    ax2, rot2, ind_mask = MAVEN_plot.VDF_2d_slip(
        ax2,v,dF; 
        normal_vectors = [[1,0,0],[0,1,0]],  
        # normal_vectors=MAVEN_STATIC.rotate_vector_with_Matrix.([[1,0,0],[0,1,0]],rotation_Q), 
        colorrange=colorrange, vbluk=vel, plot_range=plot_range, 
        vsc=vsc,magf=b0,return_rot_matrix=true, angle_range=[-30, 30])
    local masked_range_data = Dict(
        :shape => ones(size(swi_data[:diff_en_fluxes][time_ind,:,:,:])),
        :mode => "swi",
        :m_int => spice[:m_int],
    )
    for i in keys_range
        masked_range_data[i] = range_data[i][ind_mask]
    end
    ax2 = MAVEN_plot.VDF_2d_mask(ax2,masked_range_data;rot = rot2, vsc=vsc)

    ax3, rot3, ind_mask = MAVEN_plot.VDF_2d_slip(
        ax3,v,dF; 
        normal_vectors = [[1,0,0],[0,0,1]], 
        # normal_vectors=MAVEN_STATIC.rotate_vector_with_Matrix.([[1,0,0],[0,0,1]],rotation_Q), 
        colorrange=colorrange, vbluk=vel, plot_range=plot_range, 
        vsc=vsc,magf=b0, return_rot_matrix=true, angle_range=[-30, 30],xlabel=X_vdf_label)
    local masked_range_data = Dict(
        :shape => ones(size(swi_data[:diff_en_fluxes][time_ind,:,:,:])),
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
        # arrows!(ax, [Point2f(0, 0)],[Vec2f(vel_local[1],vel_local[2])]; 
        #     color=:blue,arrowsize=15,arrowtail=15,linewidth=3,align = :endpoint) # 速度矢量arrowtail
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

Colorbar(fig_VDF[1:2, 5], limits=colorrange, label=rich("i",superscript("+")), 
    colormap=:viridis, scale=log10)

save(save_path, fig_VDF)