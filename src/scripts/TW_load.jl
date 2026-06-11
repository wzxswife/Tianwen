# 读取和计算天问一号数据
# 使用data_get_from_date函数做全局读取

module TW_load
root_path = "E:/Tianwen-1/"
using Dates
include("TCWavelet.jl")
using Wavelets

export load_mag_2c_bydlm, load_minpa 
export  calculate_mag, get_spc2mso_rot_matrix_via_2c
export read_list, caculate_wavelet, find_avail_data

# function load_mag(file::String);#弃用，以dlm方法为准

#     lines=readlines(file)
#     lines = lines[20:end]

#     nums = length(lines)
#     times = Vector{DateTime}(undef,nums)
#     B = Matrix{Float32}(undef,nums,3)
#     position = Matrix{Float32}(undef,nums,3)
#     flag = Vector{Int32}(undef,nums)
    
#     @inbounds for (i,line) in enumerate(lines)
#         times[i],B[i,:],position[i,:],flag[i] = get_data_from_line(line)
#     end

#     B_total = sqrt.(sum(B.^2, dims=2)); B_total = B_total[:,1]

#     data=Dict(
#         "data staute" => true,
#         "Var name" => "time[Ntime], B_total[Ntime],B[Ntime,3],position[Ntime,3]",
#         "Vars"     => [times,B_total,B,position],
#         "flag"     => flag
#     )
#     # data = MGS_mag(true,times,B_total,B,position)
#     return data
# end
# function get_data_from_line(line::String)
#     epoch = DateTime(line[1:23], DateFormat("y-m-dTH:M:S.s"))
#     # data_t = [line[33:44],line[45:56],line[57:68],line[69:79],line[80:90],line[91:101]]
#     bx = parse(Float32,line[33:44])
#     by = parse(Float32,line[45:56])
#     bz = parse(Float32,line[57:68])
#     x = parse(Float32,line[69:79])
#     y = parse(Float32,line[80:90])
#     z = parse(Float32,line[91:101])
#     iflag = parse(Int32,line[139:140])
#     return epoch,[bx,by,bz],[x,y,z],iflag
# end
using DelimitedFiles,DataFrames
function load_mag_2c_bydlm(file::String);
    mag2c1hz = identity.(DataFrame(readdlm(file, skipstart=19), :auto))
    name = ["Time", "Sampling_Rate", "X_MSO", "Y_MSO", "Z_MSO", "Probe_Position_X_MSO", "Probe_Position_Y_MSO", "Probe_Position_Z_MSO", "Roll", "Pitch", "Yaw",  "Quality_Flags"]
    rename!(mag2c1hz, name)
    mag2c1hz[!, :Time] = map(x->DateTime(x[begin:end-4], DateFormat("y-m-dTH:M:S.s")), mag2c1hz[!, :Time])
    unique!(mag2c1hz) #remove depulicate rows
    sort!(mag2c1hz) #sorting

    times = mag2c1hz[!, :Time]
    julUTtimes = datetime2julian.(times)
    BMSO1hz = mag2c1hz[:, [:X_MSO, :Y_MSO, :Z_MSO]]
    position = mag2c1hz[:, [:Probe_Position_X_MSO, :Probe_Position_Y_MSO, :Probe_Position_Z_MSO]]
    roll = mag2c1hz[!, :Roll]
    pitch = mag2c1hz[!, :Pitch]
    yaw = mag2c1hz[!, :Yaw]
    flag = mag2c1hz[:, :Quality_Flags]
    # ind0 = findall(x-> !isnan(x), BMSO1hz[!, 1] )
    # ind1 = findall(x-> isnan(x), BMSO1hz[!, 1] )
    # if length(ind1) >= 1 && length(ind0) > 1
    #     for ib in 1:3 
    #         interp_linear = linear_interpolation(magut1hz[ind0], BMSO1hz[ind0, ib]; extrapolation_bc=Line())
    #         BMSO1hz[ind1, ib] = interp_linear(magut1hz[ind1])
    #     end
    # end

    B_total = sqrt.(BMSO1hz[:,:X_MSO].^2 + BMSO1hz[:,:Y_MSO].^2 + BMSO1hz[:,:Z_MSO].^2); B_total = B_total[:,1]

    data=Dict{Symbol,Any}(
        :data_load_flag => true,
        :time => times |> Array,
        :JulUTtime => julUTtimes |> Array,
        :epoch => times |> Array,
        :B => BMSO1hz |> Array,
        :B_total => B_total |> Array,
        :position => position |> Array,
        :flag     => flag |> Array,
        :yaw => yaw |> Array,
        :pitch => pitch |> Array,
        :roll => roll |> Array,
    )
    return data
end

function calculate_mag(position;models=["alt"])
    function c_alt(position) 
        alt = sqrt.(sum(position.^2, dims=2)) .- 3393.5
        alt = alt[:,1]
        return alt
    end
    function c_local_time(position) 
        local_time = atan.(position[:, 2], position[:, 1]) ./ π .* 12.0 .+ 12.0
        return local_time
    end
    function c_latitude(position) 
        latitude   = atan.(position[:, 3], sqrt.(sum(position[:,1:2].^2, dims=2)) ) ./ π .* 180.0
        return latitude
    end
    results = Dict{Symbol,Any}(
        :alt => c_alt,
        :local_time => c_local_time,
        :latitude => c_latitude,
    )
    for model in models
        data[model] = results[model](position)
    end
end

function get_spc2mso_rot_matrix_via_2c(dat) # 需要load_mag_2c_bydlm的返回值\
    ntime = length(dat[:time])
    rot_matrix = zeros(ntime,3,3)

    local pitch = dat[:pitch].|>deg2rad
    local yaw   = dat[:yaw]  .|>deg2rad
    local roll  = dat[:roll] .|>deg2rad
    rot_matrix = zeros(Float64, ntime, 3, 3)
    for i in 1:ntime
        local r,p,y = roll[i],pitch[i],yaw[i]
        local cp = cos(p)
        local sp = sin(p)
        local cy = cos(y)
        local sy = sin(y)
        local cr = cos(r)
        local sr = sin(r)
        rot_matrix[i,1,:] = [cp*cy, sy*cr+sr*sp*cy, sr*sy-cr*sp*cy]
        rot_matrix[i,2,:] =  [-cp*sy, cr*cy-sr*sp*sy, sr*cy+cr*sp*sy]
        rot_matrix[i,3,:] = [sp, -sr*cp, cr*cp]
    end
    return dat[:time],rot_matrix
end
function read_list()
    list_path = root_path*"List/"*"mag_list.txt"
    files = readlines(list_path)
    files = root_path.*files  
    return files
end

"""
计算小波变换（测试通过）
"""
function caculate_wavelet(mag_data, dt; mother="MORLET")
    ns =  length(mag_data[:, 1])
    wave, period, scale, coi = wavelet(reshape(mag_data[:, 1], ns), dt; pad=1, mother=mother)
    xpower = abs.(wave).^2
    wave, period, scale, coi = wavelet(reshape(mag_data[:, 2], ns), dt; pad=1, mother=mother)
    ypower = abs.(wave).^2
    wave, period, scale, coi = wavelet(reshape(mag_data[:, 3], ns), dt; pad=1, mother=mother)
    zpower = abs.(wave).^2
    Bpower =xpower+ypower+zpower
    return Bpower, period
end

"""
从数据字典中筛选出在指定时间范围内的数据（测试通过）
"""
function find_avail_data(data::Dict, time_range::Vector{DateTime}, keys)
    avail_data = Dict()
    ind = findall(minimum(time_range) .<= data[:epoch] .<= maximum(time_range))
    
    # 筛选 position 中的 NaN
    if haskey(data, :position)
        nan_rows = any(isnan.(data[:position][ind, :]), dims=2)[:]
        ind = ind[.!nan_rows]
    end
    
    # 筛选 B 中的 NaN
    if haskey(data, :B)
        nan_rows = any(isnan.(data[:B][ind, :]), dims=2)[:]
        ind = ind[.!nan_rows]
    end
    
    # 筛选 B_total 中的 NaN
    if haskey(data, :B_total)
        nan_rows = isnan.(data[:B_total][ind])
        ind = ind[.!nan_rows]
    end
    
    for key in keys
        haskey(data, key) || continue
        value = data[key]
        length(value) == 1 && continue
        avail_data[key] = ndims(value) == 2 ? value[ind, :] : value[ind]
    end
    return avail_data
end

const energymod1 = [2.81, 3.548928, 4.482167, 5.660813, 7.149401, 9.029433, 11.40385,
    14.40264, 18.19001, 22.97332, 29.01447, 36.64422, 46.28032, 58.45036, 73.82067,
    93.23282, 117.7497, 148.7135, 187.8198, 237.2095, 299.587, 378.3675, 477.8644,
    603.5253, 762.2309, 962.6694, 1215.816, 1535.532, 1939.321, 2449.292, 3093.366,
    3906.809, 4934.157, 6231.661, 7870.361, 9939.98, 12553.83, 15855.03, 20024.33, 25290.0]
const energymod4 = [2.81, 3.246924, 3.751784, 4.335145, 5.009211, 5.788088, 6.688071, 7.727991,
    8.929607, 10.31806, 11.9224, 13.77621, 15.91825, 18.39336, 21.25332, 24.55798,
    28.37647, 32.7887, 37.88697, 43.77797, 50.58496, 58.45036, 67.53873, 78.04025,
    90.17464, 104.1958, 120.3971, 139.1175, 160.7487, 185.7433, 214.6243, 247.996,
    286.5566, 331.113, 382.5974, 442.087, 510.8266, 590.2544, 682.0324, 788.0808,
    910.6185, 1052.21, 1215.816, 1404.862, 1623.303, 1875.708, 2167.36, 2504.36,
    2893.76, 3343.707, 3863.617, 4464.366, 5158.525, 5960.618, 6887.428, 7958.346, 
    9195.78, 10625.62, 12277.79, 14186.84, 16392.74, 18941.63, 21886.84, 25290.0]
const θmod1 = [11.3, 33.8, 56.2, 78.7]
const ϕmod1 = [11.25, 33.75, 56.25, 78.75, 101.25, 123.75, 146.25, 168.75,
               191.25, 213.75, 236.25, 258.75, 281.25, 303.75, 326.25, 348.75]
const massmod1 = [1, 2, 4, 16, 28, 32, 44, 64]
const TW_minpa2sc = [0.0 0.0 1.0; -1.0 0.0 0.0; 0.0 -1.0 0.0]

function switch_mod(mod)
    if mod == 1
        mass = massmod1
        minpaphi = 180.0 .- ϕmod1
        minpatheta = -θmod1
        energy = energymod1
    elseif mod == 4
        mass = massmod1
        minpaphi = 180.0 .- ϕmod1
        minpatheta = -θmod1
        energy = energymod4
    end
    return mass, energy, minpaphi, minpatheta
end

function load_minpa(minpa_file, mod)
    mass, energy, minpaphi, minpatheta = switch_mod(mod)
    nmod = length(mass) * length(minpaphi) * length(minpatheta) * length(energy)
    minpa = identity.(DataFrame(readdlm(minpa_file), :auto))
    minpaut = map(x->DateTime(x[begin:end-4], DateFormat("y-m-dTH:M:S.s")),
                  minpa[!, 1])
    minpajul = datetime2julian.(minpaut)
    minpadfi = reshape(Array(minpa[:, 60:59+nmod]),
                       length(minpaut), length(mass), length(minpaphi),
                       length(minpatheta), length(energy))
    minpadfi[findall(x->x<0, minpadfi)] .= 0.0
    minpaeflux = minpadfi .* reshape(energy, 1, 1, 1, 1, length(energy))
    println("  MINPA data: $(length(minpaut)) time points")
    data=Dict{Symbol,Any}(
        :epoch => minpaut,
        :JulUTtime => minpajul,
        :mod => mod,
        :mass => mass,
        :energy => energy,
        :phi => minpaphi,
        :theta => minpatheta,
        # :minpadfi => minpadfi,
        :eflux => minpaeflux
    )
    return data
end

# -----------以下为DataFrame版本的读取函数（已弃用）----------------
# function load_mag_TW_DF(file::String);
#     println("Reading file: ")
#     println(file)
#     local Rm = 3393.5 #火星半径，单位km
#     mag2c32hz = identity.(DataFrame(readdlm(file, skipstart=19), :auto))
#     name = ["Time", "Sampling_Rate", "X_MSO", "Y_MSO", "Z_MSO", "Probe_Position_X_MSO", "Probe_Position_Y_MSO", "Probe_Position_Z_MSO", "Roll", "Pitch", "Yaw",  "Quality_Flags"]
#     rename!(mag2c32hz, name)
#     mag2c32hz[!, :Time] = map(x->DateTime(x[begin:end-4], DateFormat("y-m-dTH:M:S.s")), mag2c32hz[!, :Time])
#     mag2c32hz[!, :JulUT] = datetime2julian.(mag2c32hz[!, :Time])
#     mag2c32hz[!,[ :Probe_Position_X_MSO, :Probe_Position_Y_MSO, :Probe_Position_Z_MSO]] = mag2c32hz[!,[ :Probe_Position_X_MSO, :Probe_Position_Y_MSO, :Probe_Position_Z_MSO]]./Rm
#     unique!(mag2c32hz) #remove depulicate rows
#     sort!(mag2c32hz) #sorting
#     magut32hz = mag2c32hz[:, :JulUT].-mag2c32hz[1, :JulUT]
#     BMSO32hz = mag2c32hz[:, [:X_MSO, :Y_MSO, :Z_MSO]]
#     ind0 = findall(x-> !isnan(x), BMSO32hz[!, 1] ) 
#     ind1 = findall(x-> isnan(x), BMSO32hz[!, 1] ) 
#     if length(ind1) >= 1 && length(ind0) > 1
#         for ib in 1:3 
#             interp_linear = linear_interpolation(magut32hz[ind0], BMSO32hz[ind0, ib]; extrapolation_bc=Line())
#             BMSO32hz[ind1, ib] = interp_linear(magut32hz[ind1])
#         end
#     end
#     println("Data read finished!")
#     return mag2c32hz, BMSO32hz
# end

# function load_mag_MAVEN_DF(file::String, year);
#     println("Reading file: ")
#     println(file)
#     skip = 1 
#     for line in eachline(file)
#         if line[1:6] != "  $year"
#             skip = skip+1
#         else
#             break
#         end
#     end

#     local Rm = 3393.5 #火星半径，单位km

#     mag2c32hz_temp = identity.(DataFrame(readdlm(file, skipstart=skip), :auto))

#     BMSO32hz = Array(mag2c32hz_temp[:, 8:10])
#     magut32hz = @. DateTime(mag2c32hz_temp[:, 1]) + Day(mag2c32hz_temp[:, 2]-1) + Hour(mag2c32hz_temp[:, 3]) + Minute(mag2c32hz_temp[:, 4]) + Second(mag2c32hz_temp[:, 5]) + Millisecond(mag2c32hz_temp[:, 6])
#     magjlut32hz = datetime2julian.(magut32hz)
#     magjlut32hz= magjlut32hz.-magjlut32hz[1]
#     PosMSO32hz = Array(mag2c32hz_temp[:, 12:14])./Rm
#     mag2c32hz = DataFrame()
#     mag2c32hz.Time = magut32hz
#     mag2c32hz.JulUT = magjlut32hz
#     mag2c32hz.X_MSO = BMSO32hz[:, 1]
#     mag2c32hz.Y_MSO = BMSO32hz[:, 2]
#     mag2c32hz.Z_MSO = BMSO32hz[:, 3]
#     mag2c32hz.Probe_Position_X_MSO = PosMSO32hz[:, 1]
#     mag2c32hz.Probe_Position_Y_MSO = PosMSO32hz[:, 2]
#     mag2c32hz.Probe_Position_Z_MSO = PosMSO32hz[:, 3]
#     println("Data read finished!")
#     return mag2c32hz, BMSO32hz, PosMSO32hz
# end


end