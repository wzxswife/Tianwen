module MVN3DIonDisttribution

# Searching File
using Glob
# Reading File
using CDFDatasets
using DelimitedFiles
using DataFrames
# Stucturing Data
using Dates
using Unitful
# Computing
using StaticArrays
using LinearAlgebra
using Statistics
using Rotations
using Polyester # parallel computing
using LoopVectorization # parallel loops
using Tullio  # Einstein summation
using AngleBetweenVectors
include("PhysicsCoordTransforms.jl") # coordinate transformation
include("SphericalSplinesOptim.jl") # Spheric Spline Interpolation
using .PhysicsCoordTransforms # coordinate transformation
using .SphericalSplinesOptim # Spheric Spline Interpolation 
using WENO4 # WENO interpolation
using BenchmarkTools
# Plotting Data
using CairoMakie # 2D
using GLMakie # 3D

# constants
import PhysicalConstants.CODATA2018: m_p, m_e, e, ε_0, μ_0, c_0, G, k_B
import UnitfulAstro: AU
using AstroLib
Mars = AstroLib.planets["mars"]
M_Mars = Mars.mass * 1u"kg"
R_Mars = Mars.radius * 1u"m"
#===========================================================================
Planet(name, radius, eqradius, polradius, mass, ecc, axis, period)
AstroLib.planets["mars"]
julia> AstroLib.planets["mars"]
Planet:                      Mars
mean radius:                 3.3895e6 m
equatorial radius:           3.3962e6 m
polar radius:                3.3762e6 m
mass:                        6.4169e23 kg
eccentricity:                0.0933941
semi-major axis:             2.2794382242757306e11 m
period:                      5.9355072e7 s
inclination:                 1.84969142 °
longitude of ascending node: 49.55953891 °
longitude of perihelion:     -23.94362959 °
mean longitude:              -4.55343205 °
============================================================================#

unitvec(v::AbstractVector{<:AbstractFloat}) = v ./ norm(v)

"""
	Reading MVN MAG data with different coordinates and time resolutions

	coordinate::"ss" for MSO coordinate system，and "pc" for geographic coordinate system
	resolution::"1s" for data with 1 Hz sampling rate, and "" for data with 64 Hz samping rate

	The data included magetic field components and spacecraft positions are exported into a DataFrame.
"""
function read_mvn_mag(
	date_str::String;
	data_root::String = "E:/Data", 
	coordinate::String = "ss", # or "pc"
	resolution::String = "1s",  # or "" representing data with 64 Hz samping rate
	Rm = ustrip(u"km", R_Mars)  # Mars' radius in units of km
)
	folder = "$(data_root)/MAG/l2/$(date_str[1:4])/$(date_str[5:6])/"	
	file = glob("*$coordinate$resolution*_*$date_str*.sts", folder)
	
	if isempty(file)
		error("No matching file found: folder=$folder, pattern=*$coordinate$resolution*_*$date_str*.sts")
	elseif length(file) > 1
		@warn "Multiple files matched, using the first one" folder=folder matched=file
	end

	file = file[1]
	lines = readlines(file)
	skip = findfirst(line -> occursin(r"^\s*\d{4}\s+\d{3}", line), lines)
	data = identity.(DataFrame(readdlm(IOBuffer(join(lines[skip:end], "\n"))), :auto))
	name = [:year,  :doy,  :hour,  :min,  :sec,  :msec,  :dddy,  :bx,  :by,  :bz,  :brange,  :px,  :py,  :pz,  :plbx,  :plby,  :plbz,  :plrange]
	rename!(data, name)
	transform!(data,
    [:year, :doy, :hour, :min, :sec, :msec] => ByRow(
        (y, d, h, m, s, ms) ->
            DateTime(y) + Day(d - 1) + Hour(h) +
            Minute(m) + Second(s) + Millisecond(ms)
    ) => :epoch)
	transform!(data,
    :epoch => ByRow(time -> datetime2julian(time)) => :julian)
	transform!(data, [:px, :py, :pz] .=> (x -> x ./ Rm) .=> [:px, :py, :pz]) 
	transform!(data, [:px, :py, :pz] => ByRow((x, y, z) -> sqrt(x^2+y^2+z^2)) => :pr) 
	select!(data, :epoch, :julian, :bx, :by, :bz, :px, :py, :pz, :pr)

	return data
end

"""
	Reading MVN KP data

	The parameters and data are exported as two DataFrames. Through searching the "param" DataFrame, 
	one can identify which colomun of "data" DataFrame contains the required data. For Example:
		filter(row -> row.instrument == "MAG", param)
		filter(row -> row.unit == "cm^-3", param)
		filter(row -> occursin(r"density"i, row.parameter), param)
		filter(row -> occursin(r"periapsis"i, row.note), param)
		subset(param, 
			:instrument => ByRow(==("STATIC")),
			:parameter  => ByRow(p -> occursin(r"temperature"i, p))
		)
	Once the related coloumn number "n" is identified, one can access the data through data[:, n]. 
"""
function read_mvn_kp(
	date_str::String; 
	data_root::String = "E:/Data/MAVEN"
)
	folder = "$data_root/KP/insitu/$(date_str[1:4])/$(date_str[5:6])/"
    file = glob("*$date_str*.tab", folder)

	if isempty(file)
		error("No matching file found: folder=$folder, pattern=*$coordinate$resolution*_*$date_str*.sts")
	elseif length(file) > 1
		@warn "Multiple files matched, using the first one" folder=folder matched=file
	end

	file = file[1]
	lines = readlines(file)

    data_start_line = parse(Int, match(r"(\d+)\s+Line on which data begins",
        lines[findfirst(l -> occursin("Line on which data begins", l), lines)])[1]
	)
    tbl_start = findfirst(l -> startswith(l, "# PARAMETER"), lines) + 2
    tbl_end = findlast(l -> startswith(l, "#") && isempty(strip(l[2:end])), lines) - 2

    param = DataFrame(
        parameter   = strip.(SubString.(lines[tbl_start:tbl_end], 2, 60)),
        instrument  = strip.(SubString.(lines[tbl_start:tbl_end], 61, 72)),
        unit        = strip.(SubString.(lines[tbl_start:tbl_end], 73, 90)),
        column      = strip.(SubString.(lines[tbl_start:tbl_end], 91, 100)),
        format      = strip.(SubString.(lines[tbl_start:tbl_end], 101, 122)),
        note        = strip.(SubString.(lines[tbl_start:tbl_end], 123, 270))
    )

	data = identity.(DataFrame(readdlm(IOBuffer(join(lines[data_start_line:end], "\n"))), :auto))
	transform!(data,
    :x1 => ByRow(time -> DateTime(time, "yyyy-mm-ddTHH:MM:SS")) => :x1)
	transform!(data,
    :x1 => ByRow(time -> datetime2julian(time)) => :julian)

    return param, data
end

"""
	Reading MVN SWIA coarsesvy3d data
	
	The deflection angle ("theta_coarse") grids is energy-dependent for the first 12 energy channels, but consistent for the rest energy channels. 
"""
function read_mvn_swia_coarsesvy3d(
	date_str::String,
	kp::DataFrame;
	data_root::String = "E:/Data/MAVEN"
)
	folder = "$data_root/SWIA/l2/$(date_str[1:4])/$(date_str[5:6])/"
	file = glob("*coarsesvy3d_*$date_str*.cdf", folder)
	
	if isempty(file)
		error("No matching file found: folder=$folder, pattern=*coarsesvy3d_*$date_str*_v02_r01.cdf")
	elseif length(file) > 1
		@warn "Multiple files matched, using the first one" folder=folder matched=file
	end
	datacdf = cdfopen(file)
	def = sanitize(datacdf["diff_en_fluxes"]) # differential energy flux phi_coarse × dindex × energy_coarse × epoch
	epoch = DateTime.(Array(datacdf["epoch"]))
	julian = datetime2julian.(epoch)
	energy = dropdims(sanitize(datacdf["energy_coarse"]), dims=(2, ))  # energy
	speed = sqrt.(energy .* Float32(2.0/ustrip(u"eV", m_p*c_0^2)) ) .* Float32(ustrip(u"km/s", c_0)) # scalar speed
	phi = dropdims(sanitize(datacdf["phi_coarse"]), dims=(2, )) # longitude deg
	psi = dropdims(sanitize(datacdf["theta_coarse"]), dims=(3, )) # latitude deg
	theta = 90.0f0 .- psi # zenith angle
	denergy = sanitize(datacdf["de_over_e_coarse"]) .* energy # energy width
	domega = stack(
		map(eachcol(psi)) do col
			#cosd.(col) .* deg2rad(col[2]-col[1]) .* deg2rad(phi[2]-phi[1])
			2.0f0 * sind((col[2]-col[1])/2.0f0) * deg2rad(phi[2]-phi[1]) .* cosd.(col) # more accurate than the previous one
		end
	)# solid angle width
	
	sc_from_swia = [0.0f0 -1.0f0  0.0f0; 
					0.0f0  0.0f0 -1.0f0; 
					1.0f0  0.0f0  0.0f0
	] # from swia to spacecraft (sc) transformation matrix
	mso_from_sc = let j = Ref(1), m = size(kp, 1)
		map(epoch) do t
			while j[] < m && abs(kp[j[]+1, 1] - t) <= abs(kp[j[], 1] - t)
				j[] += 1
			end
			transpose(reshape(Array(kp[j[], 227:235]), 3, 3)) .|> Float32
		end
	end # from sc to mso transformation matrix
	mso_from_swia = stack(mso_from_sc .* Ref(sc_from_swia))

	return (; def, epoch, julian, energy, speed, phi, psi, theta, denergy, domega, mso_from_swia)
end

"""
	Calculating MVN SWIA ion moments: density, bulk velocity, and temperature

	mvn_swia: MVN SWIA data from read_mvn_swia_coarsesvy3d
	xd: time interval for computing the moments

	average = true: averaged over the time interval
	average = false: time-dependent moments

	coordinate = "swia" or "mso"

	Return:
	(; density, bulk_velocity, temperature)
	density = ∑ def / speed / energy * denergy * domega
	bulkvelocity = (∑ def / speed * velocity / energy * denergy * domega) / density
	pressure = m (∑ def / speed / energy * denergy * domega * (velocity - bulkvelocity) ⊗ (velocity - bulkvelocity)) 
	temperature = tr(pressure) / 3 / density 
"""
function calculate_mvn_swia_ion_moments(
	mvn_swia,
	xd;
	average = false,
	coordinate = "mso"
)
	densityunit = Float32(ustrip(u"cm^-3", 1u"cm^-2*s^-1" / 1u"km/s"))
	temperatureunit = Float32(ustrip(u"eV", m_p * 1u"km/s"^2))
	phi = mvn_swia.phi
	theta = mvn_swia.theta
	speed = mvn_swia.speed
	nphi, ntheta, nspeed, nepoch = size(mvn_swia.def)
	#ind = findall(xd[begin] .≤ mvn_swia.epoch .≤ xd[end])
	ind1 = searchsortedfirst(mvn_swia.epoch, xd[begin])  
	ind2 = searchsortedlast(mvn_swia.epoch, xd[end])
	ind = ind1:ind2
	velocity = Array{Float32, 4}(undef, 3, size(mvn_swia.def)[1:3]...)
	@batch for iphi in 1:nphi
		sp, cp = sincosd(phi[iphi])
		for ispeed in 1:nspeed
			for itheta in 1:ntheta
				st, ct = sincosd(theta[itheta, ispeed])
				velocity[1, iphi, itheta, ispeed] = st * cp * speed[ispeed]
				velocity[2, iphi, itheta, ispeed] = st * sp * speed[ispeed]
				velocity[3, iphi, itheta, ispeed] = ct * speed[ispeed]
			end
		end
	end

	@tullio weight[itheta, ispeed] := mvn_swia.domega[itheta, ispeed] * mvn_swia.denergy[ispeed] / mvn_swia.energy[ispeed] / speed[ispeed]
	if average
		def = dropdims(mean(mvn_swia.def[:, :, :, ind], dims = 4), dims = 4)
		@tullio density := weight[itheta, ispeed] * def[iphi, itheta, ispeed]
		@tullio bulkvelocity[ic] := weight[itheta, ispeed] * def[iphi, itheta, ispeed] * velocity[ic, iphi, itheta, ispeed] / $density # in units of km/s
		@tullio temperature := weight[itheta, ispeed] * def[iphi, itheta, ispeed] * (velocity[ic, iphi, itheta, ispeed] - bulkvelocity[ic])^2 / $density / 3
		density *= densityunit
		temperature *= temperatureunit	
		indx = argmin(abs.(mvn_swia.epoch .- ((xd[end] - xd[begin])/2 + xd[begin])) )
		mso_from_swia = mvn_swia.mso_from_swia[:, :, indx]
		epoch = mvn_swia.epoch[indx]
		if coordinate == "mso"
			bulkvelocity = mso_from_swia * bulkvelocity
		end
	else
		def = view(mvn_swia.def, :, :, :, ind)
		@tullio density[iepoch] := weight[itheta, ispeed] * def[iphi, itheta, ispeed, iepoch]
		@tullio bulkvelocity[ic, iepoch] := weight[itheta, ispeed] * def[iphi, itheta, ispeed, iepoch] * velocity[ic, iphi, itheta, ispeed] / density[iepoch] # in units of km/s
		@tullio temperature[iepoch] := weight[itheta, ispeed] * def[iphi, itheta, ispeed, iepoch] * (velocity[ic, iphi, itheta, ispeed] - bulkvelocity[ic, iepoch])^2 / density[iepoch] / 3

		density .*= densityunit
		temperature .*= temperatureunit	
		mso_from_swia = mvn_swia.mso_from_swia[:, :, ind]
		epoch = mvn_swia.epoch[ind]
		if coordinate == "mso"
			@batch for iepoch in eachindex(ind)
				bulkvelocity[:, iepoch] = view(mso_from_swia, :, :, iepoch) * view(bulkvelocity, :, iepoch)
			end
		end
	end

	return (; epoch, density, bulkvelocity, temperature, mso_from_swia, phi, theta, speed, velocity, def)
end

"""
	Reading MVN STATIC 32e4d16a8m data

	"swp_ind" index (0-based indexing, < NSWP) identifies the energy and deflector sweep look up tables (LUT) for the sensor. SWP_IND is an index that selects the following support data arrays: ENERGY, DENERGY, THETA, DTHETA, PHI, DPHI, DOMEGA, GF and MASS_ARR.
	
	resolution = "d0" or "d1".
"""
function read_mvn_static_32e4d16a8m(
	date_str::String;
	data_root::String = "E:/Data/MAVEN",
	resolution::String = "d0"
)
	folder = "$data_root/STATIC/l2/$(date_str[1:4])/$(date_str[5:6])/"
	file = glob("*$resolution*32e4d16a8m_*$date_str*.cdf", folder)
	
	if isempty(file)
		error("No matching file found: folder=$folder, pattern=*32e4d16a8m_*$date_str*_v02_r07.cdf")
	elseif length(file) > 1
		@warn "Multiple files matched, using the first one" folder=folder matched=file
	end
	datacdf = cdfopen(file)

	swp = Array(datacdf["swp_ind"]) .+ 1 # Julia is 1-based indexing while the original data is 0-based indexing
	epoch = DateTime.(Array(datacdf["epoch"]))
	julian = datetime2julian.(epoch)
	def = PermutedDimsArray(reshape(sanitize(datacdf["eflux"]), 8, 4, 16, 32, length(epoch)), [1, 3, 2, 4, 5]) # differential energy flux compno_8 × compno_64 × compno_32 × epoch
	energy = PermutedDimsArray(reshape(dropdims(sanitize(datacdf["energy"]), dims = (5, )), 8, 4, 16, 32, 28), [1, 3, 2, 4, 5])  # energy 	(8×4×16×32×28)
	mass = PermutedDimsArray(reshape(dropdims(sanitize(datacdf["mass_arr"]), dims = (5, )), 8, 4, 16, 32, 28), [1, 3, 2, 4, 5]) # mass (8×4×16×32×28)
	psi = PermutedDimsArray(reshape(dropdims(sanitize(datacdf["theta"]), dims = (5, )), 8, 4, 16, 32, 28), [1, 3, 2, 4, 5]) # latitude deg (8×4×16×32×28)
	theta = 90.0 .- psi # zenith angle deg (8×4×16×32×28)
	phi = PermutedDimsArray(reshape(dropdims(sanitize(datacdf["phi"]), dims = (5, )), 8, 4, 16, 32, 28), [1, 3, 2, 4, 5]) # longitude deg (8×4×16×32×28)
	speed = sqrt.(energy .* (2.0/ustrip(u"eV", m_p*c_0^2)) ./ mass).*ustrip(u"km/s", c_0) # scalar speed
	denergy = PermutedDimsArray(reshape(dropdims(sanitize(datacdf["denergy"]), dims = (5, )), 8, 4, 16, 32, 28), [1, 3, 2, 4, 5]) # energy width (8×4×16×32×28)
	#domega = PermutedDimsArray(reshape(dropdims(sanitize(datacdf["domega"]), dims = (5, )), 8, 4, 16, 32, 28), [1, 3, 2, 4, 5]) # solid angle (8×4×16×32×28)
	domega = similar(denergy) # more accurate than the one above, taking into account the integration over solid angle
	@batch for iswp in 1:28
		for ienergy in 1:32
			for itheta in 1:4
				for iphi in 1:16
					for imass in 1:8
						domega[imass, iphi, itheta, ienergy, iswp] = abs(2.0 * sind((theta[imass, iphi, 2, ienergy, iswp]-theta[imass, iphi, 1, ienergy, iswp])/2.0) * deg2rad(phi[imass, 2, itheta, ienergy, iswp]-phi[imass, 1, itheta, ienergy, iswp]) * sind(theta[imass, iphi, itheta, ienergy, iswp])) 
					end
				end
			end
		end
	end

	# Doc about quaternion
	# https://danceswithcode.net/engineeringnotes/quaternions/quaternions.html
	quat_mso = sanitize(datacdf["quat_mso"]) # Quaternion elements to rotate from STATIC coordinates (same as APP coordinates) to MSO coordinates (NUM_DISTS, 3)
	mso_from_static = stack([QuatRotation(quat) for quat in eachcol(quat_mso)]) # active (right handed) rotations where applicable.

	return (; def, swp, epoch, julian, mass, energy, speed, phi, psi, theta, denergy, domega, mso_from_static)
end

"""
	Calculating MVN STATIC ion moments: density, bulk velocity, and temperature

	mvn_static: MVN static data from read_mvn_static_32e4d16a8m
	xd: time interval for computing the moments
	imass: ion species
	average: if true, average the moments over the time interval
	coordinate: "static" or "mso"
"""
function calculate_mvn_static_ion_moments(
	mvn_static,
	xd,
	imass;
	average = false,
	coordinate = "mso"
)
	densityunit = Float32(ustrip(u"cm^-3", 1u"cm^-2*s^-1" / 1u"km/s"))
	temperatureunit = Float32(ustrip(u"eV", m_p * 1u"km/s"^2))

	nmass, nphi, ntheta, nspeed, nswp = size(mvn_static.energy)
	#ind = findall(xd[begin] .≤ mvn_static.epoch .≤ xd[end])
	ind1 = searchsortedfirst(mvn_static.epoch, xd[begin])  
	ind2 = searchsortedlast(mvn_static.epoch, xd[end])
	ind = ind1:ind2

	if average
		swpu = unique(mvn_static.swp[ind])
		if  length(swpu) > 1
			error("More than one swp found during the specified time interval. Please specify a time interval that corresponds to a single swp.")
		end
		swp = swpu[1]
		def = dropdims(mean(view(mvn_static.def, imass, :, :, :, ind), dims=4), dims=4)
		speed = view(mvn_static.speed, imass, :, :, :, swp)
		theta = view(mvn_static.theta, imass, :, :, :, swp)
		phi = view(mvn_static.phi, imass, :, :, :, swp)
		domega = view(mvn_static.domega, imass, :, :, :, swp)
		energy = view(mvn_static.energy, imass, :, :, :, swp)
		denergy = view(mvn_static.denergy, imass, :, :, :, swp)
		velocity = Array{Float32, 4}(undef, 3, size(speed)...)
		@batch for ispeed in 1:nspeed
			for itheta in 1:ntheta
				@simd for iphi in 1:nphi
					st, ct = sincosd(theta[iphi, itheta, ispeed])
					sp, cp = sincosd(phi[iphi, itheta, ispeed])
					velocity[1, iphi, itheta, ispeed] = st * cp * speed[iphi, itheta, ispeed]
					velocity[2, iphi, itheta, ispeed] = st * sp * speed[iphi, itheta, ispeed]
					velocity[3, iphi, itheta, ispeed] = ct * speed[iphi, itheta, ispeed]
				end
			end
		end
		weight = domega .* denergy ./ energy ./ speed
		@tullio density := weight[iphi, itheta, ispeed] * def[iphi, itheta, ispeed]
		@tullio bulkvelocity[ic] := weight[iphi, itheta, ispeed] * def[iphi, itheta, ispeed] * velocity[ic, iphi, itheta, ispeed] / $density # in units of km/s
		@tullio temperature := weight[iphi, itheta, ispeed] * def[iphi, itheta, ispeed] * (velocity[ic, iphi, itheta, ispeed] - bulkvelocity[ic])^2 / $density / 3
		density *= densityunit
		temperature *= temperatureunit	
		indx = argmin(abs.(mvn_static.epoch .- ((xd[end] - xd[begin])/2 + xd[begin])) )
		epoch = mvn_static.epoch[indx]
		mso_from_static = mvn_static.mso_from_static[:, :, indx]
		if coordinate == "mso"
			bulkvelocity = mso_from_static * bulkvelocity
		end
	else
		nepoch = length(ind)
		def = mvn_static.def[imass, :, :, :, ind]
		swp = mvn_static.swp[ind]
		speed = view(mvn_static.speed, imass, :, :, :, swp)
		theta = view(mvn_static.theta, imass, :, :, :, swp)
		phi = view(mvn_static.phi, imass, :, :, :, swp)
		domega = view(mvn_static.domega, imass, :, :, :, swp)
		energy = view(mvn_static.energy, imass, :, :, :, swp)
		denergy = view(mvn_static.denergy, imass, :, :, :, swp)
		mso_from_static = view(mvn_static.mso_from_static, :, :, ind)
		velocity = Array{Float32, 5}(undef, 3, size(speed)...)
		@batch for iepoch in 1:nepoch
			for ispeed in 1:nspeed
				for itheta in 1:ntheta
					@simd for iphi in 1:nphi
						st, ct = sincosd(theta[iphi, itheta, ispeed, iepoch])
						sp, cp = sincosd(phi[iphi, itheta, ispeed, iepoch])
						velocity[1, iphi, itheta, ispeed, iepoch] = st * cp * speed[iphi, itheta, ispeed, iepoch]
						velocity[2, iphi, itheta, ispeed, iepoch] = st * sp * speed[iphi, itheta, ispeed, iepoch]
						velocity[3, iphi, itheta, ispeed, iepoch] = ct * speed[iphi, itheta, ispeed, iepoch]
					end
				end
			end
		end
		weight = domega .* denergy ./ energy ./ speed
		@tullio density[iepoch] := weight[iphi, itheta, ispeed, iepoch] * def[iphi, itheta, ispeed, iepoch]
		@tullio bulkvelocity[ic, iepoch] := weight[iphi, itheta, ispeed, iepoch] * def[iphi, itheta, ispeed, iepoch] * velocity[ic, iphi, itheta, ispeed, iepoch] / density[iepoch] # in units of km/s
		@tullio temperature[iepoch] := weight[iphi, itheta, ispeed, iepoch] * def[iphi, itheta, ispeed, iepoch] * (velocity[ic, iphi, itheta, ispeed, iepoch] - bulkvelocity[ic, iepoch])^2 / density[iepoch] / 3
		density *= densityunit
		temperature *= temperatureunit	
		epoch = mvn_static.epoch[ind]
		if coordinate == "mso"
			@batch for iepoch in 1:nepoch
				bulkvelocity[:, iepoch] = view(mso_from_static, :, :, iepoch) * view(bulkvelocity, :, iepoch)
			end
		end
	end

	return (; epoch, density, bulkvelocity, temperature, mso_from_static, phi, theta, speed, velocity, def)
end

"""
	Plotting MVN SWIA eflux overview

	mvn_swia: MVN SWIA data from read_mvn_swia_coarsesvy3d
	xd: x date ticks for plotting
"""
function plot_mvn_swia_eflux_overview(mvn_swia, xd)
    size_inches = (20, 30)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 20, figure_padding = (5,20, 9, 10))
    cn = cgrad(:seaborn_dark, 10, categorical = true)
    lcn = Makie.wong_colors()
    xtk = datetime2julian.(xd)
    xtk = Float32.(xtk .- mvn_swia.julian[1])

    axes = [Axis(fig[j, i], yscale=log10, xticks = (xtk, Dates.format.(xd, "HH:MM")), backgroundcolor=:gray80, xticklabelrotation=-π/3) for j in 1:16, i in 1:4]

    for j in 1:16
        for i in 1:4
            ax = axes[j, i]
            heatmap!(ax, mvn_swia.julian .- mvn_swia.julian[1], mvn_swia.energy[:, 1], transpose(mvn_swia.def[j, i, :, :]), colormap=:gist_earth, colorscale=log10, colorrange = (1e5, 1e9))
            text!(ax, 0.0, 1.0, text = rich("ψ=$i ϕ=$j"), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 20, color = :white)
            xlims!(ax, xtk[begin], xtk[end]) 
            ylims!(ax, 1e1, 1e4)  
        end
    end
    linkaxes!(axes...)
    hidexdecorations!.(axes[1:15,:], grid=false)
    hideydecorations!.(axes[:, 2:end], grid=false)
    colgap!(fig.layout, 30)
    Colorbar(fig[1:16, 5], colormap=:gist_earth, scale=log10, colorrange = (1e5, 1e9), label = rich(rich("MAVEN SWIA Ion Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")))
    Label(fig[1:16, 0], rich("Energy (eV)"), fontsize = 25, halign = :center, rotation = π/2)
    Label(fig[17, 1:5], rich("Universal Time"), fontsize = 20, halign = :center)

    fig
end

"""
	Plotting MVN SWIA eflux averaged over the specified time interval

	mvn_swia: MVN SWIA data from read_mvn_swia_coarsesvy3d
	xd: time interval for averaging
"""
function plot_mvn_swia_efluxmean(mvn_swia, xd)
	ind = findall(xd[begin] .≤ mvn_swia.epoch .≤ xd[end])
	efluxmean = dropdims(mean(mvn_swia.def[:, :, :, ind], dims = 4), dims = 4)
    size_inches = (9, 18)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 20, 9, 8))
    
    cn = cgrad(:gist_earth, 18, categorical = true)

    axs = [Axis(fig[i, 1], xminorticksvisible=true, xminorticks = IntervalsBetween(9), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich("Energy (eV)"), ylabel = rich(rich("Ion Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), xscale = log10, yscale=log10, backgroundcolor=:gray80) for i in 1:4]
    for i in 1:4
        for j in 1:16
            ind = findall(efluxmean[j, i, :] .> 0)
            scatterlines!(axs[i], mvn_swia.energy[ind, 1], efluxmean[j, i, ind], color = cn[j+1], linewidth = 3)
        end
        text!(axs[i], 0.0, 1.0, text = "MAVEN SWIA ψ=$i", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
        ax1_inset = Axis(fig[i, 1], width=Relative(0.8), height=Relative(0.05), halign=0.5, valign=0.02)
        phi = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], (16, 1))
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        y = [1]
        heatmap!(ax1_inset, x, y, phi, colormap=cn, colorrange=(0,17))
        hidedecorations!(ax1_inset, grid=false)
        for ix in x
            text!(ax1_inset, ix, 1.0, text = "$ix", align = (:center, :center), offset = (0, 0), fontsize = 20, color = :black)
        end
        text!(axs[i], 0.042, 0.07, text = "φ=", align = (:left, :top), offset = (0, 0), space = :relative, fontsize = 20, color = :black)
    end
    
    xlims!.(axs, 1e0, 1e4)
    ylims!.(axs, 1e3, maximum(filter(!isnan, efluxmean)))
    linkaxes!(axs...)
    hidexdecorations!.(axs[begin:end-1], grid=false)

    fig
end

"""
	Computing intersection points between an infinite line and a sphere.

# Arguments
- `p0`: A point on the line (3D vector)
- `d`: Direction vector of the line (3D, non-zero; need not be normalized)
- `center`: Sphere center (3D vector)
- `radius`: Sphere radius (scalar, same units as p0/center/d)

# Returns
- `nothing` if no intersection
- `(Point,)` if tangent (one point)
- `(Point, Point)` if two intersections
"""
function intersect_line_sphere(p0, d, center, radius)
    oc = p0 - center
    a = dot(d, d)
    b = 2 * dot(d, oc)
    c = dot(oc, oc) - radius^2

    Δ = b^2 - 4 * a * c

    if Δ < 0
        return nothing
    else
        sqrtΔ = sqrt(Δ)
        t1 = (-b - sqrtΔ) / (2 * a)
        t2 = (-b + sqrtΔ) / (2 * a)
        return (p0 + t1 * d, p0 + t2 * d)
    end
end

"""
	Interpolating MVN SWIA eflux profile

	mvn_swia: MVN SWIA data from read_mvn_swia_coarsesvy3d
	xd: time interval for averaging
	point: a point on the line (3D vector)
	dir: direction vector of the line (3D, non-zero; need not be normalized)
	coordinate: coordinate system of the point to be interpolated
"""
function interpolate_mvn_swia_eflux_profile(mvn_swia, xd, point, dir; coordinate="swia")
	mvn_swia_ion_moments = calculate_mvn_swia_ion_moments(
		mvn_swia,
		xd;
		average = true,
		coordinate = "swia"
	)
	def = mvn_swia_ion_moments.def
	phi = mvn_swia_ion_moments.phi
	theta = mvn_swia_ion_moments.theta
	speed = mvn_swia_ion_moments.speed
	mso_from_swia = mvn_swia_ion_moments.mso_from_swia
	swia_from_mso = transpose(mso_from_swia)

	if coordinate == "swia"
		p0 = point
		d = dir / norm(dir)
	elseif coordinate == "mso"
		p0 = swia_from_mso * point
		d = swia_from_mso * dir / norm(dir)
	else
		error("Invalid coordinate system")
	end

	nphi, ntheta, nspeed = size(def)

	def_profile = fill(NaN64, nspeed * 2)
	vel_profile = fill(NaN64, nspeed * 2)

	# Spheric Splines at all the energy channels
	@batch for ispeed in 1:nspeed
		dirs = [cartesian_from_spherical()(1.0, deg2rad(theta[itheta, ispeed]), deg2rad(phi[iphi]))
		for itheta in 1:ntheta for iphi in 1:nphi] |> x -> stack(x, dims = 1)
		vals = [Float64(def[iphi, itheta, ispeed]) for itheta in 1:ntheta for iphi in 1:nphi] .|> log10
		idx = findall(isfinite, vals)
		ssp = spherical_spline(view(dirs, idx, :), view(vals, idx))

		# Effetive theta range for interpolation
		theta_range = collect(deg2rad.(extrema(theta[:, ispeed])))
		half_width = deg2rad.((theta[1, ispeed] - theta[2, ispeed]) / 2)
		theta_range[1] -= half_width
		theta_range[2] += half_width

		center = [0.0, 0.0, 0.0]
		radius = Float64(speed[ispeed])
		pinterp = intersect_line_sphere(p0, d, center, radius)
		if !isnothing(pinterp)
			for ipt in 1:2
				pt = pinterp[ipt]
				vel_profile[ispeed*2+ipt-2] = dot(pt, d)	
				vρ, vθ, vϕ = spherical_from_cartesian()(pt[1], pt[2], pt[3])
				idx = searchsortedlast(speed, vρ, rev = true)
				if 1 < idx < 48 && theta_range[1] ≤ vθ ≤ theta_range[2]
					query_dir = collect(cartesian_from_spherical()(1.0, vθ, vϕ))
					def_profile[ispeed*2+ipt-2] = exp10(ssp(query_dir))
				end
			end
		end
	end

	idx = findall(isfinite, def_profile)
	def = def_profile[idx]
	velocity = vel_profile[idx]
	idx = sortperm(velocity)
	def = def[idx]
	velocity = velocity[idx]
	return (; def, velocity)
end

"""
	Plotting MVN SWIA eflux profile parallel and perpendicular to the background magnetic field 

	mvn_swia: MVN SWIA data from read_mvn_swia_coarsesvy3d
	mvn_mag: MVN magnetic field data from read_mvn_mag
	xd: time interval for averaging

	These profiles cross the bulk velocity. The profiles are parallel or perpendicular to the background magnetic field.
"""
function plot_mvn_swia_eflux_profile(mvn_swia, mvn_mag, xd)
	# retriving MVN SWIA ion moments
	mvn_swia_ion_moments = calculate_mvn_swia_ion_moments(
		mvn_swia,
		xd;
		average = true,
		coordinate = "swia"
	)
	# bulkvelocity vector in SWIA coordinate 
	bulkvelocity = mvn_swia_ion_moments.bulkvelocity

	# tranfromation matrices
	mso_from_swia = mvn_swia_ion_moments.mso_from_swia
	swia_from_mso = transpose(mso_from_swia)

	# magnetic field unit vector in SWIA coordinate
	bvec = (mvn_mag 
		|> df -> subset(df, :epoch => ByRow(e -> xd[begin] ≤ e ≤ xd[end])) 
		|> df -> [mean(df.bx), mean(df.by), mean(df.bz)]
		|> v -> swia_from_mso * v
		|> unitvec
	)
	# two perpandicular unit vectors
	evec = cross(bvec, bulkvelocity) |> unitvec
	vvec = cross(evec, bvec)

	size_inches = (9, 6)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 20, 9, 8))
    cn = cgrad(:gist_earth, 18, categorical = true)
	lcn = Makie.wong_colors()

    axes = [Axis(fig[1, i], xminorticksvisible=true, xminorticks = IntervalsBetween(9), yminorticksvisible=true, yminorticks = IntervalsBetween(9), yscale=log10, backgroundcolor=:gray80) for i in 1:3]

	# bvec direction
 	mvn_swia_eflux_bvec_profile = interpolate_mvn_swia_eflux_profile(mvn_swia, xd, bulkvelocity, bvec; coordinate="swia")
	scatter!(axes[1], mvn_swia_eflux_bvec_profile.velocity, mvn_swia_eflux_bvec_profile.def, markersize = 6, color = lcn[1], label = "Parallel")

	axes[1].ylabel = rich(rich("MVN SWIA H"), superscript("+"), rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")"))

	axes[1].xlabel = rich(rich("Parallel Velocity\n"), rich("(km⋅s"), superscript("-1"), rich(")"))

	# vvec direction
 	mvn_swia_eflux_bvec_profile = interpolate_mvn_swia_eflux_profile(mvn_swia, xd, bulkvelocity, vvec; coordinate="swia")
	scatter!(axes[2], mvn_swia_eflux_bvec_profile.velocity, mvn_swia_eflux_bvec_profile.def, markersize = 6, color = lcn[1], label = "Parallel")

	axes[2].xlabel = rich(rich("Perpendicular Velocity\n"), rich("(km⋅s"), superscript("-1"), rich(")"))

	# evec direction
 	mvn_swia_eflux_bvec_profile = interpolate_mvn_swia_eflux_profile(mvn_swia, xd, bulkvelocity, evec; coordinate="swia")
	scatter!(axes[3], mvn_swia_eflux_bvec_profile.velocity, mvn_swia_eflux_bvec_profile.def, markersize = 6, color = lcn[1], label = "Parallel")

	axes[3].xlabel = rich(rich("Perpendicular Velocity\n"), rich("(km⋅s"), superscript("-1"), rich(")"))

	linkyaxes!(axes...)
	hideydecorations!.(axes[2:3], grid=false)
	fig
end

"""
	Interpolating MVN SWIA eflux in the 3D velocity space with field-algined coordinate system

	mvn_swia: MVN SWIA data from read_mvn_swia_coarsesvy3d
	mvn_mag: MVN MAG data from read_mvn_mag
	xd: time interval for averaging

	bvec, vvec, and evec: swia coordinate system

	Note: Float64 data is required to ensure accuracy. 
"""
function interpolate_mvn_swia_eflux_3d(mvn_swia, mvn_mag, xd)
	# retriving MVN SWIA ion moments
	mvn_swia_ion_moments = calculate_mvn_swia_ion_moments(
		mvn_swia,
		xd;
		average = true,
		coordinate = "swia"
	)
	def = mvn_swia_ion_moments.def
	phi = mvn_swia_ion_moments.phi
	theta = mvn_swia_ion_moments.theta
	speed = mvn_swia_ion_moments.speed

	nphi, ntheta, nspeed = size(def)

	# Effetive theta range for interpolation
	theta_range = Matrix{Float64}(undef, 2, nspeed)  # rad
	for i in 1:nspeed
		theta_range[:, i] .= deg2rad.(extrema(theta[:, i]))
		half_width = deg2rad.((theta[1, i] - theta[2, i]) / 2)
		theta_range[1, i] -= half_width
		theta_range[2, i] += half_width
	end

	# Spheric Splines at all the energy channels
	ssp = map(1:nspeed) do ispeed
		dirs = [cartesian_from_spherical()(1.0, deg2rad(theta[itheta, ispeed]), deg2rad(phi[iphi]))
		for itheta in 1:ntheta for iphi in 1:nphi] |> x -> stack(x, dims = 1)
		vals = [Float64(def[iphi, itheta, ispeed]) for itheta in 1:ntheta for iphi in 1:nphi] .|> log10
		idx = findall(isfinite, vals)
		dir = view(dirs, idx, :)
		val = view(vals, idx)
		spherical_spline(dir, val)
	end

	function ssp_interp!(ssp, query_dir, logdefarr)
		for ispeed in 1:2:nspeed
			logdefarr[(ispeed + 1) ÷ 2] = length(ssp[ispeed].a) > 3 ? ssp[ispeed](query_dir) : NaN64
		end
		reverse!(logdefarr)
	end

	# bulkvelocity vector in SWIA coordinate 
	bulkvelocity = mvn_swia_ion_moments.bulkvelocity |> SVector{3, Float64}

	# tranfromation matrices
	mso_from_swia = mvn_swia_ion_moments.mso_from_swia
	swia_from_mso = transpose(mso_from_swia)

	# magnetic field unit vector in SWIA coordinate
	bvec = (mvn_mag 
		|> df -> subset(df, :epoch => ByRow(e -> xd[begin] ≤ e ≤ xd[end])) 
		|> df -> [mean(df.bx), mean(df.by), mean(df.bz)]
		|> v -> swia_from_mso * v
		|> SVector{3, Float64}
		|> unitvec
	)
	# two perpandicular unit vectors
	evec = cross(bvec, bulkvelocity) |> unitvec
	vvec = cross(evec, bvec)

	# interpolation
	step = 10.0
	vmax = 2000.0
	vmin = -2000.0
	vb = collect(range(vmin, vmax, step = step))
	ve = collect(range(vmin, vmax, step = step))
	vv = collect(range(vmin, vmax, step = step))

	def = Array{Float64, 3}(undef, length(vb), length(vv), length(ve))
	logvarr = reverse!(log.(speed[1:2:nspeed]))
	narr = length(logvarr)
	@batch per=thread for ive in 1:length(ve)
		logdefarr = Vector{Float64}(undef, narr)
		for ivv in 1:length(vv)
			for ivb in 1:length(vb)
				v = vb[ivb] * bvec + vv[ivv] * vvec + ve[ive] * evec
				logdef = NaN64
				vρ, vθ, vϕ = spherical_from_cartesian()(v[1], v[2], v[3])
				idx = searchsortedlast(speed, vρ, rev = true)
				if 1 < idx < 48 && theta_range[1, idx] ≤ vθ ≤ theta_range[2, idx]
					query_dir = collect(cartesian_from_spherical()(1.0, vθ, vϕ))
					ssp_interp!(ssp, query_dir, logdefarr)
					logdef = interpolate_weno4([log(vρ)], logvarr, logdefarr)[1]
				end
				def[ivb, ivv, ive] = exp10(logdef)
			end
		end
	end
	return (; def, vb, vv, ve, bvec, vvec, evec, mso_from_swia)
end

"""
	Plotting MVN SWIA eflux averaged over the specified time interval in the 3D velocity space

	mvn_swia_interp: MVN SWIA data interpolated to the static coordinate system
"""
function plot_mvn_swia_efluxmean_3d(mvn_swia_interp)
	def = mvn_swia_interp.def
	vb = mvn_swia_interp.vb
	vv = mvn_swia_interp.vv
	ve = mvn_swia_interp.ve

	# ============================================================
	# Unified colormap and colorrange settings
	# ============================================================
	cmap = :gist_earth
	crange = (1e3, 1e9)
	cscale = log10

	# ============================================================
	# Figure layout setup
	# ============================================================
	size_inches = (9, 20)
	size_pt = 72 .* size_inches
	fig = Figure(
		size = size_pt,
		fontsize = 15,
		figure_padding = (5, 10, 9, 8)
	)

	# ============================================================
	# Left panel: stacked heatmap slices (all energy channels)
	# ============================================================
	ls1 = LScene(fig[1, 1], show_axis = true)

	for (ive, vei) in enumerate(ve)
		heatmap!(ls1, vb, vv, def[:, :, ive],
			transformation = (:xy, vei),
			colormap = cmap,
			colorrange = crange,
			colorscale = cscale
		)
	end

	# ============================================================
	# Right panel: interactive volume slices with sliders
	# ============================================================
	ls2 = LScene(fig[1, 2], show_axis = true)

	plt = volumeslices!(ls2, vb, vv, ve, def,
		colormap = cmap,
		colorrange = crange,
		colorscale = cscale
	)

	# Unified colorbar for both panels
	cb = Colorbar(fig[2, 1], plt,
		label = rich("Differential Energy Flux (eV ⋅ cm", superscript("-2"), " s", superscript("-1"),
					" sr", superscript("-1"), " eV", superscript("-1"), ")"),
		vertical = false, flipaxis = false,
		ticklabelsize = 12,
		labelsize = 14,
		height = Relative(0.85)
	)
	

	# Slider grid to control slice positions along each axis
	sl_grid = SliderGrid(
		fig[2, 2],
		(label = "yz plane - v_B axis", range = 1:length(vb)),
		(label = "xz plane - v_V axis", range = 1:length(vv)),
		(label = "xy plane - v_E axis", range = 1:length(ve)),
	)

	lo = sl_grid.layout
	nc = ncols(lo)

	sl_yz, sl_xz, sl_xy = sl_grid.sliders

	# Connect sliders to the volume slice update methods
	on(sl_yz.value) do v
	plt[:update_yz][](v)
	end
	on(sl_xz.value) do v
	plt[:update_xz][](v)
	end
	on(sl_xy.value) do v
	plt[:update_xy][](v)
	end

	# Set initial slider positions to the middle of each axis
	set_close_to!(sl_yz, 0.5 * length(vb))
	set_close_to!(sl_xz, 0.5 * length(vv))
	set_close_to!(sl_xy, 0.5 * length(ve))

	# Add toggles to show/hide each slice plane
	hmaps = [plt[Symbol(:heatmap_, s)][] for s in (:yz, :xz, :xy)]
	toggles = [Toggle(lo[i, nc + 1], active = true) for i in 1:length(hmaps)]

	map(zip(hmaps, toggles)) do (h, t)
		on(t.active) do active
			h.visible = active
		end
	end

	return fig
end

"""
	Plotting MVN STATIC eflux overview

	mvn_static: MVN STATIC data from read_mvn_static_32e4d16a8m
	xd: x date ticks for plotting
	imass: mass index for the ion species
"""
function plot_mvn_static_eflux_overview(mvn_static, xd, imass)
	ind = findall(xd[begin] .<= mvn_static.epoch .<= xd[end])
    swpu = unique(mvn_static.swp[ind])
	if  length(swpu) > 1
		error("More than one swp found during the specified time interval. Please specify a time interval that corresponds to a single swp.")
	end
    swp = swpu[1]

    sp = [rich(rich("H"), superscript("+")), rich(rich("He"), superscript("++")), rich("m/q=4.573966"), rich("m/q=9.22163"), rich(rich("O"), superscript("+")), rich(rich("O"), subscript("2"), superscript("+")), rich(rich("CO"), subscript("2"), superscript("+")), rich("m/q=69.26344")]

    size_inches = (20, 30)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 20, figure_padding = (5,20, 9, 10))
    cn = cgrad(:seaborn_dark, 10, categorical = true)
    lcn = Makie.wong_colors()

    xtk = datetime2julian.(xd)
    xtk = Float32.(xtk .- mvn_static.julian[1])

    axes = [Axis(fig[j, i], yscale=log10, xticks = (xtk, Dates.format.(xd, "HH:MM")), backgroundcolor=:gray80, xticklabelrotation=-π/3) for j in 1:16, i in 1:4]

    for j in 1:16
        for i in 1:4
            ax = axes[j, i]
            heatmap!(ax, Float32.(mvn_static.julian.-mvn_static.julian[1]), mvn_static.energy[imass, j, i, :, swp], transpose(mvn_static.def[imass, j, i, :, :]), colormap=:gist_earth, colorscale=log10, colorrange = (1e5, 1e9))
            text!(ax, 0.0, 1.0, text = rich("ψ=$i ϕ=$j"), align = (:left, :top), offset = (4, -2), space = :relative, fontsize = 20, color = :white)
            xlims!(ax, xtk[begin], xtk[end]) 
            ylims!(ax, 1e1, 1e4)
        end
    end
    linkaxes!(axes...)
    hidexdecorations!.(axes[1:15,:], grid=false)
    hideydecorations!.(axes[:, 2:end], grid=false)
    colgap!(fig.layout, 30)
    Colorbar(fig[1:16, 5], colormap=:gist_earth, scale=log10, colorrange = (1e5, 1e9), label = rich("MAVEN STATIC ") * sp[imass] * rich(rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")))
    Label(fig[1:16, 0], rich("Energy (eV)"), fontsize = 25, halign = :center, rotation = π/2)
    Label(fig[17, 1:5], rich("Universal Time"), fontsize = 20, halign = :center)

    fig
end

"""
	Plotting MVN STATIC eflux averaged over the specified time interval

	mvn_static: MVN STATIC data from read_mvn_static_32e4d16a8m
	xd: time interval for averaging
	imass: mass index for the ion species
"""
function plot_mvn_static_efluxmean(mvn_static, xd, imass)
	ind = findall(xd[begin] .≤ mvn_static.epoch .≤ xd[end])
	efluxmean = dropdims(mean(mvn_static.def[:, :, :, :, ind], dims = 5), dims = 5)
    swpu = unique(mvn_static.swp[ind])
	if  length(swpu) > 1
		error("More than one swp found during the specified time interval. Please specify a time interval that corresponds to a single swp.")
	end
    swp = swpu[1]

    size_inches = (9, 18)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 20, 9, 8))
    cn = cgrad(:gist_earth, 18, categorical = true)

    axs = [Axis(fig[i, 1], xminorticksvisible=true, xminorticks = IntervalsBetween(9), yminorticksvisible=true, yminorticks = IntervalsBetween(9), xlabel = rich("Energy (eV)"), ylabel = rich(rich("H"), superscript("+"), rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")")), xscale = log10, yscale=log10, backgroundcolor=:gray80) for i in 1:4]
    for i in 1:4
        for j in 1:16
            ind = findall(efluxmean[imass, j, i, :] .> 0)
            scatterlines!(axs[i], mvn_static.energy[imass, j, i, ind, swp], efluxmean[imass, j, i, ind], color = cn[j+1], linewidth = 3)
        end
        text!(axs[i], 0.0, 1.0, text = "MAVEN STATIC ψ=$i", align = (:left, :top), offset = (4, -1), space = :relative, fontsize = 20, color = :black)
        ax1_inset = Axis(fig[i, 1], width=Relative(0.8), height=Relative(0.05), halign=0.5, valign=0.02)
        phi = reshape([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16], (16, 1))
        x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
        y = [1]
        heatmap!(ax1_inset, x, y, phi, colormap=cn, colorrange=(0,17))
        hidedecorations!(ax1_inset, grid=false)
        for ix in x
            text!(ax1_inset, ix, 1.0, text = "$ix", align = (:center, :center), offset = (0, 0), fontsize = 20, color = :black)
        end
        text!(axs[i], 0.042, 0.07, text = "φ=", align = (:left, :top), offset = (0, 0), space = :relative, fontsize = 20, color = :black)
    end
    hidexdecorations!.(axs[begin:end-1], grid=false)
    linkxaxes!(axs)
    xlims!.(axs, 1e0, 1e4)
    ylims!.(axs, 1e3, 1e9)

    fig
end


"""
	Interpolating MVN static eflux profile

	mvn_static: MVN static data from read_mvn_static_coarsesvy3d
	xd: time interval for averaging
	imass: mass index
	point: a point on the line (3D vector)
	dir: direction vector of the line (3D, non-zero; need not be normalized)
	coordinate: coordinate system of the point to be interpolated
"""
function interpolate_mvn_static_eflux_profile(mvn_static, xd, imass, point, dir; coordinate="static")
	mvn_static_ion_moments = calculate_mvn_static_ion_moments(
		mvn_static,
		xd,
		imass;
		average = true,
		coordinate = "static"
	)
	def = mvn_static_ion_moments.def
	phi = mvn_static_ion_moments.phi
	theta = mvn_static_ion_moments.theta
	speed = mvn_static_ion_moments.speed
	mso_from_static = mvn_static_ion_moments.mso_from_static
	static_from_mso = transpose(mso_from_static)

	#point = mvn_static_ion_moments.bulkvelocity
	#
	#dir = bvec
	# ind = argmax(def)
	# point = collect(cartesian_from_spherical()(speed[ind], deg2rad(theta[ind]), deg2rad(phi[ind])))
	# dir = collect(cartesian_from_spherical()(1.0, deg2rad(theta[ind]), deg2rad(phi[ind])))

	if coordinate == "static"
		p0 = point
		d = dir / norm(dir)
	elseif coordinate == "mso"
		p0 = static_from_mso * point
		d = static_from_mso * dir / norm(dir)
	else
		error("Invalid coordinate system")
	end

	nphi, ntheta, nspeed = size(def)

	def_profile = fill(NaN64, nspeed * 2)
	vel_profile = fill(NaN64, nspeed * 2)

	# Spheric Splines at all the energy channels
	@batch for ispeed in 1:nspeed
		#@show ispeed
		dirs = [cartesian_from_spherical()(1.0, deg2rad(theta[iphi, itheta, ispeed]), deg2rad(phi[iphi, itheta, ispeed]))
		for itheta in 1:ntheta for iphi in 1:nphi] |> x -> stack(x, dims = 1)
		vals = [Float64(def[iphi, itheta, ispeed]) for itheta in 1:ntheta for iphi in 1:nphi] .|> log10
		idx = findall(isfinite, vals)
		if length(idx) > 3
			ssp = spherical_spline(view(dirs, idx, :), view(vals, idx))
			# Effetive theta range for interpolation
			theta_range = collect(deg2rad.(extrema(theta[:, :, ispeed])))
			half_width = deg2rad.((theta[1, 2, ispeed] - theta[1, 1, ispeed]) / 2)
			theta_range[1] -= half_width
			theta_range[2] += half_width
			center = [0.0, 0.0, 0.0]
			speed1d = view(speed, 1, 1, :)
			radius = Float64(speed1d[ispeed])
			pinterp = intersect_line_sphere(p0, d, center, radius)
			if !isnothing(pinterp)
				for ipt in 1:2
					pt = pinterp[ipt]
					vel_profile[ispeed*2+ipt-2] = dot(pt, d)	
					vρ, vθ, vϕ = spherical_from_cartesian()(pt[1], pt[2], pt[3])
					idx = searchsortedlast(speed1d, vρ, rev = true)
					if 1 < idx < 32 && theta_range[1] ≤ vθ ≤ theta_range[2]
						query_dir = collect(cartesian_from_spherical()(1.0, vθ, vϕ))
						def_profile[ispeed*2+ipt-2] = exp10(ssp(query_dir))
					end
				end
			end
		end
	end
	idx = findall(isfinite, def_profile)
	def = def_profile[idx]
	velocity = vel_profile[idx]
	idx = sortperm(velocity)
	def = def[idx]
	velocity = velocity[idx]
	return (; def, velocity)
end

"""
	Plotting MVN static eflux profile

	mvn_static: MVN static data from read_mvn_static_coarsesvy3d
	mvn_mag: MVN magnetic field data from read_mvn_mag_coarsesvy3d
	xd: time interval for averaging
	imass: mass index
"""
function plot_mvn_static_eflux_profile(mvn_static, mvn_mag, xd, imass)
	# retriving MVN static ion moments
	mvn_static_ion_moments = calculate_mvn_static_ion_moments(
		mvn_static,
		xd,
		imass;
		average = true,
		coordinate = "static"
	)
	# tranfromation matrices
	mso_from_static = mvn_static_ion_moments.mso_from_static
	static_from_mso = transpose(mso_from_static)

	# magnetic field unit vector in static coordinate
	bvec = (mvn_mag 
		|> df -> subset(df, :epoch => ByRow(e -> xd[begin] ≤ e ≤ xd[end])) 
		|> df -> [mean(df.bx), mean(df.by), mean(df.bz)]
		|> v -> static_from_mso * v
		|> unitvec
	)

	# bulkvelocity vector in static coordinate defined by the largest def
	# bulkvelocity = mvn_static_ion_moments.bulkvelocity
	idx = argmax(mvn_static_ion_moments.def)
	bulkvelocity = mvn_static_ion_moments.velocity[:, idx]
	# two perpandicular unit vectors
	evec = cross(bvec, bulkvelocity) |> unitvec
	vvec = cross(evec, bvec)
	
	size_inches = (9, 6)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 15, figure_padding = (5, 20, 9, 8))
    cn = cgrad(:gist_earth, 18, categorical = true)
	lcn = Makie.wong_colors()

    axes = [Axis(fig[1, i], xminorticksvisible=true, xminorticks = IntervalsBetween(9), yminorticksvisible=true, yminorticks = IntervalsBetween(9), yscale=log10, backgroundcolor=:gray80) for i in 1:3]

	# bvec direction
 	mvn_static_eflux_bvec_profile = interpolate_mvn_static_eflux_profile(mvn_static, xd, imass, bulkvelocity, bvec; coordinate="static")
	scatter!(axes[1], mvn_static_eflux_bvec_profile.velocity, mvn_static_eflux_bvec_profile.def, markersize = 6, color = lcn[1], label = "Parallel")

	axes[1].ylabel = rich(rich("MVN Static H"), superscript("+"), rich(" Differential Energy Flux (cm"), superscript("-2"), rich("s"), superscript("-1"), rich("sr"), superscript("-1"), rich(")"))

	axes[1].xlabel = rich(rich("Parallel Velocity\n"), rich("(km⋅s"), superscript("-1"), rich(")"))

	# vvec direction
 	mvn_static_eflux_bvec_profile = interpolate_mvn_static_eflux_profile(mvn_static, xd, imass, bulkvelocity, vvec; coordinate="static")
	scatter!(axes[2], mvn_static_eflux_bvec_profile.velocity, mvn_static_eflux_bvec_profile.def, markersize = 6, color = lcn[1], label = "Parallel")

	axes[2].xlabel = rich(rich("Perpendicular Velocity\n"), rich("(km⋅s"), superscript("-1"), rich(")"))

	# evec direction
 	mvn_static_eflux_bvec_profile = interpolate_mvn_static_eflux_profile(mvn_static, xd, imass, bulkvelocity, evec; coordinate="static")
	scatter!(axes[3], mvn_static_eflux_bvec_profile.velocity, mvn_static_eflux_bvec_profile.def, markersize = 6, color = lcn[1], label = "Parallel")

	axes[3].xlabel = rich(rich("Perpendicular Velocity\n"), rich("(km⋅s"), superscript("-1"), rich(")"))

	linkyaxes!(axes...)
	hideydecorations!.(axes[2:3], grid=false)
	fig
end

"""
	Interpolating MVN STATIC eflux in the 3D velocity space with field-algined coordinate system

	mvn_static: MVN STATIC data from read_mvn_static
	mvn_mag: MVN MAG data from read_mvn_mag
	xd: time interval for averaging
	imass: ion mass index

	bvec, vvec, and evec: static coordinate system
	Note: Float64 data is required to ensure accuracy.
"""
function interpolate_mvn_static_eflux_3d(mvn_static, mvn_mag, xd, imass)
	# retriving MVN STATIC ion moments
	mvn_static_ion_moments = calculate_mvn_static_ion_moments(
		mvn_static,
		xd,
		imass;
		average = true,
		coordinate = "static"
	)
	def = mvn_static_ion_moments.def
	phi = mvn_static_ion_moments.phi
	theta = mvn_static_ion_moments.theta
	speed = mvn_static_ion_moments.speed
	# transformation matrices
	mso_from_static = mvn_static_ion_moments.mso_from_static
	static_from_mso = transpose(mso_from_static)
	# bulkvelocity vector in SWIA coordinate 
	bulkvelocity = mvn_static_ion_moments.bulkvelocity |> SVector{3, Float64}

	nphi, ntheta, nspeed = size(def)

	# Effetive theta range for interpolation
	theta_range = Matrix{Float64}(undef, 2, nspeed)  # rad
	for i in 1:nspeed
		theta_range[:, i] .= deg2rad.(extrema(view(theta, :, :, i)))
		half_width = deg2rad.((theta[1, 2, i] - theta[1, 1, i]) / 2)
		theta_range[1, i] -= half_width
		theta_range[2, i] += half_width
	end

	# Spheric Splines at all the energy channels
	ssp = map(1:nspeed) do ispeed
		dirs = [cartesian_from_spherical()(1.0, deg2rad(theta[iphi, itheta, ispeed]), deg2rad(phi[iphi, itheta, ispeed]))
		for iphi in 1:nphi for itheta in 1:ntheta] |> x -> stack(x, dims = 1)
		vals = [Float64(def[iphi, itheta, ispeed]) for iphi in 1:nphi for itheta in 1:ntheta] .|> log10
		idx = findall(isfinite, vals)
		dir = view(dirs, idx, :)
		val = view(vals, idx)
		spherical_spline(dir, val)
	end

	function ssp_interp!(ssp, query_dir, logdefarr)
		for ispeed in 1:nspeed
			logdefarr[ispeed] = length(ssp[ispeed].a) > 3 ? ssp[ispeed](query_dir) : NaN64
		end
		reverse!(logdefarr)
	end

	# magnetic field unit vector in STATIC coordinate
	bvec = (mvn_mag 
		|> df -> subset(df, :epoch => ByRow(e -> xd[begin] ≤ e ≤ xd[end])) 
		|> df -> [mean(df.bx), mean(df.by), mean(df.bz)]
		|> v -> static_from_mso * v
		|> SVector{3, Float64}
		|> unitvec
	)
	# two perpandicular unit vectors
	evec = cross(bvec, bulkvelocity) |> unitvec
	vvec = cross(evec, bvec)

	# interpolation
	step = 10.0
	vmax = 2000.0
	vmin = -2000.0
	vb = collect(range(vmin, vmax, step = step))
	ve = collect(range(vmin, vmax, step = step))
	vv = collect(range(vmin, vmax, step = step))

	def = Array{Float64, 3}(undef, length(vb), length(vv), length(ve))
	speed1d = view(speed, 1, 1, :)
	logvarr = reverse!(log.(speed1d))
	narr = length(logvarr)
	@batch per=thread for ive in 1:length(ve)
		logdefarr = Vector{Float64}(undef, narr)
		for ivv in 1:length(vv)
			for ivb in 1:length(vb)
				v = vb[ivb] * bvec + vv[ivv] * vvec + ve[ive] * evec
				logdef = NaN64
				vρ, vθ, vϕ = spherical_from_cartesian()(v[1], v[2], v[3])
				idx = searchsortedlast(speed1d, vρ, rev = true)
				if 1 < idx < 32 && theta_range[1, idx] ≤ vθ ≤ theta_range[2, idx]
					query_dir = collect(cartesian_from_spherical()(1.0, vθ, vϕ))
					ssp_interp!(ssp, query_dir, logdefarr)
					logdef = interpolate_weno4([log(vρ)], logvarr, logdefarr)[1]
				end
				def[ivb, ivv, ive] = exp10(logdef)
			end
		end
	end
	return (; def, vb, vv, ve, bvec, vvec, evec, mso_from_static)

end

"""
	Plotting MVN STATIC eflux averaged over the specified time interval in the 3D velocity space

	mvn_static_interp: MVN STATIC data interpolated to the static coordinate system
"""
function plot_mvn_static_efluxmean_3d(mvn_static_interp)
	def = mvn_static_interp.def
	vb = mvn_static_interp.vb
	vv = mvn_static_interp.vv
	ve = mvn_static_interp.ve

	# ============================================================
	# Unified colormap and colorrange settings
	# ============================================================
	cmap = :gist_earth
	crange = (1e3, 1e9)
	cscale = log10

	# ============================================================
	# Figure layout setup
	# ============================================================
	size_inches = (9, 20)
	size_pt = 72 .* size_inches
	fig = Figure(
		size = size_pt,
		fontsize = 15,
		figure_padding = (5, 10, 9, 8)
	)

	# ============================================================
	# Left panel: stacked heatmap slices (all energy channels)
	# ============================================================
	ls1 = LScene(fig[1, 1], show_axis = true)

	for (ive, vei) in enumerate(ve)
		heatmap!(ls1, vb, vv, def[:, :, ive],
			transformation = (:xy, vei),
			colormap = cmap,
			colorrange = crange,
			colorscale = cscale
		)
	end

	# ============================================================
	# Right panel: interactive volume slices with sliders
	# ============================================================
	ls2 = LScene(fig[1, 2], show_axis = true)

	plt = volumeslices!(ls2, vb, vv, ve, def,
		colormap = cmap,
		colorrange = crange,
		colorscale = cscale
	)

	# Unified colorbar for both panels
	cb = Colorbar(fig[2, 1], plt,
		label = rich("Differential Energy Flux (eV ⋅ cm", superscript("-2"), " s", superscript("-1"),
					" sr", superscript("-1"), " eV", superscript("-1"), ")"),
		vertical = false, flipaxis = false,
		ticklabelsize = 12,
		labelsize = 14,
		height = Relative(0.85)
	)
	

	# Slider grid to control slice positions along each axis
	sl_grid = SliderGrid(
		fig[2, 2],
		(label = "yz plane - v_B axis", range = 1:length(vb)),
		(label = "xz plane - v_V axis", range = 1:length(vv)),
		(label = "xy plane - v_E axis", range = 1:length(ve)),
	)

	lo = sl_grid.layout
	nc = ncols(lo)

	sl_yz, sl_xz, sl_xy = sl_grid.sliders

	# Connect sliders to the volume slice update methods
	on(sl_yz.value) do v
	plt[:update_yz][](v)
	end
	on(sl_xz.value) do v
	plt[:update_xz][](v)
	end
	on(sl_xy.value) do v
	plt[:update_xy][](v)
	end

	# Set initial slider positions to the middle of each axis
	set_close_to!(sl_yz, 0.5 * length(vb))
	set_close_to!(sl_xz, 0.5 * length(vv))
	set_close_to!(sl_xy, 0.5 * length(ve))

	# Add toggles to show/hide each slice plane
	hmaps = [plt[Symbol(:heatmap_, s)][] for s in (:yz, :xz, :xy)]
	toggles = [Toggle(lo[i, nc + 1], active = true) for i in 1:length(hmaps)]

	map(zip(hmaps, toggles)) do (h, t)
		on(t.active) do active
			h.visible = active
		end
	end

	return fig
end

"""
	Validating the moments calculated from SWIA and STATIC data. The moments should be consistent with the datasets from KP.

	Input:
		mvn_swia_ion_moments: moments calculated from SWIA data
		mvn_static_ion_moments: moments calculated from STATIC data
		mvn_kp: KP data from read_mvn_kp
		xd: x-ticks for plotting
"""
function validate_moments(mvn_swia_ion_moments, mvn_static_ion_moments, mvn_kp, xd)
	size_inches = (9, 16)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 20, figure_padding = (5,20, 9, 10))
    cn = cgrad(:seaborn_dark, 10, categorical = true)
    lcn = Makie.wong_colors()

	xst = datetime2julian(xd[1])
    xtk = datetime2julian.(xd)
    xtk = Float32.(xtk .- xst)

	axes = [Axis(fig[i, 1], backgroundcolor=:gray80) for i in 1:5]
	axes[1].ylabel = rich(rich("Density (cm"), superscript("-3"), rich(")"))
	#axes[2].ylabel = rich(rich("Bulk Velocity MSO X (km⋅s"), superscript("-1"), rich(")"))
	axes[3].ylabel = rich(rich("Bulk Velocity MSO (km⋅s"), superscript("-1"), rich(")"))
	#axes[4].ylabel = rich(rich("Bulk Velocity MSO Z (km⋅s"), superscript("-1"), rich(")"))
	axes[5].ylabel = rich("Temperature (eV)")

	lines!(axes[1], mvn_kp.julian .- xst, mvn_kp[:, 2], linewidth = 2, color = lcn[1], label = "KP LPW e⁻")
	lines!(axes[1], mvn_kp.julian .- xst, mvn_kp[:, 23], linewidth = 2, color = lcn[2], label = "KP SWEA e⁻")
	lines!(axes[1], mvn_kp.julian .- xst, mvn_kp[:, 41], linewidth = 2, color = lcn[3], label = "KP SWIA H⁺")
	lines!(axes[1], mvn_kp.julian .- xst, mvn_kp[:, 54], linewidth = 2, color = lcn[4], label = "KP STATIC H⁺")
	lines!(axes[1], datetime2julian.(mvn_swia_ion_moments.epoch) .- xst, mvn_swia_ion_moments.density, linewidth = 2, color = lcn[5], label = "SWIA H⁺")
	lines!(axes[1], datetime2julian.(mvn_static_ion_moments.epoch) .- xst, mvn_static_ion_moments.density, linewidth = 2, color = lcn[6], label = "STATIC H⁺")
	axislegend(axes[1])

	lines!(axes[2], mvn_kp.julian .- xst, mvn_kp[:, 43], linewidth = 2, color = lcn[3], label = "KP SWIA H⁺")
	lines!(axes[2], datetime2julian.(mvn_swia_ion_moments.epoch) .- xst, mvn_swia_ion_moments.bulkvelocity[1, :], linewidth = 2, color = lcn[5], label = "SWIA H⁺")
	lines!(axes[2], datetime2julian.(mvn_static_ion_moments.epoch) .- xst, mvn_static_ion_moments.bulkvelocity[1, :], linewidth = 2, color = lcn[6], label = "STATIC H⁺")

	lines!(axes[3], mvn_kp.julian .- xst, mvn_kp[:, 45], linewidth = 2, color = lcn[3], label = "KP SWIA H⁺")
	lines!(axes[3], datetime2julian.(mvn_swia_ion_moments.epoch) .- xst, mvn_swia_ion_moments.bulkvelocity[2, :], linewidth = 2, color = lcn[5], label = "SWIA H⁺")
	lines!(axes[3], datetime2julian.(mvn_static_ion_moments.epoch) .- xst, mvn_static_ion_moments.bulkvelocity[2, :], linewidth = 2, color = lcn[6], label = "STATIC H⁺")

	lines!(axes[4], mvn_kp.julian .- xst, mvn_kp[:, 47], linewidth = 2, color = lcn[3], label = "KP SWIA H⁺")
	lines!(axes[4], datetime2julian.(mvn_swia_ion_moments.epoch) .- xst, mvn_swia_ion_moments.bulkvelocity[3, :], linewidth = 2, color = lcn[5], label = "SWIA H⁺")
	lines!(axes[4], datetime2julian.(mvn_static_ion_moments.epoch) .- xst, mvn_static_ion_moments.bulkvelocity[3, :], linewidth = 2, color = lcn[6], label = "STATIC H⁺")

	lines!(axes[5], mvn_kp.julian .- xst, mvn_kp[:, 49], linewidth = 2, color = lcn[3], label = "KP SWIA H⁺")
	lines!(axes[5], mvn_kp.julian .- xst, mvn_kp[:, 60], linewidth = 2, color = lcn[4], label = "KP STATIC H⁺")
	lines!(axes[5], datetime2julian.(mvn_swia_ion_moments.epoch) .- xst, mvn_swia_ion_moments.temperature, linewidth = 2, color = lcn[5], label = "SWIA H⁺")
	lines!(axes[5], datetime2julian.(mvn_static_ion_moments.epoch) .- xst, mvn_static_ion_moments.temperature, linewidth = 2, color = lcn[6], label = "STATIC H⁺")

	axes[5].xticks = (xtk, Dates.format.(xd, "HH:MM"))
	xlims!(axes[5], xtk[begin], xtk[end])
	ylims!(axes[1], 1, 10)
	ylims!(axes[2], -400, 0)
	ylims!(axes[3], -100, 100)
	ylims!(axes[4], -100, 200)
	ylims!(axes[5], 0, 300)
	linkxaxes!(axes...)
	hidexdecorations!.(axes[1:4])
	fig
end

end

# # Specifying date
# date_str = "20220519"
# date = DateTime(date_str, "yyyymmdd")
# date_julian = datetime2julian(date)

# # Reading data
# mvn_mag = read_mvn_mag(date_str) # DataFrame
# mvn_kp_param, mvn_kp = read_mvn_kp(date_str) #(DataFrame, DataFrame)
# filter(row -> occursin(r"rotation"i, row.parameter) && occursin(r"spacecraft"i, row.parameter), mvn_kp_param)
# temp = filter(row -> occursin(r"shape"i, row.parameter), mvn_kp_param)
# filter(row -> occursin(r"density"i, row.parameter), mvn_kp_param)
# filter(row -> occursin(r"velocity"i, row.parameter), mvn_kp_param)
# filter(row -> occursin(r"temperature"i, row.parameter), mvn_kp_param)
# mvn_swia = read_mvn_swia_coarsesvy3d(date_str, mvn_kp)
# mvn_static = read_mvn_static_32e4d16a8m(date_str) 

# # Plotting 2D data
# CairoMakie.activate!()
# xd = DateTime(2022, 05, 19, 02, 50, 00) + Minute.(0:2:10)
# fig = plot_mvn_swia_eflux_overview(mvn_swia, xd)
# save("mvn_swia_eflux_overview_plot_" * Dates.format(xd[1], "yyyymmddHHMMSS") * ".pdf", fig, pt_per_unit = 1)

# xd = DateTime(2022, 05, 19, 02, 52, 00) + Minute.([0, 3])
# fig = plot_mvn_swia_efluxmean(mvn_swia, xd)
# save("mvn_swia_efluxmean_plot_" * Dates.format(xd[1], "yyyymmddHHMMSS") * ".pdf", fig, pt_per_unit = 1)

# imass = 1 # H+
# xd = DateTime(2022, 05, 19, 02, 50, 00) + Minute.(0:2:10)
# fig = plot_mvn_static_eflux_overview(mvn_static, xd, imass)
# save("mvn_static_eflux_overview_plot_" * Dates.format(xd[1], "yyyymmddHHMMSS") * ".pdf", fig, pt_per_unit = 1)

# imass = 1 # H+
# xd = DateTime(2022, 05, 19, 02, 52, 00) + Minute.([0, 3])
# fig = plot_mvn_static_efluxmean(mvn_static, xd, imass)
# save("mvn_static_efluxmean_plot_" * Dates.format(xd[1], "yyyymmddHHMMSS") * ".pdf", fig, pt_per_unit = 1)

# # Calculating moments
# xd = DateTime(2022, 05, 19, 02, 00, 00) + Minute.(0:10:60)
# mvn_swia_ion_moments = calculate_mvn_swia_ion_moments(mvn_swia, xd)
# imass = 1 # H+, caucious to use, since the mass value is not always the same for H+ in different swp. Check mvn_static.mass to find the correct mass index for H+.
# mvn_static_ion_moments = calculate_mvn_static_ion_moments(mvn_static, xd, imass)
# fig = validate_moments(mvn_swia_ion_moments, mvn_static_ion_moments, mvn_kp, xd)
# save("validate_moments_" * Dates.format(xd[1], "yyyymmddHHMMSS") * ".pdf", fig, pt_per_unit = 1)

# # Plotting 1D profiles
# xd = DateTime(2022, 05, 19, 02, 20, 00) + Minute.(0:2:10)
# fig = plot_mvn_swia_eflux_profile(mvn_swia, mvn_mag, xd)
# save("mvn_swia_eflux_profile_plot_" * Dates.format(xd[1], "yyyymmddHHMMSS") * ".pdf", fig, pt_per_unit = 1)
# imass = 1 # H+
# xd = DateTime(2022, 05, 19, 02, 50, 00) + Minute.(0:2:10)
# fig = plot_mvn_static_eflux_profile(mvn_static, mvn_mag, xd, imass)
# save("mvn_static_eflux_profile_plot_" * Dates.format(xd[1], "yyyymmddHHMMSS") * ".pdf", fig, pt_per_unit = 1)

# # Plotting 3D data
# GLMakie.activate!()
# xd = DateTime(2022, 05, 19, 02, 20, 00) + Minute.(0:2:10)
# mvn_swia_interp = interpolate_mvn_swia_eflux_3d(mvn_swia, mvn_mag, xd)
# mvn_swia_efluxmean_3d_plot = plot_mvn_swia_efluxmean_3d(mvn_swia_interp)
# imass = 1 # H+
# xd = DateTime(2022, 05, 19, 02, 50, 00) + Minute.(0:2:10)
# mvn_static_interp = interpolate_mvn_static_eflux_3d(mvn_static, mvn_mag, xd, imass)
# mvn_static_efluxmean_3d_plot = plot_mvn_static_efluxmean_3d(mvn_static_interp)