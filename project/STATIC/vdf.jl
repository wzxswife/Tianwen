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

include("../../src/scripts/MAVEN_load.jl")
include("../../src/scripts/MAVEN_STATIC.jl")
include("../../src/scripts/SphericalSplines.jl")
include("../../src/scripts/PhysicsCoordTransforms.jl")
include("../../src/scripts/SphericalSplinesOptim.jl")
include("../../src/scripts/MVN3DIonDisttribution.jl")
using .MAVEN_STATIC
using .MAVEN_load
using .SphericalSplines
using .PhysicsCoordTransforms
using .SphericalSplinesOptim
using .MVN3DIonDisttribution

date = Date(2022, 8, 24)
datestr = Dates.format.(date, "yyyymmdd")

ts = DateTime(date) + Dates.Hour(7) + Dates.Minute(26)
te = ts + Dates.Minute(3)

output_dir = "E:\\Output"
MAVEN_output_path = joinpath(output_dir, "MAVEN")
STATIC_output_path = joinpath(MAVEN_output_path, "STATIC")
VDF_path = joinpath(STATIC_output_path, "VDF", datestr[1:4], datestr[5:6])
mkpath(VDF_path)

mvn_mag = MVN3DIonDisttribution.read_mvn_mag(
    datestr; 
    data_root="E:/Data/MAVEN", 
    resolution=""
    )
mvn_static = MVN3DIonDisttribution.read_mvn_static_32e4d16a8m(
    datestr; 
    data_root="E:/Data/MAVEN"
    )

println("Data loaded.")

imass = 1
xd = ts + Minute.(1:3)

mvn_static_interp = MVN3DIonDisttribution.interpolate_mvn_static_eflux_3d(mvn_static, mvn_mag, xd, imass)
println("Data interpolated.")
fig = MVN3DIonDisttribution.plot_mvn_static_efluxmean_3d(mvn_static_interp)
println("Data plotted.")

while isopen(fig.scene)
    sleep(0.1)
end