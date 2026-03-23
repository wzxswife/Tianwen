using Dates
using DataFrames
using Statistics
using LinearAlgebra

include("./scripts/MAVEN_load.jl")
include("./scripts/wave_caculate.jl")
using .MAVEN_load
using .WaveCaculate

data_path = joinpath(@__DIR__, "..", "data", "MAVEN")
file_name = "mvn_mag_l2_2022236ss_20220824_v01_r01.sts"
mag2c32hz = load_mag_l2(joinpath(data_path, file_name))
mag2c32hz[:position] .= mag2c32hz[:position] ./ 3390.0
println("Data loaded.")

