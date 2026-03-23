using DataFrames, CSV

dir = joinpath(@__DIR__, "..", "..")
data_path = joinpath(dir, "Data")

println(dir)
println(data_path)

data_name = joinpath(data_path, "1HzWaves.csv")
df = CSV.read(data_name, DataFrame)

println("Data loaded.")