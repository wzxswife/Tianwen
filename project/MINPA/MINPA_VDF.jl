using Dates, Glob

TW_path = "E:/Data/Tianwen-1/"
MAG_path = joinpath(TW_path, "MOMAG/2C.V03/2022/08/")
MINPA_path = joinpath(TW_path, "MINPA/2B/2022/08/")
MINPA_file_list = glob("*20220824*.2B", MINPA_path)

println("MINPA files found: ", length(MINPA_file_list))
for file in MINPA_file_list
    println("  ", file)
end