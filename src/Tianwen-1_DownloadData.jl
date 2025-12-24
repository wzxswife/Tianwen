using Downloads
using Dates

# date = DateTime(2021, 11, 16, 00, 00, 00) .+ Dates.Day(1) .* range(0, 45)
# date = DateTime(2022, 01, 01, 00, 00, 00) .+ Dates.Day(1) .* range(0, 180) #20220630
# date = DateTime(2021, 11, 16, 00, 00, 00) .+ Dates.Day(1) .* range(0, 226) #20220630
# date = DateTime(2022, 7, 1, 00, 00, 00) .+ Dates.Day(1) .* range(0, 125) #20221103
date = DateTime(2023, 4, 6, 00, 00, 00) .+ Dates.Day(1) .* range(0, 54)
datestr = Dates.format.(date, "yyyymmdd")

for si in datestr
    local link = "https://space.ustc.edu.cn/dreams/tw1_momag/fetch.php?datafile=TW1_MOMAG_MSO_32Hz_"*si*"_2C_v03.dat"
    Downloads.download(link, "E:/Acode/JuliaCode/Tianwen/Data/32Hz/TW1_MOMAG_MSO_32Hz_"*si*"_2C_v03.dat")
    println(link)
end