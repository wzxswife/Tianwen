using Dates, Glob
using GLMakie, CairoMakie
include("../../src/scripts/TW_MINPA.jl")
include("../../src/scripts/TW_load.jl")
include("../../src/scripts/MINPA_plot.jl")
using .TW_MINPA
using .TW_load
using .MINPA_plot

date = Date(2022, 8, 24)
datestr = Dates.format.(date, "yyyymmdd")
year_num = Dates.format(date, "yyyy")
month_num = Dates.format(date, "mm")

data_dir = "E:\\Data"
output_dir = "E:\\Output"
TW_path = joinpath(data_dir, "Tianwen-1")
MINPA_path = joinpath(TW_path, "MINPA", "2B", year_num, month_num)
MOMAG_path = joinpath(TW_path, "MOMAG", "2C.V03", year_num, month_num)
TW_output_path = joinpath(output_dir, "Tianwen-1")
MINPA_output_path = joinpath(TW_output_path, "MINPA")
overview_path = joinpath(MINPA_output_path, "Overview", year_num, month_num)
mkpath(overview_path)
MINPA_file_list = glob("*$datestr*.2B", MINPA_path)
MOMAG_file = glob("*01Hz*$datestr*.dat", MOMAG_path)

println("MOMAG file: ", MOMAG_file)
println("MINPA files found: ", length(MINPA_file_list))
for file in MINPA_file_list
    println("  ", file)
end
println("Output path: ", overview_path)

file_list = Dict{String,Any}(
    "20220824" => begin
        date = DateTime(2022, 8, 24)
        Dict{Symbol,Any}(
            :date    => date,
            :datestr => Dates.format(date, "yyyymmdd"),
            :file1   => begin
                ts = [date + Dates.Hour(0), 
                    date + Dates.Hour(2), 
                    date + Dates.Hour(4),
                    date + Dates.Hour(6)]
                Dict{Symbol,Any}(
                    :mode => 1, 
                    :file_name => MINPA_file_list[1],
                    :ts_list => ts, 
                    :time_list => [ts[1].+ Dates.Minute(20).*(0:6), 
                    ts[2].+ Dates.Minute(20).*(0:6), 
                    ts[3].+ Dates.Minute(20).*(0:6),
                    ts[4].+ Dates.Minute(10).*(0:6)])
            end,
            :file2   => begin
                ts = [date + Dates.Hour(8), 
                    date + Dates.Hour(10), 
                    date + Dates.Hour(12)]
                Dict{Symbol,Any}(
                    :mode => 1, 
                    :file_name => MINPA_file_list[2],
                    :ts_list => ts, 
                    :time_list => [t .+ Dates.Minute(20) .* (0:6) for t in ts])
            end,
            :file4   => begin
                ts = [date + Dates.Hour(10) + Dates.Minute(53)]
                Dict{Symbol,Any}(
                    :mode => 12, 
                    :file_name => MINPA_file_list[4],
                    :ts_list => ts, 
                    :time_list => [t .+ Dates.Minute(5) .* (0:6) for t in ts])
            end,
            :file6   => begin
                ts = [date + Dates.Hour(7)+Dates.Minute(30)]
                Dict{Symbol,Any}(
                    :mode => 4, 
                    :file_name => MINPA_file_list[6],
                    :ts_list => ts, 
                    :time_list => [t .+ Dates.Minute(5) .* (0:6) for t in ts])
            end,
        )
    end
)

mag_data = TW_load.load_mag_2c_bydlm(MOMAG_file[1])
minpa_Info = file_list[datestr][:file6]
mod = minpa_Info[:mode]
minpa_data = TW_load.load_minpa(minpa_Info[:file_name], mod)

h_list = minpa_Info[:ts_list]
img_list = [1, 2, 5, 6, 7]

sp_title = ["H+", "He++", "m/q=4.573966", "m/q=9.22163", 
    "O+", "O2+", "CO2+", "m/q=69.26344"]
sp = [rich(rich("H"), superscript("+")), rich(rich("He"), superscript("++")), 
    rich("m/q=4.573966"), rich("m/q=9.22163"), rich(rich("O"), superscript("+")), 
    rich(rich("O"), subscript("2"), superscript("+")), rich(rich("CO"), 
    subscript("2"), superscript("+")), rich("m/q=69.26344")]

using CairoMakie
CairoMakie.activate!()
set_theme!(;
    Axis=(; bordercolor=:black, bordersize=15, spinewidth = 3, 
        xlabelsize=18, ylabelsize=18, xticklabelsize=16, yticklabelsize=16)
)

for img in img_list
    for i in 1:length(h_list)
        xd = minpa_Info[:time_list][i]
        name = "MINPA_mod" * string(mod) * "_" * sp_title[img] * "eflux_Overview" * Dates.format(xd[1], "yyyymmddTHH") * ".pdf"
        fig = TW_MINPA_eflux_Overview_plot(minpa_data, xd, img)
        save(joinpath(overview_path, name), fig, pt_per_unit = 1)
        println("Saved figure: ", joinpath(overview_path, name))
    end
end