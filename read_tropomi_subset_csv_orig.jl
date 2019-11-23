# $hplin$2019-11-14$
using NCDatasets
using DataFrames
using CSV

output_lons = Float32[];
output_lats = Float32[];
output_qa   = Float32[];
output_no2co= Float32[];
output_time = String[];

files = readdir("D:\\")
for file in files
    if !occursin(r"nc", file)
        continue
    end

    filename = "/d/" * file
    println("Reading ", filename)
    ds = Dataset(filename, "r")
    prod = ds.group["PRODUCT"];

    # Copy the lat-lon data into local arrays
    no2cols = convert(Array{Union{Float32, Missing}, 3}, prod["nitrogendioxide_tropospheric_column"]);
    lons = convert(Array{Float32, 3}, prod["longitude"]);
    lats = convert(Array{Float32, 3}, prod["latitude"]);
    qavals = convert(Array{Float32, 3}, prod["qa_value"]);
    utctimes = convert(Array{String, 1}, variable(prod, "time_utc")[1:end, 1]); # 1x1, JM

    IM, JM, LM = size(no2cols)

    println("Loaded file with dim'l IM = ", IM, " JM = ", JM, " LM = ", LM)

    for i in 1:IM
        for j in 1:JM
            if isequal(missing, no2cols[i, j, 1]) || qavals[i, j, 1] < 0.5
                continue
            end

            if lats[i, j, 1] >= 42 && lats[i, j, 1] <= 43 && lons[i, j, 1] >= -72 && lons[i, j, 1] <= -70
                push!(output_lons,  lons[i, j, 1])
                push!(output_lats,  lats[i, j, 1])
                push!(output_no2co, no2cols[i, j, 1])
                push!(output_qa,    qavals[i, j, 1])
                push!(output_time,  utctimes[j])
                # println("Data at i = ", i, " j = ", j, " lat = ", lats[i, j, 1], " lon = ", lons[i, j, 1], " val = ", no2cols[i, j, 1], " qa_value = ", qavals[i, j, 1])
                print("*")
            else
                continue
            end
        end
    end
    println("OK! Total of ", length(output_lons), " data points saved from this file.")
end

output = Dict{Symbol,Union{Array{Float32, 1}, Array{String, 1}}}(:lats => output_lats, :lons => output_lons, :no2 => output_no2co, :qa => output_qa, :time => output_time)
output_df = DataFrame(;output...)
CSV.write("total_output_2019_11_14.csv",  DataFrame(output_df), writeheader=true)
