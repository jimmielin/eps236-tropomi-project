# $hplin$2019-11-22$
using NCDatasets
using DataFrames
using CSV

output_lons   = Float32[];
output_lats   = Float32[];
output_no2co  = Float32[];
output_no2err = Float32[];
output_qa     = Float32[];
output_time   = String[];
output_timestamp = Int32[];
output_cloud  = Float32[];

output_cornerlon1 = Float32[];
output_cornerlon2 = Float32[];
output_cornerlon3 = Float32[];
output_cornerlon4 = Float32[];
output_cornerlat1 = Float32[];
output_cornerlat2 = Float32[];
output_cornerlat3 = Float32[];
output_cornerlat4 = Float32[];

files = readdir("D:\\")
for file in files
    if !occursin(r"nc", file) # only read netCDF files
        continue
    end

    filename = "/d/" * file
    println("Reading ", filename)
    ds = Dataset(filename, "r")
    prod = ds.group["PRODUCT"];

    # See if we can loop through lon, lat without convert, which is a huge performance hit
    IM, JM, LM = size(prod["nitrogendioxide_tropospheric_column"])

    println("Loaded file with dim'l IM = ", IM, " JM = ", JM, " LM = ", LM)

    @time for i in 1:IM
        for j in 1:JM
            if isequal(missing, prod["nitrogendioxide_tropospheric_column"][i, j, 1]) || prod["qa_value"][i, j, 1] < 0.5
                continue
            end

            if prod["latitude"][i, j, 1] >= 42 && prod["latitude"][i, j, 1] <= 43 && prod["longitude"][i, j, 1] >= -72 && prod["longitude"][i, j, 1] <= -70
                push!(output_lons,  prod["longitude"][i, j, 1])
                push!(output_lats,  prod["latitude"][i, j, 1])
                push!(output_no2co, prod["nitrogendioxide_tropospheric_column"][i, j, 1])
                push!(output_no2err, prod["nitrogendioxide_tropospheric_column_precision"][i, j, 1])
                push!(output_qa,    prod["qa_value"][i, j, 1])
                push!(output_time,  variable(prod, "time_utc")[j, 1])
                push!(output_timestamp, trunc(Int32, Dates.datetime2unix(DateTime(split(output_time[end], ".")[1], "yyyy-mm-ddTHH:MM:SS"))))
                push!(output_cloud, prod.group["SUPPORT_DATA"].group["DETAILED_RESULTS"]["cloud_fraction_crb_nitrogendioxide_window"][i, j, 1])

                push!(output_cornerlon1, prod.group["SUPPORT_DATA"].group["GEOLOCATIONS"]["longitude_bounds"][1, i, j, 1])
                push!(output_cornerlon2, prod.group["SUPPORT_DATA"].group["GEOLOCATIONS"]["longitude_bounds"][2, i, j, 1])
                push!(output_cornerlon3, prod.group["SUPPORT_DATA"].group["GEOLOCATIONS"]["longitude_bounds"][3, i, j, 1])
                push!(output_cornerlon4, prod.group["SUPPORT_DATA"].group["GEOLOCATIONS"]["longitude_bounds"][4, i, j, 1])
                push!(output_cornerlat1, prod.group["SUPPORT_DATA"].group["GEOLOCATIONS"]["latitude_bounds"][1, i, j, 1])
                push!(output_cornerlat2, prod.group["SUPPORT_DATA"].group["GEOLOCATIONS"]["latitude_bounds"][2, i, j, 1])
                push!(output_cornerlat3, prod.group["SUPPORT_DATA"].group["GEOLOCATIONS"]["latitude_bounds"][3, i, j, 1])
                push!(output_cornerlat4, prod.group["SUPPORT_DATA"].group["GEOLOCATIONS"]["latitude_bounds"][4, i, j, 1])

                println("Data at i = ", i, " j = ", j, " lat = ", output_lats[end], " lon = ", output_lons[end], " val = ", output_no2co[end], " qa_value = ", output_qa[end])
            else
                continue
            end
        end
    end

    println("OK! Total of ", length(output_lons), " data points saved.")
end

output = Dict{Symbol,Union{Array{Float32, 1}, Array{String, 1}}}(
    :lats => output_lats,
    :lons => output_lons,
    :no2 => output_no2co,
    :no2_error => output_no2err,
    :qa => output_qa,
    :time => output_time,
    :timestamp => output_timestamp,
    :cloud => output_cloud,

    :cornerlon1 => output_cornerlon1,
    :cornerlon2 => output_cornerlon2,
    :cornerlon3 => output_cornerlon3,
    :cornerlon4 => output_cornerlon4,
    :cornerlat1 => output_cornerlat1,
    :cornerlat2 => output_cornerlat2,
    :cornerlat3 => output_cornerlat3,
    :cornerlat4 => output_cornerlat4,
)

output_df = DataFrame(;output...)
CSV.write("total_output_2019_11_22.csv",  DataFrame(output_df), writeheader=true)
