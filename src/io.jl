# ==============================================================================
# RASTER I/O
# ==============================================================================

"""
    read_raster(path::String, T::Type=Float64) -> (array, wkt, transform)

Read a single-band GeoTIFF. Returns a numeric 2D array (rows × cols),
the WKT projection string, and the geotransform vector.

The array is transposed from ArchGDAL's column-major layout to match
Julia's row-major geographic convention (row = y, col = x).
"""
function read_raster(path::String, T::Type=Float64)
    !isfile(path) && error("the file \"$(path)\" does not exist")

    ArchGDAL.read(path) do raw
        transform = ArchGDAL.getgeotransform(raw)
        wkt = ArchGDAL.getproj(raw)

        band = ArchGDAL.getband(raw, 1)
        array_t = ArchGDAL.read(band)

        # Transpose: ArchGDAL returns (x, y), we want (row, col)
        array = convert(Array{T, 2}, permutedims(array_t, [2, 1]))

        array, wkt, transform
    end
end

"""
    write_raster(path::String, array::Matrix, wkt::String, transform)

Write a 2D array as a single-band GeoTIFF with LZW compression.
Sets nodata to -9999. NaN values in the array are written as -9999.
"""
function write_raster(path::String,
                      array::Array{T, 2},
                      wkt::String,
                      transform) where T <: Number

    # Replace NaN with nodata sentinel for output
    out = copy(array)
    out[isnan.(out)] .= T(-9999)

    # Transpose back to ArchGDAL's (x, y) layout
    array_t = permutedims(out, [2, 1])
    width, height = size(array_t)

    options = ["COMPRESS=LZW"]

    # Ensure path ends with .tif
    fn = endswith(path, ".tif") ? path : string(path, ".tif")

    ArchGDAL.create(fn,
                    driver = ArchGDAL.getdriver("GTiff"),
                    width = width,
                    height = height,
                    nbands = 1,
                    dtype = T,
                    options = options) do dataset
        band = ArchGDAL.getband(dataset, 1)
        ArchGDAL.write!(band, array_t)
        ArchGDAL.setnodatavalue!(band, -9999.0)
        ArchGDAL.setgeotransform!(dataset, collect(Float64, transform))
        ArchGDAL.setproj!(dataset, wkt)
    end
end

"""
    read_reclass_table(path::String) -> Dict{Float64,Float64}

Read a 2-column delimited file mapping original resistance values to new values.
Lines where the second column is `missing` or `nodata` map to -9999.0 (nodata).
"""
function read_reclass_table(path::String; nodata_value::Float64=-9999.0)
    !isfile(path) && error("the file \"$(path)\" does not exist")
    raw = readdlm(path)
    
    table = Dict{Float64, Float64}()

    for i in 1:size(raw, 1)
        key = Float64(raw[i, 1])
        val = raw[i, 2]
        if val isa AbstractString && lowercase(val) in ("missing", "nodata", "")
            table[key] = nodata_value
        else
            table[key] = Float64(val)
        end
    end
    table
end

"""
    write_cfg(cfg::Dict{String,String}, path::String)

Write the configuration dictionary to an .ini file for reproducibility.
All keys present in cfg are written, organized into sections for known parameters.
Any extra keys not belonging to a known section are written under [Other].
"""
function write_cfg(cfg::Dict{String, String}, path::String)
    # Known keys grouped by section (defines order within each section)
    input_keys   = ["resistance_file", "reclass_table", "project_name"]
    param_keys   = ["mode", "search_radius", "num_spokes", "injected_current", "spoke_aggregation"]
    compute_keys = ["solver", "connect_four_neighbors_only", "mask_nodata",
                    "nodata_value", "parallelize", "parallel_batch_size", "verbose"]
    output_keys  = ["write_as_tif"]

    known_keys = Set(vcat(input_keys, param_keys, compute_keys, output_keys))
    extra_keys = sort([k for k in keys(cfg) if k ∉ known_keys])

    open(path, "w") do f
        write(f, "[Input files]\n")
        for k in input_keys
            write(f, "$k = $(get(cfg, k, ""))\n")
        end
        write(f, "\n[MCI Parameters]\n")
        for k in param_keys
            write(f, "$k = $(get(cfg, k, ""))\n")
        end
        write(f, "\n[Computation]\n")
        for k in compute_keys
            write(f, "$k = $(get(cfg, k, ""))\n")
        end
        write(f, "\n[Output]\n")
        for k in output_keys
            write(f, "$k = $(get(cfg, k, ""))\n")
        end
        if !isempty(extra_keys)
            write(f, "\n[Other]\n")
            for k in extra_keys
                write(f, "$k = $(cfg[k])\n")
            end
        end
    end
end
