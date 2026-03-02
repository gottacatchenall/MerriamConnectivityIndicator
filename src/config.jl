"""
    init_cfg() -> Dict{String,String}

Return a dictionary of all MCI configuration keys with their default values.
"""
function init_cfg()
    cfg = Dict{String, String}()

    # Input files
    cfg["resistance_file"] = ""
    cfg["reclass_table"] = ""
    cfg["project_name"] = ""

    # MCI parameters
    cfg["mode"] = "advanced"
    cfg["search_radius"] = ""
    cfg["num_spokes"] = "8"
    cfg["injected_current"] = "1.0"

    # Computation
    cfg["solver"] = "cg+amg"
    cfg["connect_four_neighbors_only"] = "false"
    cfg["mask_nodata"] = "true"
    cfg["nodata_value"] = "-9999.0"
    cfg["parallelize"] = "true"
    cfg["parallel_batch_size"] = "100"
    cfg["verbose"] = "true"

    # Output
    cfg["write_as_tif"] = "true"

    cfg
end

"""
    parse_cfg(path::String) -> Dict{String,String}

Parse an .ini file into a `Dict{String,String}`. Section headers (`[...]`),
blank lines, and comment lines (starting with `#`) are skipped.
"""
function parse_cfg(path::String)
    cf = Dict{String, String}()
    f = open(path, "r")
    for line in eachline(f, keep = true)
        stripped = strip(line)
        isempty(stripped) && continue
        first(stripped) == '[' && continue
        first(stripped) == '#' && continue
        idx = something(findfirst(isequal('='), stripped), 0)
        idx == 0 && continue
        var = rstrip(stripped[1:idx - 1])
        val = strip(stripped[idx + 1:end])
        cf[var] = val
    end
    close(f)
    cf
end

"""
    update_cfg!(cfg, cfg_new)

Merge user-supplied keys from `cfg_new` onto `cfg`, overwriting defaults.
"""
function update_cfg!(cfg, cfg_new)
    for (key, val) in cfg_new
        cfg[key] = val
    end
end
