"""
    run_mci(path::String) -> Matrix

Run an MCI analysis from an .ini configuration file.

Reads the resistance raster, optionally reclassifies it, computes MCI using
either advanced or pairwise mode, writes the result GeoTIFF, and saves a copy
of the configuration.
"""
function run_mci(path::String)
    cfg_user = parse_cfg(path)
    cfg = init_cfg()
    update_cfg!(cfg, cfg_user)

    # Validate required keys
    cfg["resistance_file"] == "" && error("resistance_file is required in the .ini file")
    cfg["search_radius"] == "" && error("search_radius is required in the .ini file")
    cfg["project_name"] == "" && error("project_name is required in the .ini file")

    verbose = lowercase(cfg["verbose"]) in TRUELIST
    verbose && @info "Reading resistance raster: $(cfg["resistance_file"])"

    # Read input raster
    resistance, wkt, transform = read_raster(cfg["resistance_file"], Float64)

    # Apply reclassification if specified
    if cfg["reclass_table"] != ""
        verbose && @info "Applying reclassification from: $(cfg["reclass_table"])"
        nodata_val = parse(Float64, cfg["nodata_value"])
        reclass = read_reclass_table(cfg["reclass_table"]; nodata_value=nodata_val)
        resistance = map(x -> get(reclass, x, x), resistance)
    end

    # Run MCI
    result = run_mci(cfg, resistance, wkt, transform)

    return result
end

"""
    run_mci(cfg::Dict{String,String}, resistance::Matrix, wkt::String, transform) -> Matrix

Run an MCI analysis from an in-memory configuration dictionary and pre-loaded
resistance data. This is the in-memory method, analogous to Omniscape's dual
interface.

The `wkt` and `transform` arguments are used when writing output GeoTIFFs.
"""
function run_mci(
    cfg::Dict{String, String},
    resistance::Array{T, 2},
    wkt::String,
    transform
) where T <: AbstractFloat

    # Build MCIConfig from dict (reuses same keys as compute_mci_from_dict)
    mci_config = MCIConfig(
        search_radius       = parse(Int64,   cfg["search_radius"]),
        num_spokes          = parse(Int64,   get(cfg, "num_spokes", "8")),
        injected_current    = parse(Float64, get(cfg, "injected_current", "1.0")),
        solver_type         = get(cfg, "solver", "cg+amg"),
        connect_four_neighbors = lowercase(get(cfg, "connect_four_neighbors_only", "false")) in TRUELIST,
        mask_nodata         = lowercase(get(cfg, "mask_nodata", "true")) in TRUELIST,
        nodata_value        = parse(Float64, get(cfg, "nodata_value", "-9999.0")),
        parallel_batch_size = parse(Int64,   get(cfg, "parallel_batch_size", "100")),
        spoke_aggregation   = Symbol(get(cfg, "spoke_aggregation", "median"))
    )

    parallelize = lowercase(get(cfg, "parallelize", "true")) in TRUELIST
    verbose     = lowercase(get(cfg, "verbose", "true")) in TRUELIST
    mode        = lowercase(get(cfg, "mode", "advanced"))

    # Compute MCI
    if mode == "pairwise"
        verbose && @info "Running MCI in pairwise mode"
        result = compute_mci_pairwise(resistance, mci_config;
                                      parallelize=parallelize, verbose=verbose)
    elseif mode == "pairwise-direct"
        verbose && @info "Running MCI in pairwise-direct mode"
        result = compute_mci_direct(resistance, mci_config;
                                      parallelize=parallelize, verbose=verbose)
    elseif mode == "advanced"
        verbose && @info "Running MCI in advanced mode"
        result = compute_mci(resistance, mci_config;
                             parallelize=parallelize, verbose=verbose)
    else
        error("Invalid mode: $mode")
        exit(-1)
    end

    # Write outputs
    project_name = get(cfg, "project_name", "")
    write_as_tif = lowercase(get(cfg, "write_as_tif", "true")) in TRUELIST

    if project_name != "" && write_as_tif
        # Create output directory
        mkpath(project_name)

        # Write result raster
        output_path = joinpath(project_name, "mci_$(mode).tif")
        verbose && @info "Writing result to: $output_path"
        write_raster(output_path, result, wkt, transform)

        # Write config copy for reproducibility
        config_path = joinpath(project_name, "config.ini")
        write_cfg(cfg, config_path)

        verbose && @info "Outputs written to $(abspath(project_name))"
    end

    return result
end
