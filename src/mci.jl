# ==============================================================================
# SINGLE-WINDOW MCI COMPUTATION
# ==============================================================================

"""
    compute_mci_for_window(resistance_full, row, col, config) -> T

Compute the MCI value for one window centered at global position `(row, col)`.

Returns `NaN` if the center is invalid, no valid spokes exist, or the solve fails.
"""
function compute_mci_for_window(
    resistance_full::Array{T, 2},
    center_row_global::Int64,
    center_col_global::Int64,
    config::MCIConfig
)::T where T <: AbstractFloat

    grids = build_window_grids(
        resistance_full, center_row_global, center_col_global, config
    )

    # build_window_grids returns nothing for invalid centers or no valid spokes
    if grids === nothing
        return T(NaN)
    end

    conductance, source, ground = grids

    try
        return solve_window_circuit(
            conductance, source, ground;
            solver_type = config.solver_type,
            use_four_neighbors = config.connect_four_neighbors
        )
    catch e
        @warn "MCI solve failed at ($center_row_global, $center_col_global): $e"
        return T(NaN)
    end
end

# ==============================================================================
# MAIN SLIDING-WINDOW ANALYSIS
# ==============================================================================

"""
    compute_mci(resistance, config; parallelize, verbose) -> Matrix

Compute the Merriam Connectivity Indicator for every valid pixel in a resistance
raster using a sliding circular window.

Each pixel's MCI value is the effective resistance from circumference spokes to the
center, computed by Circuitscape's Advanced mode. Higher values indicate *lower*
connectivity (more resistance to current flow).

# Arguments
- `resistance::Matrix{T}`: Resistance raster (nodata as NaN or `config.nodata_value`).
- `config::MCIConfig`: Analysis configuration.

# Keyword arguments
- `parallelize::Bool`: Use multi-threaded processing (default: true if threads > 1).
- `verbose::Bool`: Print progress information (default: true).

# Returns
A `Matrix{T}` of the same size as `resistance`, with MCI values at valid pixels
and `NaN` elsewhere.
"""
function compute_mci(
    resistance::Array{T, 2},
    config::MCIConfig;
    parallelize::Bool = (nthreads() > 1),
    verbose::Bool = true
)::Array{T, 2} where T <: AbstractFloat

    grid_height, grid_width = size(resistance)
    mci_raster = fill(T(NaN), grid_height, grid_width)

    # Restrict BLAS to 1 thread per call when using Julia-level threading
    if parallelize && nthreads() > 1
        BLAS.set_num_threads(1)
        verbose && @info "Computing MCI with $(nthreads()) threads"
    else
        verbose && @info "Computing MCI with 1 thread"
    end

    # Collect all valid center pixels (positive, non-NaN, non-nodata resistance)
    valid_centers = Tuple{Int64, Int64}[]
    for col in 1:grid_width, row in 1:grid_height
        r = resistance[row, col]
        if r != config.nodata_value && !isnan(r) && r > 0
            push!(valid_centers, (row, col))
        end
    end
    num_windows = length(valid_centers)

    verbose && @info "Processing $num_windows windows (radius=$(config.search_radius), spokes=$(config.num_spokes))"
    progress = verbose ? Progress(num_windows; dt = 0.25) : nothing

    if parallelize && nthreads() > 1
        batch_size = config.parallel_batch_size
        num_batches = cld(num_windows, batch_size)

        @threads for batch_idx in 0:(num_batches - 1)
            lo = batch_idx * batch_size + 1
            hi = min(num_windows, lo + batch_size - 1)
            for idx in lo:hi
                row, col = valid_centers[idx]
                mci_raster[row, col] = compute_mci_for_window(
                    resistance, row, col, config
                )
                progress !== nothing && next!(progress)
            end
        end
    else
        for (row, col) in valid_centers
            mci_raster[row, col] = compute_mci_for_window(
                resistance, row, col, config
            )
            progress !== nothing && next!(progress)
        end
    end

    # Mask nodata pixels in the output
    if config.mask_nodata
        mci_raster[resistance .== config.nodata_value] .= T(NaN)
        mci_raster[isnan.(resistance)] .= T(NaN)
    end

    verbose && @info "MCI computation complete"
    return mci_raster
end

# ==============================================================================
# PAIRWISE MCI MODE — SINGLE-WINDOW
# ==============================================================================

"""
    compute_mci_pairwise_for_window(resistance_full, row, col, config) -> T

Compute the pairwise MCI value for one window centered at `(row, col)`.

Uses Circuitscape's native pairwise mode: the center and all valid spoke
positions are passed as focal nodes, and Circuitscape computes effective
resistances between all pairs in a single call. The MCI value is the mean
spoke→center R_eff.

Returns `NaN` if the center is invalid, no valid spokes exist, or the solve fails.
"""
function compute_mci_pairwise_for_window(
    resistance_full::Array{T, 2},
    center_row_global::Int64,
    center_col_global::Int64,
    config::MCIConfig
)::T where T <: AbstractFloat

    w = prepare_window(resistance_full, center_row_global, center_col_global, config)
    if w === nothing
        return T(NaN)
    end

    try
        return solve_window_pairwise(
            w.conductance, w.center_row, w.center_col, w.valid_positions;
            solver_type = config.solver_type,
            use_four_neighbors = config.connect_four_neighbors
        )
    catch e
        @warn "Pairwise MCI solve failed at ($center_row_global, $center_col_global): $e"
        return T(NaN)
    end
end

# ==============================================================================
# PAIRWISE MCI MODE — SLIDING WINDOW
# ==============================================================================

"""
    compute_mci_pairwise(resistance, config; parallelize, verbose) -> Matrix

Compute the **pairwise** Merriam Connectivity Indicator for every valid pixel.

Uses Circuitscape's native pairwise mode: for each window, the center and all
valid spoke positions are passed as focal nodes. Circuitscape computes effective
resistances between all pairs using an efficient shortcut (n solves instead of
n(n-1)/2), and the pixel's MCI value is the mean spoke→center R_eff.

"""
function compute_mci_pairwise(
    resistance::Array{T, 2},
    config::MCIConfig;
    parallelize::Bool = (nthreads() > 1),
    verbose::Bool = true
)::Array{T, 2} where T <: AbstractFloat

    grid_height, grid_width = size(resistance)
    mci_raster = fill(T(NaN), grid_height, grid_width)

    if parallelize && nthreads() > 1
        BLAS.set_num_threads(1)
        verbose && @info "Computing pairwise MCI with $(nthreads()) threads"
    else
        verbose && @info "Computing pairwise MCI with 1 thread"
    end

    valid_centers = Tuple{Int64, Int64}[]
    for col in 1:grid_width, row in 1:grid_height
        r = resistance[row, col]
        if r != config.nodata_value && !isnan(r) && r > 0
            push!(valid_centers, (row, col))
        end
    end
    num_windows = length(valid_centers)

    verbose && @info "Processing $num_windows windows pairwise (radius=$(config.search_radius), spokes=$(config.num_spokes))"
    progress = verbose ? Progress(num_windows; dt = 0.25) : nothing

    if parallelize && nthreads() > 1
        batch_size = config.parallel_batch_size
        num_batches = cld(num_windows, batch_size)

        @threads for batch_idx in 0:(num_batches - 1)
            lo = batch_idx * batch_size + 1
            hi = min(num_windows, lo + batch_size - 1)
            for idx in lo:hi
                row, col = valid_centers[idx]
                mci_raster[row, col] = compute_mci_pairwise_for_window(
                    resistance, row, col, config
                )
                progress !== nothing && next!(progress)
            end
        end
    else
        for (row, col) in valid_centers
            mci_raster[row, col] = compute_mci_pairwise_for_window(
                resistance, row, col, config
            )
            progress !== nothing && next!(progress)
        end
    end

    if config.mask_nodata
        mci_raster[resistance .== config.nodata_value] .= T(NaN)
        mci_raster[isnan.(resistance)] .= T(NaN)
    end

    verbose && @info "Pairwise MCI computation complete"
    return mci_raster
end

# ==============================================================================
# DICTIONARY-BASED CONVENIENCE WRAPPER
# ==============================================================================

"""
    compute_mci_from_dict(resistance, config_dict) -> Matrix

Compute MCI using a string-keyed dictionary for configuration, similar to the
Omniscape/Circuitscape INI-file interface.

# Required key
- `"search_radius"`: Window radius in pixels.

# Optional keys (with defaults)
- `"num_spokes"` (`"8"`), `"injected_current"` (`"1.0"`), `"solver"` (`"cg+amg"`),
  `"connect_four_neighbors_only"` (`"false"`), `"mask_nodata"` (`"true"`),
  `"nodata_value"` (`"-9999.0"`), `"parallel_batch_size"` (`"100"`),
  `"parallelize"` (`"true"`), `"verbose"` (`"true"`).
"""
function compute_mci_from_dict(
    resistance::Array{T, 2},
    config_dict::Dict{String, String}
)::Array{T, 2} where T <: AbstractFloat

    mci_config = MCIConfig(
        search_radius       = parse(Int64,   config_dict["search_radius"]),
        num_spokes          = parse(Int64,   get(config_dict, "num_spokes", "8")),
        injected_current    = parse(Float64, get(config_dict, "injected_current", "1.0")),
        solver_type         = get(config_dict, "solver", "cg+amg"),
        connect_four_neighbors = lowercase(get(config_dict, "connect_four_neighbors_only", "false")) in TRUELIST,
        mask_nodata         = lowercase(get(config_dict, "mask_nodata", "true")) in TRUELIST,
        nodata_value        = parse(Float64, get(config_dict, "nodata_value", "-9999.0")),
        parallel_batch_size = parse(Int64,   get(config_dict, "parallel_batch_size", "100"))
    )

    parallelize = lowercase(get(config_dict, "parallelize", "true")) in TRUELIST
    verbose     = lowercase(get(config_dict, "verbose", "true")) in TRUELIST

    return compute_mci(resistance, mci_config; parallelize, verbose)
end
