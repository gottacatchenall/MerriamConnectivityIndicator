# ==============================================================================
# WINDOW GRID CONSTRUCTION
# ==============================================================================

"""
    prepare_window(resistance_full, row, col, config) -> NamedTuple or nothing

Shared setup for both advanced and pairwise MCI modes. Extracts the resistance
subset around `(row, col)`, converts to conductance, applies the circular mask,
identifies valid spoke positions, and builds the ground grid.

Returns a NamedTuple `(conductance, valid_positions, ground, center_row, center_col)`
or `nothing` if the center is invalid or no valid spokes exist.
"""
function prepare_window(
    resistance_full::Array{T, 2},
    center_row_global::Int64,
    center_col_global::Int64,
    config::MCIConfig
) where T <: AbstractFloat

    grid_height, grid_width = size(resistance_full)
    radius = config.search_radius

    # --- Extract rectangular subset around center ---
    row_lo = max(1, center_row_global - radius)
    row_hi = min(grid_height, center_row_global + radius)
    col_lo = max(1, center_col_global - radius)
    col_hi = min(grid_width, center_col_global + radius)

    resistance_sub = resistance_full[row_lo:row_hi, col_lo:col_hi]

    # Center position in local (subset) coordinates
    cr = center_row_global - row_lo + 1
    cc = center_col_global - col_lo + 1

    # --- Build conductance grid (zero out invalid cells) ---
    conductance = 1.0 ./ resistance_sub
    conductance[resistance_sub .== config.nodata_value] .= 0.0
    conductance[isnan.(resistance_sub)] .= 0.0
    conductance[resistance_sub .<= 0] .= 0.0
    conductance[isinf.(conductance)] .= 0.0

    # Apply circular mask — zero everything outside the radius
    sub_h, sub_w = size(conductance)
    for j in 1:sub_w, i in 1:sub_h
        if sqrt((i - cr)^2 + (j - cc)^2) > radius
            conductance[i, j] = 0.0
        end
    end

    # Center must be a valid (positive conductance) pixel
    if conductance[cr, cc] <= 0
        return nothing
    end

    # --- Identify valid spoke positions ---
    spokes = generate_spoke_points(cr, cc, radius, config.num_spokes)

    valid_positions = Tuple{Int, Int}[]
    for sp in spokes
        r, c = sp.row, sp.col
        if 1 <= r <= sub_h && 1 <= c <= sub_w && conductance[r, c] > 0
            push!(valid_positions, (r, c))
        end
    end

    if isempty(valid_positions)
        return nothing
    end

    # --- Build ground grid (center = infinite ground / Dirichlet V=0) ---
    ground = zeros(T, sub_h, sub_w)
    ground[cr, cc] = T(Inf)

    return (conductance = conductance, valid_positions = valid_positions, ground = ground,
            center_row = cr, center_col = cc)
end

"""
    build_window_grids(resistance_full, row, col, config) -> (conductance, source, ground) or nothing

Extract and prepare the three grids needed for one **advanced-mode** Circuitscape solve:

1. **conductance** — the resistance subset around `(row, col)` converted to conductance
   (1/R), with nodata and out-of-circle pixels zeroed.
2. **source** — current injected at each spoke pixel (total `config.injected_current`
   divided equally among valid spokes).
3. **ground** — `Inf` at the center pixel, zero elsewhere (Dirichlet boundary).

Returns `nothing` if the center pixel is invalid or no valid spokes are found.
"""
function build_window_grids(
    resistance_full::Array{T, 2},
    center_row_global::Int64,
    center_col_global::Int64,
    config::MCIConfig
) where T <: AbstractFloat

    w = prepare_window(resistance_full, center_row_global, center_col_global, config)
    if w === nothing
        return nothing
    end

    # Distribute total current equally among all valid spokes
    source = zeros(T, size(w.conductance))
    current_per_spoke = config.injected_current / length(w.valid_positions)
    for (r, c) in w.valid_positions
        source[r, c] += current_per_spoke
    end

    return w.conductance, source, w.ground
end
