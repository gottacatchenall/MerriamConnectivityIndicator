# ==============================================================================
# WINDOW GRID CONSTRUCTION
# ==============================================================================

# Reflect index i into the valid range 1:n (mirror padding at boundaries).
# Examples for n=5: 0→2, -1→3, 6→4, 7→3
@inline function _mirror(i::Int, n::Int)::Int
    i = i < 1 ? 2 - i : (i > n ? 2 * n - i : i)
    return clamp(i, 1, n)  # safety for radius > raster size
end

# Extract a full (2r+1)×(2r+1) window centered at (cr, cc) using mirror
# padding so out-of-raster pixels reflect real landscape values rather than
# being treated as zero conductance.
function _extract_mirrored(
    resistance_full::Array{T, 2},
    cr_global::Int, cc_global::Int, radius::Int
)::Matrix{T} where T
    h, w = size(resistance_full)
    sz = 2 * radius + 1
    result = Matrix{T}(undef, sz, sz)
    for di in 1:sz
        ri = _mirror(cr_global - radius + di - 1, h)
        for dj in 1:sz
            cj = _mirror(cc_global - radius + dj - 1, w)
            result[di, dj] = resistance_full[ri, cj]
        end
    end
    return result
end

"""
    prepare_window(resistance_full, row, col, config) -> NamedTuple or nothing

Shared setup for both advanced and pairwise MCI modes. Extracts the resistance
subset around `(row, col)`, converts to conductance, applies the circular mask,
identifies valid spoke positions, and builds the ground grid.

The window is always a full (2r+1)×(2r+1) square. Pixels outside the raster
boundary are filled by mirror reflection so edge windows have the same network
size as interior windows, eliminating the edge-inflation artifact.

Returns a NamedTuple `(conductance, valid_positions, ground, center_row, center_col)`
or `nothing` if the center is invalid or no valid spokes exist.
"""
function prepare_window(
    resistance_full::Array{T, 2},
    center_row_global::Int64,
    center_col_global::Int64,
    config::MCIConfig
) where T <: AbstractFloat

    radius = config.search_radius

    # Always extract a full (2r+1)×(2r+1) window; mirror at boundaries.
    resistance_sub = _extract_mirrored(resistance_full,
                                       center_row_global, center_col_global, radius)

    # Center is always at the middle of the fixed-size window.
    cr = radius + 1
    cc = radius + 1

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
