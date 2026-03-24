# ==============================================================================
# DIRECT SOLVER MODE — BYPASS CIRCUITSCAPE PAIRWISE FRAMEWORK
#
# Computes individual spoke→center effective resistances by:
#   1. Building the graph Laplacian directly from the conductance grid
#   2. Grounding the center node (removing its row/col)
#   3. Factorizing the reduced Laplacian once (Cholesky)
#   4. Batch back-substituting for all spokes in one call
#
# This avoids all Circuitscape overhead: no NetworkData/RasterData construction,
# no n×n pairwise resistance matrices, no temp file I/O, no shortcut algebra.
# ==============================================================================

# ==============================================================================
# LAPLACIAN CONSTRUCTION
# ==============================================================================

"""
    build_laplacian(conductance, nodemap; four_neighbors) -> SparseMatrixCSC

Build the graph Laplacian directly from a conductance grid and nodemap.

Edge weights use conductance averaging to match Circuitscape's `construct_graph`:
- Orthogonal neighbors: `(C_a + C_b) / 2`
- Diagonal neighbors:   `(C_a + C_b) / (2√2)`

The Laplacian L has `L[i,i] = sum of edge weights at node i` and
`L[i,j] = -edge_weight(i,j)` for adjacent nodes.
"""
function build_laplacian(
    conductance::Matrix{T},
    nodemap::Matrix{Int64};
    four_neighbors::Bool = false
) where T

    n_nodes = maximum(nodemap)
    nrows, ncols = size(conductance)

    # Pre-allocate COO vectors (4 entries per edge: 2 off-diag + 2 diag contributions)
    max_edges = 4 * n_nodes  # upper bound for 8-neighbor
    I_idx = Int64[]
    J_idx = Int64[]
    V_val = T[]
    sizehint!(I_idx, max_edges * 4)
    sizehint!(J_idx, max_edges * 4)
    sizehint!(V_val, max_edges * 4)

    inv_sqrt2 = T(1 / √2)

    @inbounds for j in 1:ncols, i in 1:nrows
        nodemap[i, j] == 0 && continue
        ni = nodemap[i, j]

        # Right neighbor
        if j < ncols && nodemap[i, j + 1] != 0
            nj = nodemap[i, j + 1]
            w = (conductance[i, j] + conductance[i, j + 1]) / 2
            push!(I_idx, ni); push!(J_idx, nj); push!(V_val, -w)
            push!(I_idx, nj); push!(J_idx, ni); push!(V_val, -w)
            push!(I_idx, ni); push!(J_idx, ni); push!(V_val, w)
            push!(I_idx, nj); push!(J_idx, nj); push!(V_val, w)
        end

        # Bottom neighbor
        if i < nrows && nodemap[i + 1, j] != 0
            nj = nodemap[i + 1, j]
            w = (conductance[i, j] + conductance[i + 1, j]) / 2
            push!(I_idx, ni); push!(J_idx, nj); push!(V_val, -w)
            push!(I_idx, nj); push!(J_idx, ni); push!(V_val, -w)
            push!(I_idx, ni); push!(J_idx, ni); push!(V_val, w)
            push!(I_idx, nj); push!(J_idx, nj); push!(V_val, w)
        end

        if !four_neighbors
            # Bottom-right diagonal
            if i < nrows && j < ncols && nodemap[i + 1, j + 1] != 0
                nj = nodemap[i + 1, j + 1]
                w = (conductance[i, j] + conductance[i + 1, j + 1]) / 2 * inv_sqrt2
                push!(I_idx, ni); push!(J_idx, nj); push!(V_val, -w)
                push!(I_idx, nj); push!(J_idx, ni); push!(V_val, -w)
                push!(I_idx, ni); push!(J_idx, ni); push!(V_val, w)
                push!(I_idx, nj); push!(J_idx, nj); push!(V_val, w)
            end

            # Top-right diagonal
            if i > 1 && j < ncols && nodemap[i - 1, j + 1] != 0
                nj = nodemap[i - 1, j + 1]
                w = (conductance[i, j] + conductance[i - 1, j + 1]) / 2 * inv_sqrt2
                push!(I_idx, ni); push!(J_idx, nj); push!(V_val, -w)
                push!(I_idx, nj); push!(J_idx, ni); push!(V_val, -w)
                push!(I_idx, ni); push!(J_idx, ni); push!(V_val, w)
                push!(I_idx, nj); push!(J_idx, nj); push!(V_val, w)
            end
        end
    end

    sparse(I_idx, J_idx, V_val, n_nodes, n_nodes)
end

# ==============================================================================
# DIRECT WINDOW SOLVER
# ==============================================================================

"""
    solve_window_direct(conductance, center_row, center_col, valid_positions; ...) -> T

Solve a single MCI window by direct Cholesky factorization of the graph Laplacian.

Grounds the center node (removes its row/col from the Laplacian), then injects
1A at each spoke independently via batched back-substitution. Returns the median
(default) or mean spoke→center effective resistance.

This is equivalent to the pairwise mode result but avoids all Circuitscape
framework overhead.
"""
function solve_window_direct(
    conductance::Array{T, 2},
    center_row::Int64,
    center_col::Int64,
    valid_positions::Vector{Tuple{Int, Int}};
    use_four_neighbors::Bool = false,
    aggregation::Symbol = :median
)::T where T <: AbstractFloat

    # Build nodemap (reuse existing function from solver.jl)
    nodemap = build_nodemap(conductance)
    n_nodes = maximum(nodemap)

    n_nodes == 0 && return T(NaN)

    # Build graph Laplacian
    L = build_laplacian(conductance, nodemap; four_neighbors = use_four_neighbors)

    # Ground center: remove its row/col from the Laplacian
    center_node = nodemap[center_row, center_col]
    center_node == 0 && return T(NaN)

    keep = [i for i in 1:n_nodes if i != center_node]
    n_reduced = length(keep)

    # Map old node indices to new (post-removal) indices
    new_idx = zeros(Int64, n_nodes)
    for (new_id, old_id) in enumerate(keep)
        new_idx[old_id] = new_id
    end

    L_reduced = L[keep, keep]

    # Small regularization for numerical stability (matches Circuitscape's CHOLMOD path)
    L_reduced += spdiagm(0 => fill(T(10) * eps(T), n_reduced))

    # Map spoke positions to reduced node indices
    n_spokes = length(valid_positions)
    spoke_reduced = Vector{Int64}(undef, n_spokes)
    for (k, (r, c)) in enumerate(valid_positions)
        spoke_reduced[k] = new_idx[nodemap[r, c]]
    end

    # Build RHS matrix: column k has 1.0 at spoke k's node (1A injection)
    RHS = sparse(spoke_reduced, 1:n_spokes, ones(T, n_spokes), n_reduced, n_spokes)

    # Factorize once, solve all spokes in one batched back-substitution
    F = cholesky(L_reduced)
    V = F \ RHS

    # R_eff for spoke k = voltage at spoke k's node (center is grounded at V=0)
    spoke_resistances = T[]
    for k in 1:n_spokes
        r = V[spoke_reduced[k], k]
        if r >= 0 && isfinite(r)
            push!(spoke_resistances, r)
        end
    end

    isempty(spoke_resistances) && return T(NaN)
    return aggregation === :median ? median(spoke_resistances) : mean(spoke_resistances)
end

# ==============================================================================
# SINGLE-WINDOW MCI (DIRECT MODE)
# ==============================================================================

"""
    compute_mci_direct_for_window(resistance_full, row, col, config) -> T

Compute the MCI value for one window using direct Cholesky factorization.
Returns `NaN` if the center is invalid, no valid spokes exist, or the solve fails.
"""
function compute_mci_direct_for_window(
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
        return solve_window_direct(
            w.conductance, w.center_row, w.center_col, w.valid_positions;
            use_four_neighbors = config.connect_four_neighbors,
            aggregation = config.spoke_aggregation
        )
    catch e
        @warn "Direct MCI solve failed at ($center_row_global, $center_col_global): $e"
        return T(NaN)
    end
end

# ==============================================================================
# SLIDING-WINDOW ANALYSIS (DIRECT MODE)
# ==============================================================================

"""
    compute_mci_direct(resistance, config; parallelize, verbose) -> Matrix

Compute the Merriam Connectivity Indicator for every valid pixel using direct
Cholesky factorization, bypassing Circuitscape's pairwise framework.

For each window, the graph Laplacian is factorized once, then all spoke→center
effective resistances are computed via batched back-substitution. The pixel's
MCI value is the mean spoke→center R_eff (each spoke injects 1A independently).

This produces the same result as `compute_mci_pairwise` but avoids the overhead
of Circuitscape's pairwise dispatch, n×n resistance matrices, shortcut algebra,
and temporary file I/O.
"""
function compute_mci_direct(
    resistance::Array{T, 2},
    config::MCIConfig;
    parallelize::Bool = (nthreads() > 1),
    verbose::Bool = true
)::Array{T, 2} where T <: AbstractFloat

    grid_height, grid_width = size(resistance)
    mci_raster = fill(T(NaN), grid_height, grid_width)

    if parallelize && nthreads() > 1
        BLAS.set_num_threads(1)
        verbose && @info "Computing direct MCI with $(nthreads()) threads"
    else
        verbose && @info "Computing direct MCI with 1 thread"
    end

    # Collect valid center pixels
    valid_centers = Tuple{Int64, Int64}[]
    for col in 1:grid_width, row in 1:grid_height
        r = resistance[row, col]
        if r != config.nodata_value && !isnan(r) && r > 0
            push!(valid_centers, (row, col))
        end
    end
    num_windows = length(valid_centers)

    verbose && @info "Processing $num_windows windows direct (radius=$(config.search_radius), spokes=$(config.num_spokes))"
    progress = verbose ? Progress(num_windows; dt = 0.25) : nothing

    if parallelize && nthreads() > 1
        batch_size = config.parallel_batch_size
        num_batches = cld(num_windows, batch_size)

        @threads for batch_idx in 0:(num_batches - 1)
            lo = batch_idx * batch_size + 1
            hi = min(num_windows, lo + batch_size - 1)
            for idx in lo:hi
                row, col = valid_centers[idx]
                mci_raster[row, col] = compute_mci_direct_for_window(
                    resistance, row, col, config
                )
                progress !== nothing && next!(progress)
            end
        end
    else
        for (row, col) in valid_centers
            mci_raster[row, col] = compute_mci_direct_for_window(
                resistance, row, col, config
            )
            progress !== nothing && next!(progress)
        end
    end

    if config.mask_nodata
        mci_raster[resistance .== config.nodata_value] .= T(NaN)
        mci_raster[isnan.(resistance)] .= T(NaN)
    end

    verbose && @info "Direct MCI computation complete"
    return mci_raster
end
