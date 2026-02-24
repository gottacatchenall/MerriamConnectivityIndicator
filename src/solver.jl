# ==============================================================================
# CIRCUITSCAPE CIRCUIT SOLVER — ADVANCED MODE
# ==============================================================================

"""
    solve_window_circuit(conductance, source, ground; solver_type, use_four_neighbors) -> T

Solve a single MCI window using Circuitscape's Advanced mode pipeline and return
the effective resistance.

The grids are packaged into Circuitscape's `RasterData` / `RasterFlags` structs,
then processed through `compute_advanced_data` (graph construction + connected
component detection) and `multiple_solver` (linear solve).

Effective resistance is calculated via **McRae's energy formula**:

    R_eff = P / I² = Σ(sᵢ · Vᵢ) / I_total²

where `sᵢ` is the source current at node i and `Vᵢ` is its voltage (ground V=0).
"""
function solve_window_circuit(
    conductance::Array{T, 2},
    source::Array{T, 2},
    ground::Array{T, 2};
    solver_type::String = "cg+amg",
    use_four_neighbors::Bool = false
)::T where T <: AbstractFloat

    V = Int64

    # --- Package grids into Circuitscape data structures ---
    metadata = Circuitscape.RasterMeta(
        size(conductance, 2),             # ncols
        size(conductance, 1),             # nrows
        0.0, 0.0, 1.0,                    # georef placeholders (xll, yll, cellsize)
        -9999.0,                          # nodata sentinel
        Array{T, 1}(undef, 1),            # bins (unused)
        ""                                # projection (unused)
    )

    raster_data = Circuitscape.RasterData(
        conductance,                                        # resistance/conductance grid
        Matrix{V}(undef, 0, 0),                             # polymap (unused)
        source,                                             # source grid
        ground,                                             # ground grid
        (V[], V[], V[]),                                    # points_rc (unused)
        Matrix{T}(undef, 0, 0),                             # strengths (unused)
        Circuitscape.IncludeExcludePairs(
            :undef, V[], Matrix{V}(undef, 0, 0)
        ),                                                  # pairs (unused)
        metadata
    )

    cs_cfg = Dict{String, String}(
        "solver"                       => solver_type,
        "connect_four_neighbors_only"  => use_four_neighbors ? "true" : "false",
        "suppress_messages"            => "true",
        "data_type"                    => "int",
        "cholmod_batch_size"           => "1000",
    )

    flags = Circuitscape.RasterFlags(
        true,                           # is_raster
        false,                          # is_pairwise (advanced mode, not pairwise)
        true,                           # is_advanced
        false,                          # is_onetoall
        false,                          # is_alltoone
        false,                          # grnd_file_is_res
        Symbol("rmvsrc"),               # policy — remove sources co-located with ground
        use_four_neighbors,             # four_neighbors
        false,                          # avg_res (false = use conductance averaging)
        solver_type,                    # solver
        Circuitscape.OutputFlags(       # outputflags — no file output
            false, false, false, false,
            false, false, false, false
        )
    )

    # --- Build circuit graph and detect connected components ---
    data = Circuitscape.compute_advanced_data(raster_data, flags, cs_cfg)

    # --- Solve each connected component and accumulate power ---
    P_total = zero(T)   # total dissipated power: Σ(sᵢ · Vᵢ)
    I_total = zero(T)   # total current entering the circuit

    for cc_nodes in data.cc
        # Skip components that don't include the check node (if set)
        if data.check_node != -1 && !(data.check_node in cc_nodes)
            continue
        end

        # Extract sub-problem for this connected component
        G_local = data.G[cc_nodes, cc_nodes]
        s_local = data.sources[cc_nodes]
        g_local = data.grounds[cc_nodes]

        # Skip components with no sources or no grounds
        if sum(s_local) == 0 || sum(g_local) == 0
            continue
        end

        f_local = data.finitegrounds != [-9999.0] ?
            data.finitegrounds[cc_nodes] : data.finitegrounds

        # Solve for voltages (Circuitscape handles infinite-ground removal internally)
        voltages = Circuitscape.multiple_solver(
            cs_cfg, data.solver, G_local, s_local, g_local, f_local
        )

        # Map voltages back to the 2D grid to compute power
        local_nodemap = Circuitscape.construct_local_node_map(
            data.nodemap, cc_nodes, data.polymap
        )
        voltage_grid = zeros(T, size(source))
        @inbounds for j in 1:size(local_nodemap, 2), i in 1:size(local_nodemap, 1)
            idx = local_nodemap[i, j]
            if idx > 0
                voltage_grid[i, j] = voltages[idx]
            end
        end

        # Power = Σ(source current × voltage) at each node.
        # Ground nodes have V=0, so they contribute zero power.
        P_total += sum(source .* voltage_grid)
        I_total += sum(s_local)
    end

    # R_eff = P / I²  (McRae's energy formula)
    if I_total == 0
        return T(NaN)
    end
    return P_total / (I_total^2)
end

# ==============================================================================
# GRAPH CONSTRUCTION FROM CONDUCTANCE GRID
# ==============================================================================

"""
    build_nodemap(conductance) -> Matrix{Int64}

Assign sequential node IDs to each cell with positive conductance.
Cells with zero conductance get node ID 0 (no node).
"""
function build_nodemap(conductance::Matrix{T}) where T
    nodemap = zeros(Int64, size(conductance))
    idx = 1
    for j in 1:size(conductance, 2), i in 1:size(conductance, 1)
        if conductance[i, j] > 0
            nodemap[i, j] = idx
            idx += 1
        end
    end
    return nodemap
end

"""
    build_edge_list(conductance, nodemap; four_neighbors) -> (Vector{Int64}, Vector{Int64}, Vector{T})

Build a directed edge list from a conductance grid. Each edge appears once
(from lower to higher node index direction). The caller (Circuitscape's
`compute_graph_data`) symmetrizes the graph via `A + A'`.

Edge weights use conductance averaging to match Circuitscape's `construct_graph`:
- Orthogonal neighbors: `(C_a + C_b) / 2`
- Diagonal neighbors:   `(C_a + C_b) / (2√2)`  (extra √2 for longer path)
"""
function build_edge_list(
    conductance::Matrix{T},
    nodemap::Matrix{Int64};
    four_neighbors::Bool = false
) where T
    I = Int64[]
    J = Int64[]
    V = T[]
    nrows, ncols = size(conductance)

    for j in 1:ncols, i in 1:nrows
        nodemap[i, j] == 0 && continue

        # Right neighbor (horizontal)
        if j < ncols && nodemap[i, j + 1] != 0
            push!(I, nodemap[i, j])
            push!(J, nodemap[i, j + 1])
            push!(V, (conductance[i, j] + conductance[i, j + 1]) / 2)
        end

        # Bottom neighbor (vertical)
        if i < nrows && nodemap[i + 1, j] != 0
            push!(I, nodemap[i, j])
            push!(J, nodemap[i + 1, j])
            push!(V, (conductance[i, j] + conductance[i + 1, j]) / 2)
        end

        if !four_neighbors
            # Bottom-right diagonal
            if i < nrows && j < ncols && nodemap[i + 1, j + 1] != 0
                push!(I, nodemap[i, j])
                push!(J, nodemap[i + 1, j + 1])
                push!(V, (conductance[i, j] + conductance[i + 1, j + 1]) / (2 * √2))
            end

            # Top-right diagonal
            if i > 1 && j < ncols && nodemap[i - 1, j + 1] != 0
                push!(I, nodemap[i, j])
                push!(J, nodemap[i - 1, j + 1])
                push!(V, (conductance[i, j] + conductance[i - 1, j + 1]) / (2 * √2))
            end
        end
    end

    return (I, J, V)
end

# ==============================================================================
# CIRCUITSCAPE CIRCUIT SOLVER — PAIRWISE MODE (NETWORK PATH)
# ==============================================================================

"""
    solve_window_pairwise(conductance, center_row, center_col, valid_positions; ...) -> T

Solve a single MCI window using Circuitscape's **network pairwise** pipeline and
return the mean effective resistance between each spoke and the center.

The conductance grid is converted to a network (edge list + focal node IDs),
then solved via `compute_graph_data` → `single_ground_all_pairs`, which computes
effective resistances between all focal-point pairs.

Returns the mean R_eff across all spoke→center pairs, or `NaN` if no valid pairs.
"""
function solve_window_pairwise(
    conductance::Array{T, 2},
    center_row::Int64,
    center_col::Int64,
    valid_positions::Vector{Tuple{Int, Int}};
    solver_type::String = "cg+amg",
    use_four_neighbors::Bool = false
)::T where T <: AbstractFloat

    # --- Convert conductance grid to network representation ---
    nodemap = build_nodemap(conductance)
    edges = build_edge_list(conductance, nodemap; four_neighbors = use_four_neighbors)

    # --- Map focal points (center + spokes) to graph node IDs ---
    center_node = nodemap[center_row, center_col]
    n_spokes = length(valid_positions)
    fp = Vector{Int64}(undef, 1 + n_spokes)
    fp[1] = center_node
    for (k, (sr, sc)) in enumerate(valid_positions)
        fp[k + 1] = nodemap[sr, sc]
    end

    # --- Build NetworkData and solve via Circuitscape's network pairwise path ---
    networkdata = Circuitscape.NetworkData(
        edges,                            # (I, J, V) edge list
        fp,                               # focal point node IDs
        Matrix{T}(undef, 0, 0),          # source_map (unused in pairwise)
        Matrix{T}(undef, 0, 0)           # ground_map (unused in pairwise)
    )

    output_file = tempname()
    cs_cfg = Dict{String, String}(
        "solver"             => solver_type,
        "suppress_messages"  => "true",
        "cholmod_batch_size" => "1000",
        "output_file"        => output_file,
        "parallelize"        => "false",
        "data_type"          => "network",
    )

    # is_raster=true triggers Circuitscape's shortcut optimization in solve():
    # n solves instead of n(n-1)/2 when all output is disabled. The shortcut
    # math is purely algebraic and does not depend on raster data.
    flags = Circuitscape.NetworkFlags(
        true,                           # is_raster (enables shortcut optimization)
        false,                          # is_advanced
        false,                          # is_alltoone
        false,                          # is_onetoall
        false,                          # grnd_file_is_res
        Symbol("rmvsrc"),               # policy
        solver_type,                    # solver
        Circuitscape.OutputFlags(       # outputflags — all output disabled
            false, false, false, false,
            false, false, false, false
        )
    )

    graphdata = Circuitscape.compute_graph_data(networkdata, cs_cfg)

    local r
    try
        r = Circuitscape.single_ground_all_pairs(graphdata, flags, cs_cfg)
    finally
        for ext in ["_resistances.out", "_resistances_3columns.out"]
            f = output_file * ext
            isfile(f) && rm(f)
        end
    end

    # --- Extract spoke→center resistances from padded result matrix ---
    # r layout: row 1 = [0, fp1, fp2, ...], col 1 = [0, fp1, fp2, ...]
    # Center is at index 2 (first focal point); spoke k is at index k+2
    spoke_resistances = T[]
    for k in 1:n_spokes
        r_val = r[2, k + 2]   # r[center_index, spoke_k_index]
        if r_val >= 0          # -1 = not computed (different connected components)
            push!(spoke_resistances, T(r_val))
        end
    end

    return isempty(spoke_resistances) ? T(NaN) : mean(spoke_resistances)
end
