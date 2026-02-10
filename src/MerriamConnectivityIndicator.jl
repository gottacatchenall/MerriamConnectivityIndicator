"""
    MerriamConnectivityIndicator

Measures landscape connectivity using circuit theory via a sliding-window approach.

Each window places a **ground** at the center pixel (voltage = 0) and injects current
from evenly-spaced **spoke** source points on the circumference. The resulting effective
resistance (MCI value) is computed using:

    R_eff = P / I² = Σ(sᵢ · Vᵢ) / I_total²

where sᵢ is the current injected at spoke i and Vᵢ is its solved voltage. When the
total injected current is normalized to 1.0, R_eff equals the average spoke voltage.

Circuit solving is delegated to Circuitscape, which handles graph construction,
connected-component detection, and linear system solving internally. Two modes
are available:

- **Advanced mode** (`compute_mci`): All spokes inject current simultaneously into
  a single ground at the center. One solve per window.
- **Pairwise mode** (`compute_mci_pairwise`): The center and spokes are passed as
  focal nodes to Circuitscape's native pairwise solver, which computes effective
  resistance between all pairs. The MCI value is the mean spoke→center R_eff.
"""
module MerriamConnectivityIndicator

    using SparseArrays
    using LinearAlgebra.BLAS
    using ArchGDAL
    using Base.Threads
    using LinearAlgebra
    using Circuitscape
    using DelimitedFiles
    using ProgressMeter
    using Random
    using StatsBase
    using Statistics
    using Omniscape

    include("types.jl")
    include("spokes.jl")
    include("window.jl")
    include("solver.jl")
    include("mci.jl")

    # ==============================================================================
    # EXPORTS
    # ==============================================================================

    export MCIConfig, SpokePoint
    export compute_mci, compute_mci_pairwise, compute_mci_from_dict
    export generate_spoke_points, build_window_grids, solve_window_pairwise

end # module MerriamConnectivityIndicator
