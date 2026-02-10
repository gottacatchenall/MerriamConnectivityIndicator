# ==============================================================================
# CONSTANTS
# ==============================================================================

const TRUELIST = ["true", "True", "TRUE"]

# ==============================================================================
# TYPE DEFINITIONS
# ==============================================================================

"""
    MCIConfig(; search_radius, num_spokes=8, ...)

Configuration for a Merriam Connectivity Indicator analysis.

# Required keyword argument
- `search_radius::Int64`: Radius of the circular analysis window in pixels.

# Optional keyword arguments
- `num_spokes::Int64 = 8`: Number of evenly-spaced source points on the circumference.
- `injected_current::Float64 = 1.0`: Total current distributed among valid spokes.
- `solver_type::String = "cg+amg"`: Circuitscape solver — `"cg+amg"` or `"cholmod"`.
- `connect_four_neighbors::Bool = false`: Use 4-neighbor connectivity (vs 8-neighbor).
- `mask_nodata::Bool = true`: Set MCI to NaN where the input raster is nodata.
- `nodata_value::Float64 = -9999.0`: Value that encodes missing data in the raster.
- `parallel_batch_size::Int64 = 100`: Windows per batch in threaded mode.
"""
struct MCIConfig
    search_radius::Int64
    num_spokes::Int64
    injected_current::Float64
    solver_type::String
    connect_four_neighbors::Bool
    mask_nodata::Bool
    nodata_value::Float64
    parallel_batch_size::Int64
end

function MCIConfig(;
    search_radius::Int64,
    num_spokes::Int64 = 8,
    injected_current::Float64 = 1.0,
    solver_type::String = "cg+amg",
    connect_four_neighbors::Bool = false,
    mask_nodata::Bool = true,
    nodata_value::Float64 = -9999.0,
    parallel_batch_size::Int64 = 100
)
    return MCIConfig(
        search_radius, num_spokes, injected_current, solver_type,
        connect_four_neighbors, mask_nodata, nodata_value, parallel_batch_size
    )
end

"""
    SpokePoint(row, col, angle)

A source point on the window circumference at grid position `(row, col)`,
located at `angle` radians from the center.
"""
struct SpokePoint
    row::Int64
    col::Int64
    angle::Float64
end
