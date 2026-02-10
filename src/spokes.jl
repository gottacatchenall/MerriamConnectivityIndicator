"""
    generate_spoke_points(center_row, center_col, radius, num_spokes) -> Vector{SpokePoint}

Place `num_spokes` evenly-spaced points on a circle of `radius` pixels centered at
`(center_row, center_col)`. Positions are rounded to integer grid coordinates.

Angle 0 points east (increasing column); row increases downward, so the sin
component is negated.
"""
function generate_spoke_points(
    center_row::Int64,
    center_col::Int64,
    radius::Int64,
    num_spokes::Int64
)::Vector{SpokePoint}
    angle_increment = 2π / num_spokes
    spokes = Vector{SpokePoint}(undef, num_spokes)
    for i in 0:(num_spokes - 1)
        angle = i * angle_increment
        spoke_row = center_row - round(Int64, radius * sin(angle))
        spoke_col = center_col + round(Int64, radius * cos(angle))
        spokes[i + 1] = SpokePoint(spoke_row, spoke_col, angle)
    end
    return spokes
end
