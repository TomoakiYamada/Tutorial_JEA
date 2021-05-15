"""
grid_uni(grid_min::Float64, grid_max::Float64, num_grid::Int64)

Purpose:
Generate uniform grid from "grid_min" to "grid_max".
"linspace" in Matlab.

Input:
grid_min -> minimum of grid
grid_max -> maximum of grid
num_grid -> # grid

Output:
grid -> Vector of num_grid × 1
"""
function grid_uni(grid_min::Float64, grid_max::Float64, num_grid::Int64)
    grid = zeros(num_grid)
    increment = (grid_max - grid_min) / (num_grid-1)
    for i in 1:num_grid
        grid[i] = (i-1)*increment + grid_min
    end
    # avoid rounding error
    if grid[num_grid] != grid_max
        grid[num_grid] = grid_max
    end
    return grid
end