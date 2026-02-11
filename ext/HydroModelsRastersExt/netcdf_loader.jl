"""
NetCDF Data Loader for HydroModels

This module provides functionality to load forcing data from NetCDF files
for hydrological model execution.

# Main Functions
- `load_netcdf_data(file_path, var_mapping, positions)`: Load NetCDF data with variable mapping
- `prepare_netcdf_forcing(data, var_order, positions)`: Prepare forcing matrix from NetCDF data
"""

using Rasters
using NCDatasets

"""
    load_netcdf_data(file_path::String, var_mapping::Dict, positions::Vector{Tuple{Int,Int}})

Load data from NetCDF file with variable name mapping and spatial positions.

# Arguments
- `file_path`: Path to NetCDF file
- `var_mapping`: Dictionary mapping model variables to NetCDF variable names
- `positions`: Vector of (row, col) positions to extract data from

# Returns
- Dictionary with model variable names as keys and data arrays as values
  Each array has shape (n_timesteps, n_positions)
"""
function load_netcdf_data(file_path::String, var_mapping::Dict, positions::Vector{Tuple{Int,Int}})
    # Load NetCDF as Raster stack
    raster_stack = Raster(file_path)

    data = Dict{Symbol,Matrix{Float64}}()

    for (model_var, nc_var) in var_mapping
        model_var_sym = Symbol(model_var)
        nc_var_str = String(nc_var)

        # Load the specific variable
        if haskey(raster_stack, Symbol(nc_var_str))
            var_raster = raster_stack[Symbol(nc_var_str)]
        else
            error("Variable '$nc_var_str' not found in NetCDF file")
        end

        # Extract data at specified positions
        n_times = size(var_raster, 3)  # Assuming time is 3rd dimension
        n_positions = length(positions)

        extracted_data = zeros(Float64, n_times, n_positions)

        for (pos_idx, (row, col)) in enumerate(positions)
            for t in 1:n_times
                extracted_data[t, pos_idx] = var_raster[row, col, t]
            end
        end

        data[model_var_sym] = extracted_data
    end

    return data
end

"""
    load_netcdf_data_simple(file_path::String, var_mapping::Dict)

Load data from NetCDF file for lumped models (single spatial location).

# Arguments
- `file_path`: Path to NetCDF file
- `var_mapping`: Dictionary mapping model variables to NetCDF variable names

# Returns
- Dictionary with model variable names as keys and time series arrays as values
"""
function load_netcdf_data_simple(file_path::String, var_mapping::Dict)
    NCDataset(file_path) do ds
        data = Dict{Symbol,Vector{Float64}}()

        for (model_var, nc_var) in var_mapping
            model_var_sym = Symbol(model_var)
            nc_var_str = String(nc_var)

            if !haskey(ds, nc_var_str)
                error("Variable '$nc_var_str' not found in NetCDF file. Available: $(keys(ds))")
            end

            # Read variable data
            var_data = ds[nc_var_str][:]

            # Handle different dimensions
            if ndims(var_data) == 1
                # Time series data
                data[model_var_sym] = Float64.(var_data)
            elseif ndims(var_data) == 3
                # Spatial data (lat, lon, time) - take mean over space
                data[model_var_sym] = Float64.(vec(mean(var_data, dims=(1,2))))
            else
                error("Unsupported data dimensions for variable '$nc_var_str'")
            end
        end

        return data
    end
end

"""
    prepare_netcdf_forcing(data::Dict{Symbol,Matrix{Float64}}, var_order::Vector{Symbol})

Prepare forcing data matrix from NetCDF data with specified variable order.

# Arguments
- `data`: Dictionary of variable data (n_timesteps, n_positions)
- `var_order`: Ordered list of variable names

# Returns
- 3D array with shape (n_timesteps, n_variables, n_positions)
"""
function prepare_netcdf_forcing(data::Dict{Symbol,Matrix{Float64}}, var_order::Vector{Symbol})
    # Get dimensions
    first_var = data[var_order[1]]
    n_steps, n_positions = size(first_var)
    n_vars = length(var_order)

    # Create 3D array
    forcing = zeros(Float64, n_steps, n_vars, n_positions)

    for (i, var) in enumerate(var_order)
        if !haskey(data, var)
            error("Variable $var not found in data. Available: $(keys(data))")
        end
        forcing[:, i, :] = data[var]
    end

    return forcing
end

"""
    extract_positions_from_netcdf(file_path::String, var_name::String; threshold::Float64=0.0)

Extract valid spatial positions from NetCDF file based on a variable threshold.

# Arguments
- `file_path`: Path to NetCDF file
- `var_name`: Variable name to use for masking
- `threshold`: Minimum value threshold for valid positions

# Returns
- Vector of (row, col) positions
"""
function extract_positions_from_netcdf(file_path::String, var_name::String; threshold::Float64=0.0)
    raster = Raster(file_path, name=Symbol(var_name))

    # Take first time step if 3D
    if ndims(raster) == 3
        raster_2d = raster[:, :, 1]
    else
        raster_2d = raster
    end

    positions = Tuple{Int,Int}[]
    rows, cols = size(raster_2d)

    for i in 1:rows
        for j in 1:cols
            val = raster_2d[i, j]
            if !ismissing(val) && !isnan(val) && val > threshold
                push!(positions, (i, j))
            end
        end
    end

    return positions
end
