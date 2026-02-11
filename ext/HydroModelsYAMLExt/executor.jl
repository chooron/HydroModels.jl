"""
YAML Executor for HydroModels

This module provides functionality to execute hydrological models directly from YAML configuration,
including data loading from CSV and NetCDF files.

# Main Functions
- `execute_from_yaml(yaml_file::String)`: Execute model from YAML configuration
- `load_data_from_config(data_config::Dict)`: Load data from configuration

# Example YAML with Data Configuration
```yaml
version: "1.0"
schema: "hydromodels"

parameters:
  Smax:
    description: "Maximum soil moisture"
    units: "mm"
    default: 100.0

components:
  - type: HydroBucket
    name: soil
    fluxes:
      - formula: "baseflow ~ Qmax * exp(-f * (Smax - soilwater))"

model:
  type: HydroModel
  name: simple_model
  components: [soil]

config:
  solver: MutableSolver
  interpolator: DirectInterpolation

data:
  type: csv
  path: "forcing.csv"
  variables:
    temp: Temp        # model_var: file_var
    prcp: Precipitation
    pet: PET
  time_column: date

parameters_data:
  type: csv
  path: "params.csv"

initial_states:
  soilwater: 50.0
```
"""

using CSV
using DataFrames

"""
    load_csv_data(file_path::AbstractString, var_mapping::Dict, time_col::Union{String,Nothing}=nothing)

Load data from CSV file with variable name mapping.

# Arguments
- `file_path`: Path to CSV file
- `var_mapping`: Dictionary mapping model variables to file column names
- `time_col`: Name of time column (optional)

# Returns
- Dictionary with model variable names as keys and data arrays as values
"""
function load_csv_data(file_path::AbstractString, var_mapping::Dict, time_col::Union{String,Nothing}=nothing)
    # Read CSV file
    df = CSV.read(file_path, DataFrame)

    # Extract data with variable mapping
    data = Dict{Symbol,Vector{Float64}}()

    for (model_var, file_var) in var_mapping
        model_var_sym = Symbol(model_var)
        file_var_str = String(file_var)

        if !hasproperty(df, file_var_str)
            error("Column '$file_var_str' not found in CSV file. Available columns: $(names(df))")
        end

        data[model_var_sym] = Float64.(df[!, file_var_str])
    end

    # Extract time if specified
    if !isnothing(time_col) && hasproperty(df, time_col)
        data[:time] = df[!, time_col]
    end

    return data
end

"""
    load_data_from_config(data_config::Dict)

Load forcing data from configuration dictionary.

# Arguments
- `data_config`: Dictionary containing data configuration

# Returns
- Dictionary with variable names as keys and data arrays as values
"""
function load_data_from_config(data_config::Dict)
    data_type = get(data_config, "type", "csv")

    if data_type == "csv"
        file_path = data_config["path"]
        var_mapping = get(data_config, "variables", Dict())
        time_col = get(data_config, "time_column", nothing)

        return load_csv_data(file_path, var_mapping, time_col)
    elseif data_type == "netcdf" || data_type == "nc"
        # NetCDF support requires Rasters extension
        if !isdefined(Base, :get_extension)
            error("NetCDF support requires Julia 1.9+ with package extensions")
        end

        file_path = data_config["path"]
        var_mapping = get(data_config, "variables", Dict())

        # Check if spatial positions are provided
        if haskey(data_config, "positions")
            positions = [(p["row"], p["col"]) for p in data_config["positions"]]
            return load_netcdf_data(file_path, var_mapping, positions)
        else
            # Lumped model - simple NetCDF loading
            return load_netcdf_data_simple(file_path, var_mapping)
        end
    else
        error("Unsupported data type: $data_type. Supported types: 'csv', 'netcdf', 'nc'")
    end
end

"""
    load_parameters_data(params_config::Dict)

Load parameter values from configuration.

# Arguments
- `params_config`: Dictionary containing parameter data configuration

# Returns
- Dictionary with parameter names as keys and values
"""
function load_parameters_data(params_config::Dict)
    if haskey(params_config, "type") && params_config["type"] == "csv"
        file_path = params_config["path"]
        df = CSV.read(file_path, DataFrame)

        # Convert first row to parameter dictionary
        params = Dict{Symbol,Float64}()
        for col in names(df)
            params[Symbol(col)] = Float64(df[1, col])
        end

        return params
    elseif haskey(params_config, "values")
        # Direct parameter values
        params = Dict{Symbol,Float64}()
        for (k, v) in params_config["values"]
            params[Symbol(k)] = Float64(v)
        end
        return params
    else
        error("Invalid parameter configuration")
    end
end

"""
    prepare_forcing_matrix(data::Dict{Symbol,Vector{Float64}}, var_order::Vector{Symbol})

Prepare forcing data matrix with specified variable order.

# Arguments
- `data`: Dictionary of variable data
- `var_order`: Ordered list of variable names

# Returns
- Matrix with shape (n_timesteps, n_variables)
"""
function prepare_forcing_matrix(data::Dict{Symbol,Vector{Float64}}, var_order::Vector{Symbol})
    # Get number of timesteps from first variable
    n_steps = length(first(values(data)))
    n_vars = length(var_order)

    # Create matrix
    forcing = zeros(Float64, n_steps, n_vars)

    for (i, var) in enumerate(var_order)
        if !haskey(data, var)
            error("Variable $var not found in data. Available: $(keys(data))")
        end
        forcing[:, i] = data[var]
    end

    return forcing
end

"""
    HydroModels.execute_from_yaml(yaml_file::AbstractString; return_components::Bool=false)

Execute a hydrological model from YAML configuration file.

# Arguments
- `yaml_file`: Path to YAML configuration file
- `return_components`: If true, return (output, model, config, data) tuple

# Returns
- Model output, or tuple of (output, model, config, data) if return_components=true

# Example
```julia
using HydroModels
using YAML

# Execute model
output = execute_from_yaml("model_config.yaml")

# Or get all components
output, model, config, data = execute_from_yaml("model_config.yaml", return_components=true)
```
"""
function HydroModels.execute_from_yaml(yaml_file::AbstractString; return_components::Bool=false)
    # Load YAML file
    yaml_dict = YAML.load_file(yaml_file)

    # Load model
    model = HydroModels.load_model_from_yaml(yaml_file)

    # Load configuration
    config = HydroModels.load_config_from_yaml(yaml_file)

    # Load forcing data
    if !haskey(yaml_dict, "data")
        error("YAML file missing 'data' section for forcing data")
    end
    forcing_data = load_data_from_config(yaml_dict["data"])

    # Load parameters
    params = if haskey(yaml_dict, "parameters_data")
        load_parameters_data(yaml_dict["parameters_data"])
    elseif haskey(yaml_dict, "parameters")
        # Use default values from parameter definitions
        params_dict = Dict{Symbol,Float64}()
        for (param_name, param_def) in yaml_dict["parameters"]
            if haskey(param_def, "default")
                params_dict[Symbol(param_name)] = Float64(param_def["default"])
            end
        end
        params_dict
    else
        error("No parameter values provided")
    end

    # Get variable order from data configuration
    var_order = if haskey(yaml_dict["data"], "variables")
        [Symbol(k) for k in keys(yaml_dict["data"]["variables"])]
    else
        # Use all non-time variables
        [k for k in keys(forcing_data) if k != :time]
    end

    # Prepare forcing matrix
    forcing = prepare_forcing_matrix(forcing_data, var_order)

    # Get initial states
    initial_states = if haskey(yaml_dict, "initial_states")
        ComponentVector(; (Symbol(k) => Float64(v) for (k, v) in yaml_dict["initial_states"])...)
    else
        nothing
    end

    # Convert parameters to ComponentVector
    params_cv = ComponentVector(; params...)

    # Execute model
    output = if isnothing(initial_states)
        model(forcing, params_cv, config)
    else
        model(forcing, params_cv, config, initial_states)
    end

    if return_components
        return output, model, config, forcing_data
    else
        return output
    end
end
