"""
YAML Loader for HydroModels

This module provides functionality to load hydrological models from YAML configuration files.

# Main Functions
- `load_model_from_yaml(yaml_file::String)`: Load a complete model from YAML file
- `load_config_from_yaml(yaml_file::String)`: Load configuration from YAML file
- `load_parameters_from_yaml(yaml_file::String)`: Load parameter metadata from YAML file

# Example
```julia
using HydroModels
using YAML

# Load model from YAML
model = load_model_from_yaml("exphydro.yaml")

# Run model
output = model(forcing, params, config)
```
"""

"""
    HydroModels.load_model_from_yaml(yaml_file::AbstractString)

Load a complete hydrological model from a YAML configuration file.

# Arguments
- `yaml_file`: Path to YAML configuration file

# Returns
- HydroModel object

# YAML Format
The YAML file should contain:
- `parameters`: Parameter definitions with metadata
- `components`: List of component definitions (HydroBucket, etc.)
- `model`: Model composition specification
- `config`: Optional configuration settings

# Example YAML
```yaml
version: "1.0"
schema: "hydromodels"

parameters:
  Smax:
    description: "Maximum soil moisture"
    units: "mm"
    default: 100.0
    bounds: [50.0, 500.0]

components:
  - type: HydroBucket
    name: soil
    fluxes:
      - name: baseflow
        formula: "baseflow ~ Qmax * exp(-f * (max(0.0, Smax - soilwater)))"
    state_fluxes:
      - name: soil_balance
        formula: "soilwater ~ rainfall - evap - baseflow"

model:
  type: HydroModel
  name: simple_model
  components: [soil]

config:
  solver: MutableSolver
  interpolator: DirectInterpolation
  min_value: 1.0e-6
```
"""
function HydroModels.load_model_from_yaml(yaml_file::AbstractString)
    # Load YAML file
    yaml_dict = YAML.load_file(yaml_file)

    # Extract global parameters
    global_params = Dict{Symbol,Num}()
    if haskey(yaml_dict, "parameters")
        for (param_name, param_def) in yaml_dict["parameters"]
            # Create parameter symbol
            param_sym = Symbol(param_name)
            # Create a Symbolics parameter (so it is not treated as an input)
            param_var = HydroModels.toparam(Symbolics.variable(param_sym))
            global_params[param_sym] = param_var
        end
    end

    # Extract all variables from formulas
    global_vars = Dict{Symbol,Num}()
    if haskey(yaml_dict, "components")
        for comp_def in yaml_dict["components"]
            # Extract variables from fluxes
            if haskey(comp_def, "fluxes")
                for flux_def in comp_def["fluxes"]
                    if haskey(flux_def, "formula")
                        vars = extract_variables_from_formula(flux_def["formula"])
                        for var in vars
                            if !haskey(global_vars, var) && !haskey(global_params, var)
                                # Create variable dynamically using Symbolics functional API
                                var_sym = Symbolics.variable(var)
                                global_vars[var] = var_sym
                            end
                        end
                    end
                end
            end

            # Extract variables from state fluxes
            if haskey(comp_def, "state_fluxes")
                for dflux_def in comp_def["state_fluxes"]
                    if haskey(dflux_def, "formula")
                        vars = extract_variables_from_formula(dflux_def["formula"])
                        for var in vars
                            if !haskey(global_vars, var) && !haskey(global_params, var)
                                # Create variable dynamically using Symbolics functional API
                                var_sym = Symbolics.variable(var)
                                global_vars[var] = var_sym
                            end
                        end
                    end
                end
            end
        end
    end

    # Build components
    components = []
    component_names = Symbol[]
    if haskey(yaml_dict, "components")
        for comp_def in yaml_dict["components"]
            comp_type = get(comp_def, "type", "")

            if comp_type == "HydroBucket"
                comp = build_bucket_from_yaml(comp_def, global_vars, global_params)
                push!(components, comp)
                push!(component_names, comp.name)
            else
                @warn "Unknown component type: $comp_type"
            end
        end
    end

    # Build model
    if haskey(yaml_dict, "model")
        model_def = yaml_dict["model"]
        model_name = Symbol(get(model_def, "name", "yaml_model"))

        # Create HydroModel
        model = HydroModels.HydroModel(name=model_name, components=Tuple(components))
        return model
    else
        error("YAML file missing 'model' section")
    end
end

"""
    HydroModels.load_config_from_yaml(yaml_file::AbstractString)

Load configuration from YAML file.

# Arguments
- `yaml_file`: Path to YAML configuration file

# Returns
- HydroConfig object
"""
function HydroModels.load_config_from_yaml(yaml_file::AbstractString)
    yaml_dict = YAML.load_file(yaml_file)

    if haskey(yaml_dict, "config")
        return build_config_from_yaml(yaml_dict["config"])
    else
        return HydroModels.default_config()
    end
end

"""
    HydroModels.load_parameters_from_yaml(yaml_file::AbstractString)

Load parameter metadata from YAML file.

# Arguments
- `yaml_file`: Path to YAML configuration file

# Returns
- Dictionary of parameter metadata
"""
function HydroModels.load_parameters_from_yaml(yaml_file::AbstractString)
    yaml_dict = YAML.load_file(yaml_file)

    params_metadata = Dict{Symbol,Dict{Symbol,Any}}()

    if haskey(yaml_dict, "parameters")
        for (param_name, param_def) in yaml_dict["parameters"]
            param_sym = Symbol(param_name)
            metadata = Dict{Symbol,Any}()

            if haskey(param_def, "description")
                metadata[:description] = param_def["description"]
            end
            if haskey(param_def, "units")
                metadata[:units] = param_def["units"]
            end
            if haskey(param_def, "default")
                metadata[:default] = param_def["default"]
            end
            if haskey(param_def, "bounds")
                metadata[:bounds] = Tuple(param_def["bounds"])
            end

            params_metadata[param_sym] = metadata
        end
    end

    return params_metadata
end
