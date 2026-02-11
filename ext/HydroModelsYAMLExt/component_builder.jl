"""
YAML component builder - constructs HydroModels components from YAML dictionaries.

This module provides functions to build HydroFlux, StateFlux, HydroBucket,
and other components from parsed YAML configuration.
"""

"""
    build_flux_from_yaml(flux_dict::Dict, global_vars::Dict, global_params::Dict)

Build a HydroFlux component from YAML dictionary.

# Arguments
- `flux_dict`: Dictionary containing flux configuration
- `global_vars`: Global variables dictionary
- `global_params`: Global parameters dictionary

# Returns
- HydroFlux component
"""
function build_flux_from_yaml(flux_dict::Dict, global_vars::Dict{Symbol,Num}, global_params::Dict{Symbol,Num})
    # Extract formula
    if !haskey(flux_dict, "formula")
        error("Flux definition missing 'formula' field")
    end

    formula_str = flux_dict["formula"]

    # Parse formula
    equation = parse_formula(formula_str, global_vars, global_params)

    # Create HydroFlux using symbolic constructor
    return HydroModels.HydroFlux(exprs=[equation])
end

"""
    build_stateflux_from_yaml(stateflux_dict::Dict, global_vars::Dict, global_params::Dict)

Build a StateFlux component from YAML dictionary.

# Arguments
- `stateflux_dict`: Dictionary containing state flux configuration
- `global_vars`: Global variables dictionary
- `global_params`: Global parameters dictionary

# Returns
- StateFlux component
"""
function build_stateflux_from_yaml(stateflux_dict::Dict, global_vars::Dict{Symbol,Num}, global_params::Dict{Symbol,Num})
    # Extract formula
    if !haskey(stateflux_dict, "formula")
        error("StateFlux definition missing 'formula' field")
    end

    formula_str = stateflux_dict["formula"]

    # Parse formula
    equation = parse_formula(formula_str, global_vars, global_params)

    # Create StateFlux
    return HydroModels.StateFlux(exprs=[equation])
end

"""
    build_bucket_from_yaml(bucket_dict::Dict, global_vars::Dict, global_params::Dict)

Build a HydroBucket component from YAML dictionary.

# Arguments
- `bucket_dict`: Dictionary containing bucket configuration
- `global_vars`: Global variables dictionary
- `global_params`: Global parameters dictionary

# Returns
- HydroBucket component
"""
function build_bucket_from_yaml(bucket_dict::Dict, global_vars::Dict{Symbol,Num}, global_params::Dict{Symbol,Num})
    # Extract name
    name = Symbol(get(bucket_dict, "name", "unnamed_bucket"))

    # Build fluxes
    fluxes = HydroModels.HydroFlux[]
    if haskey(bucket_dict, "fluxes")
        for flux_def in bucket_dict["fluxes"]
            flux = build_flux_from_yaml(flux_def, global_vars, global_params)
            push!(fluxes, flux)
        end
    end

    # Build state fluxes
    dfluxes = HydroModels.StateFlux[]
    if haskey(bucket_dict, "state_fluxes")
        for dflux_def in bucket_dict["state_fluxes"]
            dflux = build_stateflux_from_yaml(dflux_def, global_vars, global_params)
            push!(dfluxes, dflux)
        end
    end

    # Extract HRU types
    htypes = nothing
    if haskey(bucket_dict, "htypes")
        htypes = convert(Vector{Int}, bucket_dict["htypes"])
    end

    # Create HydroBucket
    return HydroModels.HydroBucket(name=name, fluxes=fluxes, dfluxes=dfluxes, htypes=htypes)
end

"""
    build_config_from_yaml(config_dict::Dict)

Build a HydroConfig from YAML dictionary.

# Arguments
- `config_dict`: Dictionary containing configuration

# Returns
- HydroConfig object
"""
function build_config_from_yaml(config_dict::Dict)
    # Parse solver
    solver_str = get(config_dict, "solver", "MutableSolver")
    solver = if solver_str == "MutableSolver"
        HydroModels.MutableSolver
    elseif solver_str == "ImmutableSolver"
        HydroModels.ImmutableSolver
    elseif solver_str == "ODESolver"
        HydroModels.ODESolver
    elseif solver_str == "DiscreteSolver"
        HydroModels.DiscreteSolver
    else
        error("Unknown solver type: $solver_str")
    end

    # Parse interpolator
    interp_str = get(config_dict, "interpolator", "ConstantInterpolation")
    interpolator = if interp_str in ("ConstantInterpolation", "DirectInterpolation")
        Val(HydroModels.ConstantInterpolation)
    elseif interp_str in ("LinearInterpolation", "EnzymeCompatibleInterpolation", "EnzymeInterpolation")
        Val(HydroModels.LinearInterpolation)
    else
        error("Unknown interpolator type: $interp_str")
    end

    # Other config parameters
    min_value = get(config_dict, "min_value", 1e-6)
    parallel = get(config_dict, "parallel", false)

    return HydroModels.HydroConfig(
        solver=solver,
        interpolator=interpolator,
        min_value=min_value,
        parallel=parallel
    )
end
