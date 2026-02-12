"""
Configuration module - provides type-stable configuration system for compiler optimization.
"""

"""
    HydroConfig{S<:SolverType, I, D}

Type-stable hydrological model configuration for better performance and type stability.

# Fields
- `solver::S`: Solver type (MutableSolver, ImmutableSolver, ODESolver, DiscreteSolver)
- `interpolator::I`: Interpolator type, wrapped in Val for type stability
- `timeidx::Vector{Int}`: Time index vector
- `device::D`: Device function (e.g., for GPU acceleration)
- `min_value::Float64`: Minimum value threshold for numerical stability
- `parallel::Bool`: Whether to enable parallel computation

# Examples
```jldoctest
julia> config = HydroConfig(
           solver=MutableSolver,
           interpolator=Val(ConstantInterpolation),
           timeidx=1:100,
           device=identity,
           min_value=1e-6,
           parallel=false
       )
```
"""
struct HydroConfig{S<:SolverType,I,D}
    solver::S
    interpolator::I
    timeidx::Vector{Int}
    device::D
    min_value::Float64
    parallel::Bool
    
    function HydroConfig(;
        solver::SolverType=MutableSolver,
        interpolator::Val=Val(ConstantInterpolation),
        timeidx::Union{Vector{Int},AbstractVector{<:Integer}}=Int[],
        device=identity,
        min_value::Real=1e-6,
        parallel::Bool=false
    )
        min_value > 0 || throw(ArgumentError("min_value must be positive, got $min_value"))
        timeidx_vec = collect(Int, timeidx)
        new{typeof(solver),typeof(interpolator),typeof(device)}(
            solver, interpolator, timeidx_vec, device, Float64(min_value), parallel
        )
    end
end

"""
    default_config()

Create a default configuration instance.

# Returns
- `HydroConfig`: Default configuration object
"""
function default_config()
    HydroConfig(
        solver=MutableSolver,
        interpolator=Val(ConstantInterpolation),
        timeidx=Int[],
        device=identity,
        min_value=1e-6,
        parallel=false
    )
end

"""
    to_namedtuple(config::HydroConfig)

Convert HydroConfig to NamedTuple for backward compatibility.
"""
function to_namedtuple(config::HydroConfig)
    (
        solver=config.solver,
        interpolator=config.interpolator,
        timeidx=config.timeidx,
        device=config.device,
        min_value=config.min_value,
        parallel=config.parallel
    )
end

"""
    merge_config(base::HydroConfig; kwargs...)

Create a new configuration based on existing config, modifying only specified fields.
"""
function merge_config(base::HydroConfig; kwargs...)
    nt = to_namedtuple(base)
    new_nt = merge(nt, kwargs)
    HydroConfig(; new_nt...)
end

# Backward compatibility: allow NamedTuple as configuration
const ConfigType = Union{HydroConfig,NamedTuple}

"""
    normalize_config(config::ConfigType)

Normalize configuration to HydroConfig type.
"""
normalize_config(config::HydroConfig) = config
function normalize_config(config::NamedTuple)
    # Check if this is a nested config with hydro_config field (from OptimizationExt)
    if haskey(config, :hydro_config)
        # Extract the actual HydroConfig, ignoring optimization-specific fields
        return normalize_config(config.hydro_config)
    end
    # Filter out unsupported keys and only keep HydroConfig-compatible fields
    supported_keys = (:solver, :interpolator, :timeidx, :device, :min_value, :parallel)
    filtered_config = NamedTuple{Tuple(k for k in keys(config) if k in supported_keys)}(
        Tuple(config[k] for k in keys(config) if k in supported_keys)
    )
    HydroConfig(; filtered_config...)
end
normalize_config(::Nothing) = default_config()

"""
    get_config_value(config::ConfigType, key::Symbol, default)

Get value from configuration with type-safe default fallback.
Handles nested config structures (e.g., from OptimizationExt).
"""
@inline function get_config_value(config::HydroConfig, key::Symbol, default)
    hasproperty(config, key) ? getproperty(config, key) : default
end

@inline function get_config_value(config::NamedTuple, key::Symbol, default)
    # Check if key exists directly in config
    if haskey(config, key)
        return config[key]
    end
    # Check if this is a nested config with hydro_config field
    if haskey(config, :hydro_config)
        hydro_cfg = config.hydro_config
        if hydro_cfg isa HydroConfig && hasproperty(hydro_cfg, key)
            return getproperty(hydro_cfg, key)
        elseif hydro_cfg isa NamedTuple && haskey(hydro_cfg, key)
            return hydro_cfg[key]
        end
    end
    return default
end

"""
    extract_hydro_config(config::ConfigType)

Extract the core HydroConfig from any config structure.
This is useful when working with nested configs from optimization.
"""
extract_hydro_config(config::HydroConfig) = config
function extract_hydro_config(config::NamedTuple)
    if haskey(config, :hydro_config)
        return extract_hydro_config(config.hydro_config)
    else
        return normalize_config(config)
    end
end

# Export main interfaces
export HydroConfig, default_config, merge_config, normalize_config, get_config_value, ConfigType
export extract_hydro_config

