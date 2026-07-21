"""
Configuration module - provides type-stable configuration system for compiler optimization.
"""

"""
    HydroConfig{S<:SolverType, I, TI, D, SA, SEA, CB}

Type-stable hydrological model configuration for better performance and type stability.

# Fields
- `solver::S`: Solver type (MutableSolver, ImmutableSolver, ODESolver, DiscreteSolver)
- `interpolator::I`: An interpolator type (constructed from forcing data) or a
  pre-built callable `t -> value`
- `timeidx::TI`: Time index vector, preserving the caller's element type
- `device::D`: Device function (e.g., for GPU acceleration)
- `min_value::Float64`: Minimum value threshold for numerical stability
- `parallel::Bool`: Whether to enable parallel computation
- `solve_alg::SA`: Optional SciML solve algorithm override
- `sense_alg::SEA`: Optional SciML sensitivity algorithm override
- `solve_cb::CB`: Optional SciML callback override

# Examples
```jldoctest
julia> config = HydroConfig(
           solver=MutableSolver,
           interpolator=ConstantInterpolation,
           timeidx=1:100,
           device=identity,
           min_value=1e-6,
           parallel=false
       )
```
"""
struct HydroConfig{S<:SolverType,I,TI,D,SA,SEA,CB}
    solver::S
    interpolator::I
    timeidx::TI
    device::D
    min_value::Float64
    parallel::Bool
    solve_alg::SA
    sense_alg::SEA
    solve_cb::CB
    
    function HydroConfig(;
        solver::SolverType=MutableSolver,
        interpolator=ConstantInterpolation,
        timeidx::Union{Vector{Int},AbstractVector{<:Integer}}=Int[],
        device=identity,
        min_value::Real=1e-6,
        parallel::Bool=false,
        solve_alg=nothing,
        sense_alg=nothing,
        solve_cb=nothing,
        solvealg=nothing,
        sensealg=nothing,
        callback=nothing
    )
        min_value > 0 || throw(ArgumentError("min_value must be positive, got $min_value"))
        interpolator isa Val && throw(ArgumentError(
            "`interpolator=Val(...)` was removed in HydroModels v0.7; " *
            "pass the interpolator type or a pre-built callable directly",
        ))
        isnothing(solve_alg) || isnothing(solvealg) ||
            throw(ArgumentError("Both solve_alg and solvealg were provided; use only one"))
        isnothing(sense_alg) || isnothing(sensealg) ||
            throw(ArgumentError("Both sense_alg and sensealg were provided; use only one"))
        isnothing(solve_cb) || isnothing(callback) ||
            throw(ArgumentError("Both solve_cb and callback were provided; use only one"))

        solve_alg_ = isnothing(solve_alg) ? solvealg : solve_alg
        sense_alg_ = isnothing(sense_alg) ? sensealg : sense_alg
        solve_cb_ = isnothing(solve_cb) ? callback : solve_cb

        timeidx_vec = collect(timeidx)
        new{typeof(solver),typeof(interpolator),typeof(timeidx_vec),typeof(device),typeof(solve_alg_),typeof(sense_alg_),typeof(solve_cb_)}(
            solver,
            interpolator,
            timeidx_vec,
            device,
            Float64(min_value),
            parallel,
            solve_alg_,
            sense_alg_,
            solve_cb_,
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
        interpolator=ConstantInterpolation,
        timeidx=Int[],
        device=identity,
        min_value=1e-6,
        parallel=false,
        solve_alg=nothing,
        sense_alg=nothing,
        solve_cb=nothing,
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
        parallel=config.parallel,
        solve_alg=config.solve_alg,
        sense_alg=config.sense_alg,
        solve_cb=config.solve_cb,
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
const HYDRO_CONFIG_SUPPORTED_KEYS = (
    :solver,
    :interpolator,
    :timeidx,
    :device,
    :min_value,
    :parallel,
    :solve_alg,
    :sense_alg,
    :solve_cb,
)

@inline function _canonical_config_key(key::Symbol)
    key === :solvealg && return :solve_alg
    key === :sensealg && return :sense_alg
    key === :callback && return :solve_cb
    return key
end

@inline function _alias_config_key(key::Symbol)
    key === :solve_alg && return :solvealg
    key === :sense_alg && return :sensealg
    key === :solve_cb && return :callback
    return nothing
end

@inline function _extract_supported_overrides(config::NamedTuple)
    overrides = NamedTuple()
    for key in keys(config)
        canonical = _canonical_config_key(key)
        canonical in HYDRO_CONFIG_SUPPORTED_KEYS || continue
        overrides = merge(overrides, (; canonical => config[key]))
    end
    return overrides
end

"""
    normalize_config(config::ConfigType)

Normalize configuration to HydroConfig type.
"""
normalize_config(config::HydroConfig) = config
function normalize_config(config::NamedTuple)
    overrides = _extract_supported_overrides(config)

    # Nested config from extensions (e.g., OptimizationExt)
    if haskey(config, :hydro_config)
        base_cfg = normalize_config(config.hydro_config)
        return length(overrides) == 0 ? base_cfg : merge_config(base_cfg; overrides...)
    end

    return length(overrides) == 0 ? default_config() : HydroConfig(; overrides...)
end
normalize_config(::Nothing) = default_config()

"""
    get_config_value(config::ConfigType, key::Symbol, default)

Get value from configuration with type-safe default fallback.
Handles nested config structures (e.g., from OptimizationExt).
"""
@inline function get_config_value(config::HydroConfig, key::Symbol, default)
    canonical_key = _canonical_config_key(key)
    hasproperty(config, canonical_key) ? getproperty(config, canonical_key) : default
end

@inline function get_config_value(config::NamedTuple, key::Symbol, default)
    canonical_key = _canonical_config_key(key)

    # Check if key exists directly in config
    if haskey(config, canonical_key)
        return config[canonical_key]
    end

    alias_key = _alias_config_key(canonical_key)
    if !isnothing(alias_key) && haskey(config, alias_key)
        return config[alias_key]
    end

    # Check if this is a nested config with hydro_config field
    if haskey(config, :hydro_config)
        hydro_cfg = config.hydro_config
        if hydro_cfg isa HydroConfig && hasproperty(hydro_cfg, canonical_key)
            return getproperty(hydro_cfg, canonical_key)
        elseif hydro_cfg isa NamedTuple
            if haskey(hydro_cfg, canonical_key)
                return hydro_cfg[canonical_key]
            end
            if !isnothing(alias_key) && haskey(hydro_cfg, alias_key)
                return hydro_cfg[alias_key]
            end
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

