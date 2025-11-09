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
           interpolator=Val(DirectInterpolation),
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
        interpolator::Val=Val(DirectInterpolation),
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
        interpolator=Val(DirectInterpolation),
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
normalize_config(config::NamedTuple) = HydroConfig(; config...)
normalize_config(::Nothing) = default_config()

"""
    get_config_value(config::ConfigType, key::Symbol, default)

Get value from configuration with type-safe default fallback.
"""
@inline function get_config_value(config::HydroConfig, key::Symbol, default)
    hasproperty(config, key) ? getproperty(config, key) : default
end

@inline function get_config_value(config::NamedTuple, key::Symbol, default)
    get(config, key, default)
end

# Export main interfaces
export HydroConfig, default_config, merge_config, normalize_config, get_config_value, ConfigType

