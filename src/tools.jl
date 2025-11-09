"""
Tools module - provides interpolation and solver functionality, fully compatible with Zygote AD.
"""

"""
    DirectInterpolation{N,T,V}

Lightweight direct indexing interpolator for optimal performance, avoiding complex interpolation calculations.

# Type Parameters
- `N`: Data dimensionality
- `T`: Data element type
- `V`: Time index type

# Fields
- `data::AbstractArray{T,N}`: Data array
- `ts::V`: Time index vector

# Interface
Implements the same callable interface as DataInterpolations.jl, but uses direct indexing instead of interpolation algorithms.

# Examples
```jldoctest
julia> data = rand(10, 100);
julia> ts = 1:100;
julia> interp = DirectInterpolation(data, ts);
julia> value = interp(1.7);  # Returns data[:, 2] (ceiling of 1.7)
```
"""
struct DirectInterpolation{N,T,V<:AbstractVector{<:Integer}}
    data::AbstractArray{T,N}
    ts::V
    
    function DirectInterpolation(data::AbstractArray{T,N}, ts::AbstractVector{<:Integer}) where {T,N}
        @assert size(data, N) == length(ts) "Last dimension of data must match length of time index"
        new{N,T,typeof(ts)}(data, ts)
    end
end

# 1D interpolation
@inline (interp::DirectInterpolation{1})(t::Integer) = interp.data[t]
@inline (interp::DirectInterpolation{1})(t::Real) = interp.data[ceil(Int, t)]

# 2D interpolation - returns column vector
@inline (interp::DirectInterpolation{2})(t::Integer) = @view interp.data[:, t]
@inline (interp::DirectInterpolation{2})(t::Real) = @view interp.data[:, ceil(Int, t)]

"""
    hydrointerp(::Val{I}, input, timeidx) where {I}

Factory function for creating interpolators with type-stable dispatch.

# Arguments
- `::Val{I}`: Interpolator type tag
- `input`: Input data
- `timeidx`: Time index

# Returns
- Interpolator instance

# Note
Uses Val for type-stable dispatch, ensuring compiler optimization.
"""
@inline hydrointerp(::Val{DirectInterpolation}, input, timeidx) = DirectInterpolation(input, timeidx)
@inline hydrointerp(::Val{I}, input, timeidx) where {I} = I(input, timeidx)

"""
    hydrosolve(solver, du_func, params, initstates, timeidx, config)

Unified ODE/difference equation solver interface.

# Arguments
- `solver::SolverType`: Solver type
- `du_func`: Derivative/update function (u, p, t) -> du
- `params`: Parameter vector or ComponentVector
- `initstates`: Initial states
- `timeidx`: Time index
- `config`: Configuration object

# Returns
- State evolution array with dimensions (state_dims..., time_length)
"""
@inline function hydrosolve(solver::SolverType, du_func, params, initstates, timeidx, config)
    hydrosolve(Val(solver), du_func, params, initstates, timeidx, config)
end

"""
    hydrosolve(::Val{MutableSolver}, ...)

Mutable solver implementation - uses iterative updates with high memory efficiency.
Note: This implementation completely avoids in-place modifications to ensure Zygote compatibility.
"""
function hydrosolve(
    ::Val{MutableSolver},
    du_func,
    params,
    initstates::AbstractArray{T,N},
    timeidx,
    config
) where {T,N}
    device = get_config_value(config, :device, identity)
    min_val = get_config_value(config, :min_value, 1e-6)
    
    T1 = promote_type(eltype(params), T)
    time_len = length(timeidx)
    
    # Preallocate result array
    state_size = size(initstates)
    results = zeros(T1, state_size..., time_len) |> device
    
    # Use functional programming to avoid in-place modifications
    current_state = T1.(initstates)
    
    for (i, t) in enumerate(timeidx)
        du = du_func(current_state, params, t)
        # Use pure functional update to avoid in-place modification
        current_state = max.(T1(min_val), current_state .+ du)
        results[ntuple(Returns(Colon()), N)..., i] .= current_state
    end
    
    results
end

"""
    hydrosolve(::Val{ImmutableSolver}, ...)

Immutable solver implementation - uses functional accumulate, fully compatible with Zygote.
This implementation is better suited for automatic differentiation but may use more memory.
"""
function hydrosolve(
    ::Val{ImmutableSolver},
    du_func,
    params,
    initstates::AbstractArray{T,N},
    timeidx,
    config
) where {T,N}
    device = get_config_value(config, :device, identity)
    min_val = get_config_value(config, :min_value, 1e-6)
    
    # Use accumulate for pure functional solving
    states_vec = accumulate(timeidx; init=initstates) do last_state, t
        du = du_func(last_state, params, t)
        max.(min_val, last_state .+ du)
    end
    
    # Stack vector into array
    stack(states_vec; dims=N + 1) |> device
end

"""
    hydrosolve(::Val{DiscreteSolver}, ...)

Discrete solver implementation - for pure algebraic equations (no state evolution).
"""
function hydrosolve(
    ::Val{DiscreteSolver},
    compute_func,
    params,
    input,
    timeidx,
    config
) where {T,N}
    device = get_config_value(config, :device, identity)
    
    # Direct computation for each time step
    results = map(timeidx) do t
        compute_func(input, params, t)
    end
    
    stack(results; dims=ndims(first(results)) + 1) |> device
end

"""
    safe_clamp(x, min_val, max_val)

Numerically safe clamping function that handles NaN and Inf values.
"""
@inline function safe_clamp(x::T, min_val::T, max_val::T) where {T<:Real}
    isnan(x) && return min_val
    isinf(x) && return signbit(x) ? min_val : max_val
    clamp(x, min_val, max_val)
end

"""
    safe_max(x, threshold)

Numerically safe maximum function to ensure numerical stability.
"""
@inline function safe_max(x::T, threshold::T) where {T<:Real}
    isnan(x) && return threshold
    max(x, threshold)
end

# Vectorized version - using @. for better performance
@inline safe_max(x::AbstractArray, threshold) = @. safe_max(x, threshold)

export DirectInterpolation, hydrointerp, hydrosolve, safe_clamp, safe_max
