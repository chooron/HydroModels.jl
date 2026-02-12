"""
Solver module - provides ODE/difference equation solvers and numerical utilities.
Fully compatible with Zygote AD.
"""

"""
    hydrosolve(solver, du_func, params, initstates, timeidx, config)

Unified ODE/difference equation solver interface.

# Arguments
- `solver::SolverType`: Solver type
- `du_func`: Derivative/update function (u, p, t) -> du
- `params`: Parameter vector (ComponentVector or AbstractVector)
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

Mutable solver - iterative updates without in-place modifications for Zygote compatibility.
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
    state_size = size(initstates)
    results = zeros(T1, state_size..., time_len) |> device

    current_state = T1.(initstates)
    for (i, t) in enumerate(timeidx)
        du = du_func(current_state, params, t)
        current_state = max.(T1(min_val), current_state .+ du)
        results[ntuple(Returns(Colon()), N)..., i] .= current_state
    end

    results
end

"""
    hydrosolve(::Val{ImmutableSolver}, ...)

Immutable solver - functional accumulate, fully compatible with Zygote.
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

    states_vec = accumulate(timeidx; init=initstates) do last_state, t
        du = du_func(last_state, params, t)
        max.(min_val, last_state .+ du)
    end

    stack(states_vec; dims=N + 1) |> device
end

"""
    hydrosolve(::Val{DiscreteSolver}, ...)

Discrete solver - alias for MutableSolver for backward compatibility.
"""
function hydrosolve(
    ::Val{DiscreteSolver},
    du_func,
    params,
    initstates::AbstractArray{T,N},
    timeidx,
    config
) where {T,N}
    # DiscreteSolver is an alias for MutableSolver
    hydrosolve(Val(MutableSolver), du_func, params, initstates, timeidx, config)
end

# ============================================================================
# Numerical utilities
# ============================================================================

"""
    safe_clamp(x, min_val, max_val)

Numerically safe clamping that handles NaN and Inf.
"""
@inline function safe_clamp(x::T, min_val::T, max_val::T) where {T<:Real}
    isnan(x) && return min_val
    isinf(x) && return signbit(x) ? min_val : max_val
    clamp(x, min_val, max_val)
end

"""
    safe_max(x, threshold)

Numerically safe maximum for numerical stability.
"""
@inline function safe_max(x::T, threshold::T) where {T<:Real}
    isnan(x) && return threshold
    max(x, threshold)
end

@inline safe_max(x::AbstractArray, threshold) = @. safe_max(x, threshold)
