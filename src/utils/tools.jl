"""
    DirectInterpolation{D}(data::AbstractArray, ts::AbstractVector{<:Integer})

A lightweight interpolation type that provides the same interface as DataInterpolations.jl but uses direct indexing instead of interpolation algorithms.

# Arguments
- `data::AbstractArray`: The data array to be interpolated
- `ts::AbstractVector{<:Integer}`: The time points corresponding to the data

# Interface
Implements the same callable interface as DataInterpolations.jl interpolation types:
- `(interp::DirectInterpolation)(t)`: Get the value at time `t`

# Implementation Details
Rather than performing interpolation calculations, this type simply returns the data value at the ceiling of the requested time index.
For non-integer time points `t`, it uses `ceil(Int, t)` to get the next integer index.

# Example
```julia
data = rand(100)
ts = 1:100
interp = DirectInterpolation(data, ts)
value = interp(1.7)  # Returns data[2] (ceiling of 1.7)
```
"""
struct DirectInterpolation{D}
    data::AbstractArray
    ts::AbstractVector

    function DirectInterpolation(data::AbstractArray, ts::AbstractVector{<:Integer})
        @assert size(data)[end] == length(ts) "The last dimension of data must match the length of ts"
        return new{length(size(data))}(data, ts)
    end
end

@inline (interpolater::DirectInterpolation{1})(t::Integer) = interpolater.data[t]
@inline (interpolater::DirectInterpolation{1})(t::Number) = interpolater.data[ceil(Int, t)]
@inline (interpolater::DirectInterpolation{2})(t::Integer) = interpolater.data[:, t]
@inline (interpolater::DirectInterpolation{2})(t::Number) = interpolater.data[:, ceil(Int, t)]

"""
    ManualSolver{mutable} <: AbstractHydroSolver

A lightweight, type-stable ODE solver optimized for hydrological modeling.

# Type Parameters
- `mutable::Bool`: Controls array mutability and performance characteristics
    - `true`: Uses mutable arrays for in-place updates (30% faster)
    - `false`: Uses immutable arrays for functional programming style

# Description
ManualSolver implements a fixed-step explicit Euler method for solving ordinary 
differential equations (ODEs). It is specifically designed for hydrological models 
where stability and computational efficiency are prioritized over high-order accuracy.

## Performance Characteristics
- Type-stable implementation for predictable performance
- Optional in-place operations via the `mutable` parameter
- Linear time complexity with respect to simulation length
- Constant space complexity when `mutable=true`

## Memory Management
- `mutable=true`: 
    - Modifies arrays in-place
- `mutable=false`:
    - Creates new arrays at each step

"""
struct ManualSolver{mutable}
    dev

    function ManualSolver(; dev=identity, mutable::Bool=false)
        return new{mutable}(dev)
    end
end

function (solver::ManualSolver{true})(
    du_func::Function,
    pas::AbstractVector,
    initstates::AbstractArray{<:Number,N},
    timeidx::AbstractVector;
    kwargs...
) where N
    T1 = promote_type(eltype(pas), eltype(initstates))
    states_results = zeros(T1, size(initstates)..., length(timeidx)) |> solver.dev
    tmp_initstates = copy(initstates)
    for (i, t) in enumerate(timeidx)
        tmp_du = du_func(tmp_initstates, pas, t)
        tmp_initstates = tmp_initstates .+ tmp_du
        states_results[ntuple(_ -> Colon(), N)..., i] .= tmp_initstates
    end
    states_results
end

function (solver::ManualSolver{false})(
    du_func::Function,
    pas::AbstractVector,
    initstates::AbstractArray{T,N},
    timeidx::AbstractVector;
    kwargs...
) where {T,N}
    function recur_op(::Nothing, t)
        new_states = du_func(initstates, pas, t) .+ initstates
        return [new_states], new_states
    end
    function recur_op((states_list, last_state), t)
        new_states = du_func(last_state, pas, t) .+ last_state
        return vcat(states_list, [new_states]), new_states
    end
    states_vec, _ = foldl_init(recur_op, timeidx)
    stack(states_vec, dims=N + 1)
end