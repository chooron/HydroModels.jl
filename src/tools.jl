@kwdef struct HydroConfig{S,I}
    solver::S=ManualSolver(mutable=true)
    interpolater::I=DirectInterpolation
    ptyidx::Vector{Int}=Int[]
    styidx::Vector{Int}=Int[]
end

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
struct ManualSolver{M,D,NZ}
    dev::D
    nearzero::NZ

    function ManualSolver(; dev=identity, nearzero=1e-6, mutable::Bool=false)
        return new{mutable,typeof(dev),typeof(nearzero)}(dev, nearzero)
    end
end

function (solver::ManualSolver{true,D,NZ})(
    du_func::Function,
    pas::AbstractVector,
    initstates::AbstractArray{T,N},
    timeidx::AbstractVector;
    kwargs...
) where {T,D,N,NZ}
    T1 = promote_type(eltype(pas), eltype(initstates))
    states_results = zeros(T1, size(initstates)..., length(timeidx)) |> solver.dev
    tmp_initstates = copy(initstates)
    nearzero = solver.nearzero |> eltype(initstates)
    for (i, t) in enumerate(timeidx)
        tmp_du = du_func(tmp_initstates, pas, t)
        tmp_initstates = max.(nearzero, tmp_initstates .+ tmp_du)
        states_results[ntuple(Returns(Colon()), N)..., i] .= tmp_initstates
    end
    states_results
end

function (solver::ManualSolver{false,D,NZ})(
    du_func::Function,
    pas::AbstractVector,
    initstates::AbstractArray{T,N},
    timeidx::AbstractVector;
    kwargs...
) where {T,D,N,NZ}
    nearzero = solver.nearzero |> eltype(initstates)
    states_vec = accumulate(timeidx; init=initstates) do last_state, t
        max.(nearzero, du_func(last_state, pas, t) .+ last_state)
    end
    return stack(states_vec; dims=N + 1)
end

struct FunctionElement <: AbstractElement
    inputs::Vector{Num}
    func::Function
end

function (ele::FunctionElement)(input::AbstractArray, params::ComponentVector; kwargs...)
    ele.func(input, params; kwargs...)
end

function SummationElement(inputs::Vector{Num})
    return FunctionElement(inputs, (input, params) -> sum(input, dims=2))
end