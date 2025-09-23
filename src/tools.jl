"""
$(SIGNATURES)

A lightweight interpolation type that provides the same interface as DataInterpolations.jl but uses direct indexing instead of interpolation algorithms.
"""
@inline hydrointerp(::Val{I}, input, timeidx) where {I} = I(input, timeidx)

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

    function DirectInterpolation(data::AbstractArray, ts::AbstractVector{<:Integer}; kwargs...)
        @assert size(data)[end] == length(ts) "The last dimension of data must match the length of ts"
        return new{length(size(data))}(data, ts)
    end
end

@inline (interpolater::DirectInterpolation{1})(t::Integer) = interpolater.data[t]
@inline (interpolater::DirectInterpolation{1})(t::Number) = interpolater.data[ceil(Int, t)]
@inline (interpolater::DirectInterpolation{2})(t::Integer) = interpolater.data[:, t]
@inline (interpolater::DirectInterpolation{2})(t::Number) = interpolater.data[:, ceil(Int, t)]

"""
$(SIGNATURES)

Solve the ODEs for the states.
"""
hydrosolve(solver::SolverType, du_func, pas, initstates, timeidx, config) = hydrosolve(Val(solver), du_func, pas, initstates, timeidx, config)

function hydrosolve(::Val{MutableSolver}, du_func, pas, initstates::AbstractArray{T,N}, timeidx, config) where {T,N}
    dev = get(config, :device, identity)
    T1 = promote_type(eltype(pas), eltype(initstates))
    states_results = zeros(T1, size(initstates)..., length(timeidx)) |> dev
    tmp_initstates = copy(initstates)
    for (i, t) in enumerate(timeidx)
        tmp_du = du_func(tmp_initstates, pas, t)
        # todo custom the near zero value
        tmp_initstates = max.(1e-6, tmp_initstates .+ tmp_du)
        states_results[ntuple(Returns(Colon()), N)..., i] .= tmp_initstates
    end
    states_results
end

function hydrosolve(::Val{ImmutableSolver}, du_func, pas, initstates, timeidx, config)
    dev = get(config, :device, identity)
    states_vec = accumulate(timeidx; init=initstates) do last_state, t
        max.(1e-6, du_func(last_state, pas, t) .+ last_state)
    end
    return stack(states_vec; dims=N + 1) |> dev
end