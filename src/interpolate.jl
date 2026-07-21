"""
Interpolation module - lightweight interpolation implementations.

Provides `ConstantInterpolation` (step/ceiling lookup) and `LinearInterpolation`
(linear interpolation between adjacent points).
"""

# ============================================================================
# ConstantInterpolation
# ============================================================================

"""
    ConstantInterpolation{N,T,V}

Lightweight constant (step) interpolator using ceiling-based indexing.

For non-integer time `t`, returns the value at `ceil(Int, t)`.
No external dependencies.

# Type Parameters
- `N`: Data dimensionality
- `T`: Data element type
- `V`: Time index type

# Examples
```julia
data = rand(10, 100)
ts = 1:100
interp = ConstantInterpolation(data, ts)
value = interp(1.7)  # Returns data[:, 2]
```
"""
struct ConstantInterpolation{N,T,V<:AbstractVector{<:Integer}}
    data::AbstractArray{T,N}
    ts::V

    function ConstantInterpolation(data::AbstractArray{T,N}, ts::AbstractVector{<:Integer}) where {T,N}
        @assert size(data, N) == length(ts) "Last dimension of data must match length of time index"
        new{N,T,typeof(ts)}(data, ts)
    end
end

# 1D: scalar time series
@inline (interp::ConstantInterpolation{1})(t::Integer) = interp.data[t]
@inline (interp::ConstantInterpolation{1})(t::Real) = interp.data[ceil(Int, t)]

# 2D: variables × time
@inline (interp::ConstantInterpolation{2})(t::Integer) = @view interp.data[:, t]
@inline (interp::ConstantInterpolation{2})(t::Real) = @view interp.data[:, ceil(Int, t)]

# ============================================================================
# LinearInterpolation
# ============================================================================

"""
    LinearInterpolation{N,T,V}

Linear interpolation between adjacent time points.
No external dependencies.

# Type Parameters
- `N`: Data dimensionality
- `T`: Data element type
- `V`: Time index type

# Algorithm
- `t <= ts[1]`: return first value (boundary)
- `t >= ts[end]`: return last value (boundary)
- Otherwise: linear interpolation `(1-α) * data[idx-1] + α * data[idx]`

# Examples
```julia
data = rand(10, 100)
ts = 1:100
interp = LinearInterpolation(data, ts)
value = interp(1.5)  # Linearly interpolated between t=1 and t=2
```
"""
struct LinearInterpolation{N,T,V}
    data::AbstractArray{T,N}
    ts::V

    function LinearInterpolation(data::AbstractArray{T,N}, ts::V) where {T,N,V}
        @assert size(data, N) == length(ts) "Last dimension of data must match length of time index"
        new{N,T,V}(data, ts)
    end
end

# 2D: variables × time
@inline function (interp::LinearInterpolation{2})(t::Real)
    idx = searchsortedfirst(interp.ts, t)

    if idx == 1
        return @view interp.data[:, 1]
    elseif idx > length(interp.ts)
        return @view interp.data[:, end]
    end

    t1, t2 = interp.ts[idx-1], interp.ts[idx]
    α = (t - t1) / (t2 - t1)
    return (1 - α) .* @view(interp.data[:, idx-1]) .+ α .* @view(interp.data[:, idx])
end

# 1D: scalar time series
@inline function (interp::LinearInterpolation{1})(t::Real)
    idx = searchsortedfirst(interp.ts, t)

    if idx == 1
        return interp.data[1]
    elseif idx > length(interp.ts)
        return interp.data[end]
    end

    t1, t2 = interp.ts[idx-1], interp.ts[idx]
    α = (t - t1) / (t2 - t1)
    return (1 - α) * interp.data[idx-1] + α * interp.data[idx]
end

# ============================================================================
# Factory function
# ============================================================================

"""
    hydrointerp(interpolator, input, timeidx)

Resolve an interpolation specification.

Pass an interpolator type to construct it from `input` and `timeidx`, or pass
an already-constructed callable when forcing data and its time axis are owned
outside HydroModels.  A pre-built callable must implement `itp(t)`.

# Examples
```julia
interp = hydrointerp(ConstantInterpolation, data, timeidx)
interp = hydrointerp(LinearInterpolation, data, timeidx)

# Preserve an externally prepared interpolator
interp = hydrointerp(LinearInterpolation(data, timeidx), data, timeidx)
```
"""
@inline hydrointerp(interpolator::Type, input, timeidx) = interpolator(input, timeidx)
@inline function hydrointerp(::Val, input, timeidx)
    throw(ArgumentError(
        "`hydrointerp(Val(...), ...)` was removed in HydroModels v0.7; " *
        "pass the interpolator type or a pre-built callable directly",
    ))
end
@inline hydrointerp(interpolator, input, timeidx) = interpolator

# ============================================================================
# Backward compatibility aliases
# ============================================================================

const DirectInterpolation = ConstantInterpolation
