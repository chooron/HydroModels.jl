"""
Interpolation module - Enzyme-compatible interpolation implementations.

Provides `ConstantInterpolation` (step/ceiling lookup) and `LinearInterpolation`
(linear interpolation between adjacent points), both fully compatible with
Enzyme.jl automatic differentiation.
"""

# ============================================================================
# ConstantInterpolation
# ============================================================================

"""
    ConstantInterpolation{N,T,V}

Lightweight constant (step) interpolator using ceiling-based indexing.

For non-integer time `t`, returns the value at `ceil(Int, t)`.
Enzyme-compatible, no external dependencies.

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
Enzyme-compatible, no external dependencies.

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
    hydrointerp(::Val{I}, input, timeidx) where {I}

Factory function for creating interpolators with type-stable dispatch via Val.

# Examples
```julia
interp = hydrointerp(Val(ConstantInterpolation), data, timeidx)
interp = hydrointerp(Val(LinearInterpolation), data, timeidx)
```
"""
@inline hydrointerp(::Val{ConstantInterpolation}, input, timeidx) = ConstantInterpolation(input, timeidx)
@inline hydrointerp(::Val{LinearInterpolation}, input, timeidx) = LinearInterpolation(input, timeidx)
@inline hydrointerp(::Val{I}, input, timeidx) where {I} = I(input, timeidx)

# ============================================================================
# Backward compatibility aliases
# ============================================================================

const DirectInterpolation = ConstantInterpolation
const EnzymeInterpolation = ConstantInterpolation
const EnzymeCompatibleInterpolation = ConstantInterpolation
