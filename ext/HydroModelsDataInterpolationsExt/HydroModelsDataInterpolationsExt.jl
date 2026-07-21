module HydroModelsDataInterpolationsExt

using HydroModels
using DataInterpolations

"""
    hydrointerp(::Type{DataInterpolations.LinearInterpolation}, input::AbstractMatrix, timeidx)

Create a DataInterpolations.LinearInterpolation-based interpolator.

DataInterpolations supports 2D data directly, returning interpolated vectors.
"""
function HydroModels.hydrointerp(::Type{DataInterpolations.LinearInterpolation}, input::AbstractMatrix, timeidx)
    ts = Float64.(timeidx)
    return DataInterpolations.LinearInterpolation(input, ts)
end

"""
    hydrointerp(::Type{DataInterpolations.CubicSpline}, input::AbstractMatrix, timeidx)

Create a DataInterpolations.CubicSpline-based interpolator.

DataInterpolations supports 2D data directly, returning interpolated vectors.
"""
function HydroModels.hydrointerp(::Type{DataInterpolations.CubicSpline}, input::AbstractMatrix, timeidx)
    ts = Float64.(timeidx)
    return DataInterpolations.CubicSpline(input, ts)
end

end # module
