module HydroModelsDataInterpolationsExt

using HydroModels
using DataInterpolations

"""
    hydrointerp(::Val{DataInterpolations.LinearInterpolation}, input::AbstractMatrix, timeidx)

Create a DataInterpolations.LinearInterpolation-based interpolator.

DataInterpolations supports 2D data directly, returning interpolated vectors.
"""
function HydroModels.hydrointerp(::Val{DataInterpolations.LinearInterpolation}, input::AbstractMatrix, timeidx)
    ts = Float64.(timeidx)
    return DataInterpolations.LinearInterpolation(input, ts)
end

"""
    hydrointerp(::Val{DataInterpolations.CubicSpline}, input::AbstractMatrix, timeidx)

Create a DataInterpolations.CubicSpline-based interpolator.

DataInterpolations supports 2D data directly, returning interpolated vectors.
"""
function HydroModels.hydrointerp(::Val{DataInterpolations.CubicSpline}, input::AbstractMatrix, timeidx)
    ts = Float64.(timeidx)
    return DataInterpolations.CubicSpline(input, ts)
end

end # module
