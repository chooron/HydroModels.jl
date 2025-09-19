module HydroModels

## External packages
# common packages
using Random
using Reexport
using SpecialFunctions
using ComponentArrays
using DocStringExtensions
# runtime generated functions
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

# Symbolic building
using Symbolics
using Symbolics: tosymbol, Num, @variables, get_variables
using SymbolicUtils.Code
import SymbolicUtils: symtype, term, hasmetadata, issym

# graph compute
using Graphs

# deep learning
using Lux
using NNlib

# define Optional type
const Optional{T} = Union{T,Nothing}

using HydroModelCore

# utils
export @variables, @parameters

include("utils.jl")
export sort_components, sort_fluxes
#! When constructing an ODE problem to solve, use DataInterpolations.jl
include("tools.jl")
export ManualSolver, DirectInterpolation
# default hydro config
DEFAULT_CONFIG = (
    solver=ManualSolver(mutable=true),
    interpolator=DirectInterpolation,
    timeidx=Int[],
    device=identity
)
# framework build
include("flux.jl")
export HydroFlux, StateFlux, @hydroflux, @stateflux, NeuralFlux, @neuralflux
include("bucket.jl")
export HydroBucket, @hydrobucket
include("route.jl")
export HydroRoute, @hydroroute
include("uh.jl")
export UnitHydrograph, @unithydro
include("model.jl")
export HydroModel, @hydromodel

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5
export step_func

smoothlogistic_func(S, Smax, r=0.01, e=5.0) = 1 / (1 + exp((S - r * e * Smax) / (r * Smax)))
export smoothlogistic_func


end # module HydroModels
