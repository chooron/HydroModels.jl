module HydroModels

## External packages
# common packages
using Reexport

@reexport using ComponentArrays
using LinearAlgebra
using SparseArrays
using StableRNGs
using Random
using TOML

## package version
const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

# runtime generated functions
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

# Symbolic building
using Symbolics
using Symbolics: tosymbol, unwrap
using SymbolicUtils.Code
import SymbolicUtils: BasicSymbolic, Sym, Term, iscall, operation, arguments, issym, symtype, sorted_arguments
@reexport using ModelingToolkit: @variables, @parameters
using ModelingToolkit: isparameter
# graph compute
using Graphs

# deep learning
using Lux
using Lux: foldl_init
using NNlib
using MLUtils

@reexport using HydroModelCore:
    AbstractComponent,
    AbstractFlux,
    AbstractHydroFlux,
    AbstractNeuralFlux,
    AbstractStateFlux,
    AbstractElement,
    AbstractBucket,
    AbstractHydrograph,
    AbstractRoute,
    AbstractHydroRoute,
    AbstractModel,
    AbstractNNLayer,
    AbstractNNModel

using HydroModelCore: get_var_names, get_state_names, get_output_names, get_input_names, get_param_names, get_nn_names, get_name

# utils
include("utils/expression.jl")
include("utils/miscellaneous.jl")
include("utils/display.jl")
include("utils/build.jl")
include("utils/tools.jl")
include("utils/macros.jl")
export @hydroflux, @stateflux, @neuralflux, @hydrobucket, @hydroroute, @hydromodel
#! A discrete ODE solver, if want to use more efficient solver, please import HydroModelSolvers.jl
#! When constructing an ODE problem to solve, use DataInterpolations.jl
export ManualSolver, DirectInterpolation

# framework build
include("flux.jl")
export HydroFlux, StateFlux, NeuralFlux
include("bucket.jl")
export HydroBucket
include("route.jl")
export GridRoute, VectorRoute, HydroRoute, RapidRoute
include("uh.jl")
export UHFunction, UnitHydrograph
include("model.jl")
export HydroModel
include("nn.jl")
export HydroNNLayer, HydroNNModel

end # module HydroModels
