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

# runtime generated functions
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

# Symbolic building
using Symbolics
using Symbolics: tosymbol, unwrap
using SymbolicUtils.Code
import SymbolicUtils: BasicSymbolic, Sym, Term, iscall, operation, arguments, issym, symtype, sorted_arguments
@reexport using ModelingToolkit: @variables, @parameters, Num, isparameter, get_variables
# graph compute
using Graphs

# deep learning
using Lux
using Lux: foldl_init
using NNlib
using MLUtils

@reexport using HydroModelCore:
    AbstractComponent, AbstractFlux, AbstractHydroFlux, AbstractNeuralFlux, AbstractStateFlux,
    AbstractElement, AbstractBucket, AbstractHydrograph, AbstractRoute, AbstractHydroRoute, AbstractModel

@reexport using HydroModelCore:
    get_name, get_input_names, get_output_names, get_param_names,
    get_state_names, get_nn_names, get_var_names,
    get_exprs, get_inputs, get_outputs, get_params, get_nns, get_vars

# utils
include("utils/check.jl")
include("utils/aggregation.jl")
include("utils/expression.jl")
include("utils/tools.jl")
include("utils/miscellaneous.jl")
include("utils/build.jl")

#! A discrete ODE solver, if want to use more efficient solver, please import HydroModelTools.jl
#! When constructing an ODE problem to solve, use DataInterpolations.jl
export ManualSolver, DirectInterpolation

# framework build
include("flux.jl")
export HydroFlux, StateFlux, @hydroflux, @stateflux
include("bucket.jl")
export HydroBucket, @hydrobucket
include("nn.jl")
export NeuralFlux, NeuralBucket, @neuralflux, @neuralbucket
include("route.jl")
export HydroRoute, RapidRoute, @hydroroute
include("uh.jl")
export UnitHydrograph, @unithydro
include("model.jl")
export HydroModel, @hydromodel

## package version
const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

end # module HydroModels
