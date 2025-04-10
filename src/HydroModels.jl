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

abstract type AbstractComponent end

abstract type AbstractFlux <: AbstractComponent end
abstract type AbstractHydroFlux <: AbstractFlux end
abstract type AbstractNeuralFlux <: AbstractHydroFlux end
abstract type AbstractStateFlux <: AbstractFlux end

abstract type AbstractElement <: AbstractComponent end
abstract type AbstractBucket <: AbstractElement end
abstract type AbstractHydrograph <: AbstractElement end
abstract type AbstractRoute <: AbstractElement end
abstract type AbstractHydroRoute <: AbstractRoute end
abstract type AbstractModel <: AbstractComponent end

abstract type AbstractNNLayer <: AbstractComponent end
abstract type AbstractNNModel <: AbstractComponent end

# utils
include("utils/attribute.jl")
include("utils/check.jl")
include("utils/display.jl")
include("utils/aggregation.jl")
include("utils/expression.jl")
include("utils/tools.jl")
include("utils/miscellaneous.jl")
include("utils/build.jl")
#! A discrete ODE solver, if want to use more efficient solver, please import HydroModelSolvers.jl
#! When constructing an ODE problem to solve, use DataInterpolations.jl
export ManualSolver, DirectInterpolation

# framework build
include("flux.jl")
export HydroFlux, StateFlux, NeuralFlux
include("bucket.jl")
export HydroBucket
include("route.jl")
export HydroRoute, RapidRoute
include("uh.jl")
export UHFunction, UnitHydrograph
include("model.jl")
export HydroModel
include("nn.jl")
export HydroNNLayer, HydroNNModel
export @hydroflux, @stateflux, @neuralflux, @hydrobucket, @hydroroute, @hydromodel

## package version
const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

end # module HydroModels
