module HydroModels

## External packages
# common packages
using Reexport

@reexport using ComponentArrays
using LinearAlgebra
using SparseArrays
using StableRNGs
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
using ModelingToolkit: t_nounits as t
# graph compute
using Graphs

# data interpolataion
using DataInterpolations
using DataInterpolations: AbstractInterpolation

# deep learning
using Lux
using NNlib

## Abstract Component Types
abstract type AbstractComponent end
abstract type AbstractIOAdapter end
abstract type AbstractHydroWrapper <: AbstractComponent end
abstract type AbstractNeuralWrapper <: AbstractComponent end

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

export AbstractComponent, AbstractHydroWrapper, AbstractNeuralWrapper
export AbstractFlux, AbstractHydroFlux, AbstractNeuralFlux, AbstractStateFlux
export AbstractElement, AbstractBucket, AbstractHydrograph, AbstractRoute, AbstractHydroRoute, AbstractModel

# utils
include("utils/expression.jl")
include("utils/attribute.jl")
include("utils/display.jl")
include("utils/build.jl")
include("utils/sort.jl")
include("utils/check.jl")
#! A discrete ODE solver, if want to use more efficient solver, please import HydroModelTools.jl
include("utils/solver.jl")
export ManualSolver

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
include("utils/macros.jl")
export @hydroflux, @stateflux, @neuralflux, @hydrobucket
end # module HydroModels
