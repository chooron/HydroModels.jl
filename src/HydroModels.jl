module HydroModels

## External packages
# common packages
using Accessors
using Reexport

@reexport using ComponentArrays
using ComponentArrays: indexmap, getval
using Dates
using DataFrames
using IterTools: ncycle
using LinearAlgebra
using NamedTupleTools
using ProgressMeter
using SparseArrays
using StableRNGs
using Statistics
using TOML

# runtime generated functions
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

# Symbolic building
using Symbolics
using Symbolics: tosymbol
using SymbolicUtils
using SymbolicUtils.Code
@reexport using ModelingToolkit: @variables, @parameters
using ModelingToolkit: isparameter
using ModelingToolkit: t_nounits as t
# graph compute
using Graphs

# data interpolataion
using DataInterpolations
using DataInterpolations: AbstractInterpolation

# integral
using Integrals

# solve ODEProblem
using SciMLBase
using OrdinaryDiffEq
using SciMLSensitivity

# deep learning
using Lux
using LuxCore
using NNlib

# parameters Optimization
using Optimization
using OptimizationBBO
using OptimizationOptimisers

## package version
const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])
const InputType = Union{AbstractArray,AbstractMatrix}

struct HydroEquation end
## Abstract Component Types
abstract type AbstractComponent end
abstract type AbstractHydroSolver end
abstract type AbstractHydroOptimizer end
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
abstract type AbstractRapidRoute <: AbstractRoute end
abstract type AbstractModel <: AbstractComponent end

export AbstractComponent, AbstractHydroSolver, AbstractHydroOptimizer, AbstractHydroWrapper, AbstractNeuralWrapper
export AbstractFlux, AbstractHydroFlux, AbstractNeuralFlux, AbstractStateFlux
export AbstractElement, AbstractBucket, AbstractHydrograph, AbstractRoute, AbstractHydroRoute, AbstractRapidRoute, AbstractModel

# utils
include("utils/attr.jl")
include("utils/ca.jl")
include("utils/name.jl")
include("utils/show.jl")
include("utils/build.jl")
include("utils/callback.jl")
include("utils/sort.jl")
include("utils/io.jl")
export NamedTupleIOAdapter

include("optimizer.jl")
export BatchOptimizer, HydroOptimizer, GradOptimizer

include("solver.jl")
export ODESolver, DiscreteSolver, ManualSolver

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
include("wrapper.jl")
export RecordComponentState, EstimateComponentParams, WeightSumComponentOutlet, ComputeComponentOutlet
end # module HydroModels
