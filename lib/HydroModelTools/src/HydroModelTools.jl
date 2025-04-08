module HydroModelTools

# base dependencies
using Reexport
using Dates
using DataFrames
using ProgressMeter
using IterTools: ncycle
using NamedTupleTools
using Statistics
using TOML
@reexport using ComponentArrays
using ComponentArrays: indexmap, getval

# solve ODEProblem
using OrdinaryDiffEq
using SciMLSensitivity

# Optimization algorithms
using Optimization
using OptimizationBBO
using OptimizationOptimisers

const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

abstract type AbstractHydroSolver end
abstract type AbstractHydroOptimizer end

abstract type AbstractHydroWrapper end
abstract type AbstractDataPreprocessor <: AbstractHydroWrapper end
abstract type AbstractDataPostprocessor <: AbstractHydroWrapper end
abstract type AbstractComponentDecorator <: AbstractHydroWrapper end

include("utils/attribute.jl")
include("utils/ca.jl")
export update_ca, merge_ca
include("utils/callback.jl")
export get_callback_func, get_batch_callback_func
include("wrappers/namedtuple_wrapper.jl")
export NamedTuplePreprocessor, NamedTuplePostprocessor
include("wrappers/stats_outlet.jl")
export SelectComponentOutlet

# include("wappers/neural_wrapper.jl")
# export NeuralWrapper, NeuralWrapper3

# include("wappers/estimate_params.jl")
# export EstimateComponentParams

# include("wappers/record_states.jl")
# export RecordStates

# include("wappers/stats_outlet.jl")
# export StatsOutlet

include("optimizer.jl")
export BatchOptimizer, HydroOptimizer, GradOptimizer

include("solver.jl")
export ODESolver, DiscreteSolver, ManualSolver

end # module HydroTools
