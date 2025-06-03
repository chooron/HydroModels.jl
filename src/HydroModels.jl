module HydroModels

## External packages
# common packages
using Reexport
using LinearAlgebra
using SparseArrays
using Random

# runtime generated functions
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

# Symbolic building
using Symbolics
using Symbolics: tosymbol, unwrap, wrap, Num, Symbolic, @variables, get_variables
using SymbolicUtils.Code
import SymbolicUtils: symtype, term, hasmetadata, issym

# graph compute
using Graphs

# deep learning
using Lux
using NNlib

# SciML ecosystem
@reexport using ComponentArrays
@reexport using DataInterpolations
using OrdinaryDiffEq
using SciMLSensitivity

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

export AbstractComponent, # base type
    AbstractHydroFlux, AbstractNeuralFlux, AbstractStateFlux, # flux types
    AbstractElement, # element types
    AbstractBucket, # bucket types
    AbstractHydrograph, # hydrograph types
    AbstractRoute, AbstractHydroRoute, # route types
    AbstractModel # model types

export get_name, get_input_names, get_output_names, get_param_names, get_state_names, get_nn_names
export get_exprs, get_inputs, get_outputs, get_params, get_nns, get_vars

# utils
include("utils/parameters.jl") # parameters.jl in ModelingToolkit
export @variables, @parameters

include("utils/attribute.jl")
include("utils/display.jl")
include("utils/check.jl")
include("utils/aggregation.jl")
include("utils/miscellaneous.jl")
include("utils/build.jl")
include("utils/smooth.jl")
export step_func, smoothlogistic_func
#! When constructing an ODE problem to solve, use DataInterpolations.jl
include("tools.jl")
export ManualSolver, ODESolver, DiscreteSolver, DirectInterpolation

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

AVAILABLE_MODELS = [
    :alpine1,
    :alpine2,
    :australia,
    :collie1,
    :collie2,
    :collie3,
    :echo,
    :exphydro,
    :flexb,
    :gr4j,
    :gsfb,
    :gsmsocont,
    :hbv_edu,
    :hbv,
    :hillslope,
    :hymod,
    :ihacres,
    :ihm19,
    :lascam,
    :mcrm,
    :modhydrolog,
    :mopex1,
    :mopex3,
    :mopex4,
    :mopex5,
    :nam,
    :newzealand1,
    :newzealand2,
    :penman,
    :plateau,
    :prms,
    :sacramento,
    :smar,
    :smart,
    :susannah1,
    :susannah2,
    :tank,
    :tcm,
    :unitedstates,
    :wetland,
    :xinanjiang
]

AVAILABLE_ROUTERS = [
    :rapid
]
export AVAILABLE_MODELS, AVAILABLE_ROUTERS

function load_model(model_name::Symbol)
    if model_name in AVAILABLE_MODELS
        if !isdefined(HydroModels, model_name)
            model_path = joinpath(@__DIR__, "models", "$(model_name).jl")
            include(model_path)
        end
        return getfield(getfield(HydroModels, model_name), :model)
    else
        throw(ArgumentError("Model $model_name is not available"))
    end
end

function load_router(router_name::Symbol)
    if router_name in AVAILABLE_ROUTERS
        if !isdefined(HydroModels, router_name)
            router_path = joinpath(@__DIR__, "routers", "$(router_name).jl")
            include(router_path)
        end
        return getfield(getfield(HydroModels, router_name), :router)
    else
        throw(ArgumentError("Router $router_name is not available"))
    end
end
export load_model, load_router

end # module HydroModels
