module HydroModels

## External packages
# common packages
using Reexport
@reexport using ComponentArrays
using LinearAlgebra
using SparseArrays
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

# SciML ecosystem
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

export AbstractComponent, AbstractHydroFlux, AbstractNeuralFlux, AbstractStateFlux,
    AbstractElement, AbstractBucket, AbstractHydrograph,
    AbstractRoute, AbstractHydroRoute, AbstractModel

export get_name, get_input_names, get_output_names, get_param_names, get_state_names, get_nn_names
export get_exprs, get_inputs, get_outputs, get_params, get_nns, get_vars

# utils
include("utils/attribute.jl")
include("utils/display.jl")
include("utils/check.jl")
include("utils/aggregation.jl")
include("utils/expression.jl")
include("utils/miscellaneous.jl")
include("utils/build.jl")
include("utils/smooth.jl")
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
export AVAILABLE_MODELS

# 定义一个函数来按需加载模型
function load_model(model_name::Symbol)
    if model_name in AVAILABLE_MODELS
        # 检查模块是否已经加载，如果没有则加载
        if !isdefined(HydroModels, model_name)
            model_path = joinpath(@__DIR__, "models", "$(model_name).jl")
            include(model_path)
        end
        # 返回模块
        return getfield(getfield(HydroModels, model_name), :model)
    else
        throw(ArgumentError("Model $model_name is not available"))
    end
end

## package version
const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

end # module HydroModels
