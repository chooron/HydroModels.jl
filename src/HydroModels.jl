"""
    HydroModels

Modern hydrological modeling framework supporting both symbolic and functional modeling approaches.

# Key Features
- Type-stable configuration system
- Fully compatible with Zygote and Enzyme automatic differentiation
- Flexible symbolic and functional modeling interfaces
- Unified single-node (2D) and multi-node (3D) components via `htypes`
- Pluggable interpolators (supports DataInterpolations.jl via extension)
- Optimized ODE solvers

# Core Components
- `HydroFlux`: Flux calculation component (2D/3D via htypes)
- `StateFlux`: State derivative component
- `NeuralFlux`: Neural network-driven flux component
- `HydroBucket`: Hydrological bucket container (2D/3D via htypes)
- `NeuralBucket`: Neural network-based bucket (RNN-style, embeds into HydroModel)
- `HydroRoute`: Spatial routing component (multi-node only)
- `UnitHydrograph`: Unit hydrograph convolution component
- `HydroModel`: Complete hydrological model

# Basic Usage

```julia
using HydroModels

# Define symbolic variables
@variables P, ET, Q, S
@parameters Smax, f, Qmax

# Create flux component
flux1 = @hydroflux begin
    Q ~ Qmax * (S / Smax)^f
end

# Create state flux
dflux1 = @stateflux S ~ P - ET - Q

# Create bucket
bucket = @hydrobucket :my_bucket begin
    fluxes = begin
        flux1
    end
    dfluxes = begin
        dflux1
    end
end

# Create model
model = @hydromodel :my_model begin
    bucket
end
```
"""
module HydroModels

# ============================================================================
# External Dependencies
# ============================================================================

# Common packages
using Random
using Reexport
using SpecialFunctions
using ComponentArrays
using ComponentArrays: getaxes, Axis
using DocStringExtensions

# Graph computation
using Graphs
using Graphs: SimpleDiGraph, add_edge!, topological_sort, adjacency_matrix

# Deep learning
using Lux
using LuxCore
using NNlib

# ============================================================================
# Type Definitions
# ============================================================================

"""Optional type: value or Nothing"""
const Optional{T} = Union{T,Nothing}

"""Solver type enumeration"""
Base.@enum SolverType begin
    MutableSolver      # Mutable solver (efficient)
    ImmutableSolver    # Immutable solver (Zygote-friendly)
    ODESolver          # ODE solver (interface reserved)
    DiscreteSolver     # Discrete solver
end

# ============================================================================
# Core Dependencies
# ============================================================================

using HydroModelCore
using HydroModelCore: AbstractComponent, AbstractFlux, AbstractStateFlux
using HydroModelCore: AbstractHydroFlux, AbstractNeuralFlux
using HydroModelCore: AbstractHydroBucket, AbstractHydroRoute, AbstractHydrograph
using HydroModelCore: AbstractModel
using HydroModelCore: HydroInfos
using HydroModelCore: get_var_names, get_input_names, get_output_names
using HydroModelCore: get_state_names, get_param_names, get_nn_names
using HydroModelCore: isparameter, toparam
using HydroModelCore: build_flux_func, build_bucket_func, build_route_func, build_uh_func
using HydroModelCore: tosymbol, Num, @variables, get_variables

# ============================================================================
# Export Symbolic Computation Tools
# ============================================================================

export @variables, @parameters
export Num, tosymbol, get_variables, isparameter, toparam

# ============================================================================
# Export Solver Types
# ============================================================================

export SolverType, MutableSolver, ImmutableSolver, ODESolver, DiscreteSolver

# ============================================================================
# Module Includes
# ============================================================================

# Error types
include("errors.jl")
export HydroModelsError, ConfigurationError, DimensionMismatchError
export MacroSyntaxError, ParameterError, VariableError

# Configuration system
include("config.jl")
export HydroConfig, default_config, merge_config, normalize_config
export get_config_value, ConfigType

# Interpolation
include("interpolate.jl")
export ConstantInterpolation, LinearInterpolation
export DirectInterpolation, EnzymeInterpolation, EnzymeCompatibleInterpolation  # backward compat
export hydrointerp

# Solver and utilities
include("solve.jl")
export hydrosolve, safe_clamp, safe_max

# Utility functions
include("utils.jl")
export sort_components, sort_fluxes
export expand_component_params, get_default_states
export extract_variables

# Flux components
include("flux.jl")
# exports: HydroFlux, HydroMultiFlux (alias), StateFlux, @hydroflux, @hydromultiflux (deprecated), @stateflux

# Neural network components
include("nn.jl")
# exports: NeuralFlux, NeuralBucket, @neuralflux, create_neural_bucket, create_simple_neural_bucket

# Bucket components
include("bucket.jl")
# exports: HydroBucket, HydroMultiBucket (alias), @hydrobucket, @hydromultibucket (deprecated)

# Routing components
include("route.jl")
export ChannelRoute, @channelroute
export HydroRoute, @hydroroute, build_aggr_func
export RouteIRF, SparseRouteKernel, SparseRouteConvolution
export build_irf_kernels, aggregate_route_kernel

# Unit hydrograph components
include("uh.jl")
export UnitHydrograph, @unithydro


# Model components
include("model.jl")
export HydroModel, @hydromodel

# ============================================================================
# Common Hydrological Functions
# ============================================================================

"""
    step_func(x)

Smooth step function implemented using tanh.

# Arguments
- `x`: Input value

# Returns
- Smooth step value in range [0, 1]

# Examples
```jldoctest
julia> y = step_func(0.0)  # ~0.5
julia> y = step_func(2.0)  # ~1.0
julia> y = step_func(-2.0) # ~0.0
```
"""
@inline step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

"""
    smoothlogistic_func(S, Smax; r=0.01, e=5.0)

Smooth logistic function for simulating soil moisture changes.

# Arguments
- `S`: Current state value
- `Smax`: Maximum state value
- `r`: Smoothing parameter (default 0.01)
- `e`: Exponential parameter (default 5.0)

# Returns
- Smooth logistic value

# Examples
```jldoctest
julia> y = smoothlogistic_func(50.0, 100.0)
```
"""
@inline function smoothlogistic_func(S, Smax; r=0.01, e=5.0)
    1.0 / (1.0 + exp((S - r * e * Smax) / (r * Smax)))
end

# Export hydrological functions
export step_func, smoothlogistic_func

# ============================================================================
# YAML Extension Function Stubs
# ============================================================================

"""
    load_model_from_yaml(yaml_file::String)

Load a hydrological model from a YAML configuration file.

This function is provided by the HydroModelsYAMLExt extension.
To use it, you must first load the YAML package:

```julia
using HydroModels
using YAML  # This triggers the extension

model = load_model_from_yaml("model.yaml")
```

# Arguments
- `yaml_file`: Path to YAML configuration file

# Returns
- HydroModel object
"""
function load_model_from_yaml end

"""
    load_config_from_yaml(yaml_file::String)

Load configuration from a YAML file.

This function is provided by the HydroModelsYAMLExt extension.
To use it, you must first load the YAML package.

# Arguments
- `yaml_file`: Path to YAML configuration file

# Returns
- HydroConfig object
"""
function load_config_from_yaml end

"""
    load_parameters_from_yaml(yaml_file::String)

Load parameter metadata from a YAML file.

This function is provided by the HydroModelsYAMLExt extension.
To use it, you must first load the YAML package.

# Arguments
- `yaml_file`: Path to YAML configuration file

# Returns
- Dictionary of parameter metadata
"""
function load_parameters_from_yaml end

"""
    execute_from_yaml(yaml_file::String; return_components::Bool=false)

Execute a hydrological model from a YAML configuration file.

This function is provided by the HydroModelsYAMLExt extension.
To use it, you must first load the YAML package:

```julia
using HydroModels
using YAML  # This triggers the extension

# Execute model
output = execute_from_yaml("model_config.yaml")

# Or get all components
output, model, config, data = execute_from_yaml("model_config.yaml", return_components=true)
```

# Arguments
- `yaml_file`: Path to YAML configuration file
- `return_components`: If true, return (output, model, config, data) tuple

# Returns
- Model output, or tuple of (output, model, config, data) if return_components=true
"""
function execute_from_yaml end

# Export YAML functions (will be extended by HydroModelsYAMLExt)
export load_model_from_yaml, load_config_from_yaml, load_parameters_from_yaml, execute_from_yaml

# ============================================================================
# Module Version Information
# ============================================================================

"""
Module version: v0.6.2

# Major Updates (v0.6.2)
- Unified HydroFlux/HydroMultiFlux into single HydroFlux (htypes dispatch)
- Unified HydroBucket/HydroMultiBucket into single HydroBucket (htypes dispatch)
- Redesigned NeuralBucket to implement AbstractHydroBucket (embeds into HydroModel)
- Unified UnitHydrograph htypes pattern
- Added DataInterpolations.jl extension for pluggable interpolators
- Backward-compatible aliases: HydroMultiFlux, HydroMultiBucket, @hydromultiflux, @hydromultibucket

# Julia Version Requirements
- Julia >= 1.10 (recommended 1.12+)
"""
const HYDROMODELS_VERSION = v"0.6.2"

end # module HydroModels








