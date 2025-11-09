"""
    HydroModels

Modern hydrological modeling framework supporting both symbolic and functional modeling approaches.

# Key Features
- Type-stable configuration system
- Fully compatible with Zygote automatic differentiation
- Flexible symbolic and functional modeling interfaces
- Efficient computation without in-place modifications
- Support for single-node and multi-node spatial simulations
- Optimized ODE solvers

# Core Components
- `HydroFlux`: Flux calculation component
- `StateFlux`: State derivative component
- `NeuralFlux`: Neural network-driven flux component
- `HydroBucket`: Hydrological bucket container
- `NeuralBucket`: RNN-style neural bucket component
- `HydroRoute`: Spatial routing component
- `UnitHydrograph`: Unit hydrograph convolution component
- `HydroModel`: Complete hydrological model

# Basic Usage

```julia
using HydroModels
using Symbolics

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

# Functional Interface

```julia
# Define flux using pure Julia function
flux_func = HydroFlux(
    (inputs, params) -> [params.k * inputs[1]];
    inputs=[:P],
    outputs=[:Q],
    params=[:k],
    name=:simple_flux
)
```

# Configuration System

```julia
# Create configuration
config = HydroConfig(
    solver=MutableSolver,
    interpolator=Val(DirectInterpolation),
    min_value=1e-6,
    parallel=false
)

# Run model
output = model(input, params, config)
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

# Symbolic computation
using Symbolics
using Symbolics: tosymbol, Num, @variables, get_variables

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
using HydroModelCore: isparameter
using HydroModelCore: build_flux_func, build_bucket_func, build_route_func, build_uh_func

# ============================================================================
# Export Symbolic Computation Tools
# ============================================================================

export @variables, @parameters
export Num, tosymbol, get_variables, isparameter

# ============================================================================
# Export Solver Types
# ============================================================================

export SolverType, MutableSolver, ImmutableSolver, ODESolver, DiscreteSolver

# ============================================================================
# Module Includes
# ============================================================================

# Configuration system
include("config.jl")
export HydroConfig, default_config, merge_config, normalize_config
export get_config_value, ConfigType

# Tool functions
include("tools.jl")
export DirectInterpolation, hydrointerp, hydrosolve
export safe_clamp, safe_max

# Utility functions
include("utils.jl")
export sort_components, sort_fluxes
export expand_component_params, get_default_states
export extract_variables
export @hydroflux_for

# Flux components
include("flux.jl")
export HydroFlux, StateFlux
export @hydroflux, @stateflux

# Neural network components
include("nn.jl")
export NeuralFlux, NeuralBucket
export @neuralflux
export create_neural_bucket, create_simple_neural_bucket

# Bucket components
include("bucket.jl")
export HydroBucket, @hydrobucket

# Routing components
include("route.jl")
export HydroRoute, @hydroroute, build_aggr_func

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
julia> y = step_func(0.0)  # ≈ 0.5
julia> y = step_func(2.0)  # ≈ 1.0
julia> y = step_func(-2.0) # ≈ 0.0
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
# Module Version Information
# ============================================================================

"""
Module version: v0.4.5

# Major Updates
- ✨ Complete refactoring with type-stable design
- 🚀 100% Zygote automatic differentiation compatibility
- 🎯 Support for both functional and symbolic interfaces
- 📈 Significant performance improvements
- 🧹 Removed redundant code
- 📚 Improved documentation and type annotations

# Backward Compatibility
- ✅ All original computational interfaces preserved
- ✅ NamedTuple configuration supported (auto-conversion)
- ✅ Original macros and symbolic interfaces fully compatible

# Julia Version Requirements
- Julia >= 1.10 (recommended 1.12+)
"""
const HYDROMODELS_VERSION = v"0.4.5"

end # module HydroModels
