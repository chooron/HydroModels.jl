module HydroModelCore

using DocStringExtensions
using ComponentArrays
using Printf
using Symbolics
using Symbolics: tosymbol, unwrap, wrap, Num, @variables, get_variables, symbolic_type, ArraySymbolic
using SymbolicUtils.Code
import SymbolicUtils: symtype, term, hasmetadata, issym, BasicSymbolic

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

abstract type AbstractComponent end
abstract type AbstractNetwork end
abstract type AbstractConfig end
abstract type AbstractInfos end

abstract type AbstractFlux <: AbstractComponent end
abstract type AbstractHydroFlux <: AbstractFlux end
abstract type AbstractNeuralFlux <: AbstractFlux end
abstract type AbstractStateFlux <: AbstractFlux end

abstract type AbstractElement <: AbstractComponent end
abstract type AbstractBucket <: AbstractElement end
abstract type AbstractHydroBucket <: AbstractBucket end
abstract type AbstractNeuralBucket <: AbstractBucket end
abstract type AbstractRoute <: AbstractElement end
abstract type AbstractHydroRoute <: AbstractRoute end
abstract type AbstractHydrograph <: AbstractRoute end
abstract type AbstractModel <: AbstractComponent end


export AbstractComponent, AbstractConfig, AbstractInfos, # base type
    AbstractFlux, AbstractHydroFlux, AbstractNeuralFlux, AbstractStateFlux, # flux types
    AbstractElement, # element types
    AbstractBucket, AbstractHydroBucket, AbstractNeuralBucket, # bucket types
    AbstractRoute, AbstractHydroRoute, AbstractHydrograph, # route types
    AbstractModel, # model types
    AbstractNetwork # network types

include("attribute.jl")
include("check.jl")
include("display.jl")
include("parameters.jl")
include("variables.jl")
include("build.jl")  # Unified build system for ForwardDiff/Mooncake workflows
include("symbol_build.jl")

# ================================================================================================
# Core data structures
# ================================================================================================
export HydroInfos

# ================================================================================================
# Attribute accessor functions (attribute.jl)
# ================================================================================================
# Basic accessors - get names from components
export get_name           # Get component name
export get_input_names, get_output_names, get_param_names, get_state_names, get_nn_names
export get_var_names      # Get (inputs, outputs, states) tuple
export get_all_names      # Get (inputs, outputs, states, params, nns) tuple
export get_exprs          # Get expressions from flux components

# Attribute counting - count number of variables/parameters
export count_inputs, count_outputs, count_states, count_params, count_nns

# Existence checking - check if component has attributes
export has_inputs, has_outputs, has_states, has_params, has_nns

# Batch collection - collect unique attributes from multiple components
export collect_all_inputs, collect_all_outputs, collect_all_states
export collect_all_params, collect_all_nns

# ================================================================================================
# Validation functions (check.jl)
# ================================================================================================
# Core validation functions
export check              # Main validation function
export check_input        # Validate input dimensions
export check_params       # Validate parameters exist
export check_initstates   # Validate initial states exist
export check_nns          # Validate neural networks exist

# Value validation functions
export check_param_values # Check parameter values (NaN/Inf/negative)
export check_state_values # Check state values (NaN/Inf)

# Batch validation functions
export check_multiple     # Check multiple components
export check_all          # Check all components pass validation

# Dependency validation functions
export check_dependencies      # Check component dependencies
export check_component_chain   # Check component chain validity

# ================================================================================================
# Display functions (display.jl)
# ================================================================================================
# Pretty printing utility functions
export print_component_table    # Print summary table of components
export print_variable_details   # Print detailed variable information

# ================================================================================================
# Variable system (variables.jl)
# ================================================================================================
export @variables, @parameters, isparameter
export tosymbol, Num, get_variables
export getdescription, getbounds, getunit, getguess

# ================================================================================================
# Build system (build.jl)
# ================================================================================================
# Core build functions
export build_flux_func, build_bucket_func, build_route_func, build_uh_func
export build_symbolic_flux_func, build_symbolic_bucket_func, build_symbolic_route_func, build_symbolic_uh_func

# Build configuration types and constants
export Dim0, Dim1, Dim2, AbstractDimConfig
export BuildConfig, PerformanceMode, Safe, Fast, AutoDiff
export SAFE_CONFIG, FAST_CONFIG, AUTODIFF_CONFIG, DEBUG_CONFIG

# Build system utility functions
export preview_function_expr, analyze_generated_code

end # module HydroModelCore
