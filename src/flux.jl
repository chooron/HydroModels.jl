"""
Flux module - defines hydrological flux components, supporting both symbolic and functional construction approaches.
"""

"""
    HydroFlux{MS, E, F, HT, I} <: AbstractHydroFlux

Represents a simple flux component defined by mathematical expressions or Julia functions.

It automatically parses expressions to determine inputs, outputs, and parameters, and compiles them into efficient callable functions.

$(FIELDS)

# Type Parameters
- `MS`: Whether in multi-node mode
- `E`: Expression type
- `F`: Compiled function type
- `HT`: HRU type vector type
- `I`: Metadata type
"""
struct HydroFlux{MS,E,F,HT,I} <: AbstractHydroFlux
    "flux name"
    name::Symbol
    "Vector of expressions describing the formulas for output variables"
    exprs::E
    "Compiled function that calculates the flux"
    func::F
    "nodes type"
    hru_types::HT
    "Metadata about the flux, including input, output, and parameter names"
    infos::I
    
    # Symbolic constructor
    function HydroFlux(;
        exprs::E,
        hru_types::Vector{Int}=Int[],
        name::Optional{Symbol}=nothing,
    ) where {E}
        # Parse expressions and extract variables
        outputs = Num.([eq.lhs for eq in exprs])
        eqs = Num.([eq.rhs for eq in exprs])
        all_vars = Num.(mapreduce(get_variables, union, eqs, init=Set{Num}()))
        inputs = setdiff(Num.(filter(x -> !isparameter(x), collect(all_vars))), outputs)
        params = Num.(filter(x -> isparameter(x), collect(all_vars)))
        
        @assert length(exprs) == length(outputs) "Number of expressions must match number of outputs"
        
        # Build flux calculation function
        infos = HydroInfos(
            inputs=!isempty(inputs) ? tosymbol.(inputs) : Symbol[],
            params=!isempty(params) ? tosymbol.(params) : Symbol[],
            outputs=!isempty(outputs) ? tosymbol.(outputs) : Symbol[]
        )
        flux_func = build_flux_func(eqs, infos)
        flux_name = isnothing(name) ? Symbol("##hydro_flux#", hash(infos)) : name
        
        return new{length(hru_types) > 1,typeof(eqs),typeof(flux_func),typeof(hru_types),typeof(infos)}(
            flux_name, eqs, flux_func, hru_types, infos
        )
    end
    
    # Functional constructor - directly use Julia functions
    function HydroFlux(
        func::Function;
        inputs::Vector{Symbol},
        outputs::Vector{Symbol},
        params::Vector{Symbol}=Symbol[],
        hru_types::Vector{Int}=Int[],
        name::Optional{Symbol}=nothing,
    )
        infos = HydroInfos(
            inputs=inputs,
            params=params,
            outputs=outputs
        )
        flux_name = isnothing(name) ? Symbol("##hydro_flux_func#", hash(infos)) : name
        
        # Wrap user function to match internal interface
        wrapped_func = (inp_slices, pas) -> begin
            # inp_slices: iterator of input slices
            # pas: ComponentVector parameters
            # Use Tuple instead of collect for type stability and zero allocation
            inp_vals = Tuple(inp_slices)
            param_vals = haskey(pas, :params) ? pas[:params] : NamedTuple()
            func(inp_vals, param_vals)
        end
        
        return new{length(hru_types) > 1,Nothing,typeof(wrapped_func),typeof(hru_types),typeof(infos)}(
            flux_name, nothing, wrapped_func, hru_types, infos
        )
    end
end

"""
    @hydroflux(name, eqs...)

Macro to conveniently create a HydroFlux component from a set of equations.

It parses the given equations to identify output variables (left-hand sides), input variables, 
and parameters (right-hand sides). Parameters are distinguished from variables using 
ModelingToolkit's `isparameter` function.

# Arguments
- `name`: Optional Symbol to name the flux component
- `eqs...`: One or more equations defining the flux

# Examples
```jldoctest
julia> @variables P, ET, Q
julia> @parameters k
julia> @hydroflux :runoff Q ~ k * (P - ET)
```
"""
macro hydroflux(args...)
    name = length(args) == 1 ? nothing : args[1]
    eqs_expr = length(args) == 1 ? args[1] : args[2]
    
    # Handle begin...end blocks
    vect_eqs_expr = if Meta.isexpr(eqs_expr, :block)
        Expr(:tuple, filter(arg -> !(arg isa LineNumberNode) && !Meta.isexpr(arg, :line) && !Meta.isexpr(arg, :(=)), eqs_expr.args)...)
    else
        Expr(:tuple, eqs_expr)
    end
    
    # Extract hru_types (if any)
    hru_types_val = :(Int[])
    if Meta.isexpr(eqs_expr, :block)
        for assign in filter(x -> Meta.isexpr(x, :(=)), eqs_expr.args)
            lhs, rhs = assign.args
            if lhs == :hru_types
                hru_types_val = rhs
            end
        end
    end
    
    # Check variable definitions
    for var_name in extract_variables(vect_eqs_expr)
        if !@isdefined(var_name)
            expr_str = string(vect_eqs_expr)
            return :(error("Undefined variable '", $(string(var_name)), "' detected in expression: `", $expr_str, "`"))
        end
    end
    
    return esc(:(HydroFlux(exprs=$vect_eqs_expr, name=$name, hru_types=$hru_types_val)))
end

# Flux computation for 2D input (single node)
function (flux::HydroFlux)(
    input::AbstractArray{T,2},
    params::ComponentVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,2} where {T}
    result = flux.func(eachslice(input, dims=1), params)
    stack(result, dims=1)
end

# Flux computation for 3D input (multi-node)
function (flux::HydroFlux{true})(
    input::AbstractArray{T,3},
    params::ComponentVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,3} where {T}
    expand_params = expand_component_params(params, get_param_names(flux), flux.hru_types)
    output = flux.func(eachslice(input, dims=1), expand_params)
    stack(output, dims=1)
end

"""
    StateFlux{N, E, I} <: AbstractStateFlux

Represents a state flux, defining the rate of change (derivative) for a state variable.

This is a declarative component used within a HydroBucket or HydroRoute to define 
the system's differential equations. It does not perform calculations itself.

$(FIELDS)

# Examples
```jldoctest
julia> @variables S, P, Q
julia> @stateflux :storage_change S ~ P - Q
```
"""
struct StateFlux{N,E,I} <: AbstractStateFlux
    "Vector of expressions defining the state's rate of change"
    exprs::E
    "Metadata about the component, including input, state, and parameter names"
    infos::I
    
    function StateFlux(;
        exprs::E,
        name::Optional{Symbol}=nothing
    ) where {E}
        states = Num.([eq.lhs for eq in exprs])
        eqs = Num.([eq.rhs for eq in exprs])
        all_vars = Num.(mapreduce(get_variables, union, eqs, init=Set{Num}()))
        inputs = setdiff(Num.(filter(x -> !isparameter(x), collect(all_vars))), states)
        params = Num.(filter(x -> isparameter(x), collect(all_vars)))
        
        infos = HydroInfos(
            inputs=!isempty(inputs) ? tosymbol.(inputs) : Symbol[],
            states=!isempty(states) ? tosymbol.(states) : Symbol[],
            params=!isempty(params) ? tosymbol.(params) : Symbol[]
        )
        flux_name = isnothing(name) ? Symbol("##state_flux#", hash(infos)) : name
        
        return new{flux_name,typeof(eqs),typeof(infos)}(eqs, infos)
    end
end

"""
    @stateflux [name] eq

Macro to conveniently create a StateFlux component from a single equation.

# Usage
The left side of the equation is the state variable, and the right side is the expression 
for its rate of change. The `~` operator is typically used.

```jldoctest
julia> @variables S, P, Q
julia> @stateflux :storage_change S ~ P - Q
```
"""
macro stateflux(args...)
    name = length(args) == 1 ? nothing : args[1]
    eq_expr = length(args) == 1 ? args[1] : args[2]
    
    for var_name in extract_variables(eq_expr)
        if !@isdefined(var_name)
            expr_str = string(eq_expr)
            return :(error("Undefined variable '", $(string(var_name)), "' detected in expression: `", $expr_str, "`"))
        end
    end
    
    return esc(:(StateFlux(exprs=[$eq_expr], name=$name)))
end

# StateFlux cannot be run directly
(::StateFlux)(::AbstractArray, ::ComponentVector; kwargs...) = 
    error("StateFlux cannot be run directly, please use HydroBucket or HydroRoute")

# Export interfaces
export HydroFlux, StateFlux
export @hydroflux, @stateflux
