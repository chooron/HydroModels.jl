"""
Flux module - defines hydrological flux components, supporting both symbolic and functional construction approaches.
"""

"""
    HydroFlux{E, F, HT, I} <: AbstractHydroFlux

Represents a flux component defined by mathematical expressions or Julia functions.

When `htypes` is `nothing`, this is a single-node component that accepts 2D input (variables × time).
When `htypes` is a `Vector{Int}`, this is a multi-node component that accepts 3D input (variables × nodes × time)
with parameter sharing via htypes indexing.

$(FIELDS)

# Type Parameters
- `E`: Expression type
- `F`: Compiled function type
- `HT`: HRU type (`Nothing` for 2D, `Vector{Int}` for 3D)
- `I`: Metadata type
"""
struct HydroFlux{E,F,HT,I} <: AbstractHydroFlux
    "flux name"
    name::Symbol
    "Vector of expressions describing the formulas for output variables"
    exprs::E
    "Compiled function that calculates the flux"
    func::F
    "HRU type indices for parameter sharing (Nothing = single-node 2D, Vector{Int} = multi-node 3D)"
    htypes::HT
    "Metadata about the flux, including input, output, and parameter names"
    infos::I

    # Symbolic constructor
    function HydroFlux(;
        exprs::E,
        name::Optional{Symbol}=nothing,
        htypes::Optional{Vector{Int}}=nothing,
    ) where {E}
        outputs = Num.([eq.lhs for eq in exprs])
        eqs = Num.([eq.rhs for eq in exprs])
        all_vars = Num.(mapreduce(get_variables, union, eqs, init=Set{Num}()))
        inputs = setdiff(Num.(filter(x -> !isparameter(x), collect(all_vars))), outputs)
        params = Num.(filter(x -> isparameter(x), collect(all_vars)))
        @assert length(exprs) == length(outputs) "Number of expressions must match number of outputs"

        infos = HydroInfos(
            inputs=!isempty(inputs) ? tosymbol.(inputs) : Symbol[],
            params=!isempty(params) ? tosymbol.(params) : Symbol[],
            outputs=!isempty(outputs) ? tosymbol.(outputs) : Symbol[]
        )
        flux_func = build_flux_func(eqs, infos)
        flux_name = isnothing(name) ? Symbol("##hydro_flux#", hash(infos)) : name

        return new{typeof(eqs),typeof(flux_func),typeof(htypes),typeof(infos)}(
            flux_name, eqs, flux_func, htypes, infos
        )
    end

    # Functional constructor - directly use Julia functions
    function HydroFlux(
        func::Function;
        inputs::Vector{Symbol},
        outputs::Vector{Symbol},
        params::Vector{Symbol}=Symbol[],
        name::Optional{Symbol}=nothing,
        htypes::Optional{Vector{Int}}=nothing,
    )
        infos = HydroInfos(
            inputs=inputs,
            params=params,
            outputs=outputs
        )
        flux_name = isnothing(name) ? Symbol("##hydro_flux_func#", hash(infos)) : name

        wrapped_func = (inp_slices, pas) -> begin
            inp_vals = Tuple(inp_slices)
            param_vals = haskey(pas, :params) ? pas[:params] : NamedTuple()
            func(inp_vals, param_vals)
        end

        return new{Nothing,typeof(wrapped_func),typeof(htypes),typeof(infos)}(
            flux_name, nothing, wrapped_func, htypes, infos
        )
    end
end

# ============================================================================
# @hydroflux macro - unified for single-node and multi-node
# ============================================================================

"""
    @hydroflux [name] expr

Macro to conveniently create a HydroFlux component from equations.

When `htypes` is specified in the begin block, creates a multi-node flux.
Otherwise creates a single-node flux.

# Examples
```julia
# Single-node flux
@hydroflux :runoff Q ~ k * (P - ET)

# Multi-node flux with parameter sharing
@hydroflux :runoff begin
    Q ~ k * (P - ET)
    htypes = [1, 1, 2, 2, 3]
end
```
"""
macro hydroflux(args...)
    name = length(args) == 1 ? nothing : args[1]
    eqs_expr = length(args) == 1 ? args[1] : args[2]

    htypes_val = nothing

    vect_eqs_expr = if Meta.isexpr(eqs_expr, :block)
        eqs_items = []
        for item in filter(x -> !(x isa LineNumberNode) && !Meta.isexpr(x, :line), eqs_expr.args)
            if Meta.isexpr(item, :(=))
                lhs, rhs = item.args
                if lhs == :htypes
                    htypes_val = rhs
                else
                    push!(eqs_items, item)
                end
            else
                push!(eqs_items, item)
            end
        end
        Expr(:tuple, eqs_items...)
    else
        Expr(:tuple, eqs_expr)
    end

    # Check variable definitions
    for var_name in extract_variables(vect_eqs_expr)
        if !@isdefined(var_name)
            expr_str = string(vect_eqs_expr)
            return :(error("Undefined variable '", $(string(var_name)), "' detected in expression: `", $expr_str, "`"))
        end
    end

    if isnothing(htypes_val)
        return esc(:(HydroFlux(exprs=$vect_eqs_expr, name=$name)))
    else
        return esc(:(HydroFlux(exprs=$vect_eqs_expr, htypes=$htypes_val, name=$name)))
    end
end

# ============================================================================
# Functor methods - dispatch on HT type parameter
# ============================================================================

# 2D computation (single-node, htypes = Nothing)
function (flux::HydroFlux{E,F,Nothing,I})(
    input::AbstractArray{T,2},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,2} where {E,F,I,T}
    params = _as_componentvector(params)
    result = flux.func(eachslice(input, dims=1), params)
    stack(result, dims=1)
end

# Error: single-node flux receiving 3D input
function (flux::HydroFlux{E,F,Nothing,I})(
    input::AbstractArray{T,3},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
) where {E,F,I,T}
    error("HydroFlux without htypes only accepts 2D input (variables × time).\n" *
          "For multi-node computation, provide htypes.\n" *
          "Got input shape: $(size(input))")
end

# 3D computation (multi-node, htypes = Vector{Int})
function (flux::HydroFlux{E,F,Vector{Int},I})(
    input::AbstractArray{T,3},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,3} where {E,F,I,T}
    params = _as_componentvector(params)
    expand_params = expand_component_params(params, get_param_names(flux), flux.htypes)
    output = flux.func(eachslice(input, dims=1), expand_params)
    stack(output, dims=1)
end

# Fallback: multi-node flux receiving 2D input → treat as single-node
function (flux::HydroFlux{E,F,Vector{Int},I})(
    input::AbstractArray{T,2},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,2} where {E,F,I,T}
    params = _as_componentvector(params)
    result = flux.func(eachslice(input, dims=1), params)
    stack(result, dims=1)
end

# ============================================================================
# StateFlux - State derivative component (unchanged)
# ============================================================================

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
(::StateFlux)(::AbstractArray, ::AbstractVector; kwargs...) =
    error("StateFlux cannot be run directly, please use HydroBucket or HydroRoute")

# ============================================================================
# Deprecated aliases for backward compatibility
# ============================================================================

const HydroMultiFlux = HydroFlux

macro hydromultiflux(args...)
    @warn "@hydromultiflux is deprecated, use @hydroflux with htypes instead" maxlog=1
    esc(:(@hydroflux $(args...)))
end

# Export interfaces
export HydroFlux, HydroMultiFlux, StateFlux
export @hydroflux, @hydromultiflux, @stateflux
