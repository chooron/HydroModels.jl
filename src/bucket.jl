"""
Bucket module - defines hydrological bucket components, supporting both symbolic and functional construction approaches.
"""

"""
    HydroBucket{S, FF, OF, HT, I} <: AbstractHydroBucket

A container component that organizes and executes a collection of hydrological flux components.

When `htypes` is `nothing`, this is a single-node component that accepts 2D input (variables × time).
When `htypes` is a `Vector{Int}`, this is a multi-node component that accepts 3D input (variables × nodes × time)
with parameter sharing via htypes indexing.

$(FIELDS)

# Type Parameters
- `S`: Whether contains state variables
- `FF`: flux function type
- `OF`: ODE function type
- `HT`: HRU type (`Nothing` for 2D, `Vector{Int}` for 3D)
- `I`: Metadata type
"""
struct HydroBucket{S,FF,OF,HT,I} <: AbstractHydroBucket
    "bucket name"
    name::Symbol
    "Generated function for calculating all hydrological fluxes"
    flux_func::FF
    "Generated function for ODE calculations, or nothing if no ODE calculations needed"
    ode_func::OF
    "HRU type indices for parameter sharing (Nothing = single-node 2D, Vector{Int} = multi-node 3D)"
    htypes::HT
    "Metadata about the bucket"
    infos::I

    # Symbolic constructor
    function HydroBucket(;
        name::Optional{Symbol}=nothing,
        fluxes::Vector{<:AbstractFlux},
        dfluxes::Vector{<:AbstractStateFlux}=StateFlux[],
        htypes::Optional{Vector{Int}}=nothing,
    )
        inputs, outputs, states = get_var_names(vcat(fluxes, dfluxes))
        params = reduce(union, get_param_names.(vcat(fluxes, dfluxes)); init=Symbol[])
        nns = reduce(union, get_nn_names.(fluxes); init=Symbol[])

        infos = HydroInfos(
            inputs=inputs, states=states, outputs=outputs,
            params=params, nns=nns
        )

        is_multi = !isnothing(htypes)
        flux_func, ode_func = build_bucket_func(fluxes, dfluxes, infos, is_multi)
        bucket_name = isnothing(name) ? Symbol("##bucket#", hash(infos)) : name
        hasstates = !isempty(states)

        return new{hasstates,typeof(flux_func),typeof(ode_func),typeof(htypes),typeof(infos)}(
            bucket_name, flux_func, ode_func, htypes, infos
        )
    end

    # Functional constructor
    function HydroBucket(
        flux_func::Function,
        ode_func::Union{Function,Nothing};
        name::Symbol,
        inputs::Vector{Symbol},
        outputs::Vector{Symbol},
        states::Vector{Symbol}=Symbol[],
        params::Vector{Symbol}=Symbol[],
        htypes::Optional{Vector{Int}}=nothing,
    )
        infos = HydroInfos(
            inputs=inputs, states=states, outputs=outputs,
            params=params, nns=Symbol[]
        )
        hasstates = !isempty(states)

        return new{hasstates,typeof(flux_func),typeof(ode_func),typeof(htypes),typeof(infos)}(
            name, flux_func, ode_func, htypes, infos
        )
    end
end

# ============================================================================
# @hydrobucket macro - unified for single-node and multi-node
# ============================================================================

"""
    @hydrobucket [name] begin ... end

Macro to conveniently create a HydroBucket component.

When `htypes` is specified in the begin block, creates a multi-node bucket.
Otherwise creates a single-node bucket.

# Examples
```julia
# Single-node bucket
@hydrobucket :my_bucket begin
    fluxes = begin
        flux1
        flux2
    end
    dfluxes = begin
        state_flux1
    end
end

# Multi-node bucket with parameter sharing
@hydrobucket :my_bucket begin
    fluxes = begin
        flux1
    end
    dfluxes = begin
        state_flux1
    end
    htypes = [1, 1, 2, 2, 3]
end
```
"""
macro hydrobucket(args...)
    name = length(args) == 1 ? nothing : args[1]
    expr = length(args) == 1 ? args[1] : args[2]

    @assert Meta.isexpr(expr, :block) "Expected a begin...end block after bucket name"

    fluxes_expr, dfluxes_expr, htypes_val = nothing, nothing, nothing

    for assign in filter(x -> !(x isa LineNumberNode), expr.args)
        @assert Meta.isexpr(assign, :(=)) "Expected assignments in the form 'fluxes = begin...end'"
        lhs, rhs = assign.args

        if lhs == :fluxes
            @assert Meta.isexpr(rhs, :block) "Expected 'fluxes' to be defined in a begin...end block"
            fluxes_expr = Expr(:vect, filter(x -> !(x isa LineNumberNode), rhs.args)...)
        elseif lhs == :dfluxes
            @assert Meta.isexpr(rhs, :block) "Expected 'dfluxes' to be defined in a begin...end block"
            dfluxes_expr = Expr(:vect, filter(x -> !(x isa LineNumberNode), rhs.args)...)
        elseif lhs == :htypes
            htypes_val = rhs
        else
            error("Unknown assignment: $(lhs). Expected 'fluxes', 'dfluxes', or 'htypes'")
        end
    end

    kwargs = [:($:(name=$name)), :($:(fluxes=$fluxes_expr))]
    if !isnothing(dfluxes_expr)
        push!(kwargs, :($:(dfluxes=$dfluxes_expr)))
    end
    if !isnothing(htypes_val)
        push!(kwargs, :($:(htypes=$htypes_val)))
    end

    return esc(:(HydroBucket(; $(kwargs...))))
end

# ============================================================================
# Functor methods - dispatch on S (has states) and HT (2D vs 3D)
# ============================================================================

# 2D with states (single-node)
function (bucket::HydroBucket{true,FF,OF,Nothing,I})(
    input::AbstractArray{T,2},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,2} where {FF,OF,I,T}
    params = _as_componentvector(params)
    config_norm = normalize_config(config)
    solve_type = config_norm.solver
    interp_type = config_norm.interpolator

    initstates = get(kwargs, :initstates, zeros(T, length(get_state_names(bucket))))
    timeidx = isempty(config_norm.timeidx) ? collect(1:size(input, 2)) : config_norm.timeidx
    itpfuncs = hydrointerp(interp_type, input, timeidx)

    solved_states = hydrosolve(
        solve_type,
        (u, p, t) -> bucket.ode_func(itpfuncs(t), u, p),
        params, Vector(initstates), timeidx, config_norm
    )

    flux_output = bucket.flux_func(input, solved_states, params)
    flux_output_arr = stack(flux_output, dims=1)
    vcat(solved_states, flux_output_arr)
end

# 2D without states (single-node)
function (bucket::HydroBucket{false,FF,OF,Nothing,I})(
    input::AbstractArray{T,2},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,2} where {FF,OF,I,T}
    params = _as_componentvector(params)
    stack(bucket.flux_func(input, nothing, params), dims=1)
end

# 3D with states (multi-node)
function (bucket::HydroBucket{true,FF,OF,Vector{Int},I})(
    input::AbstractArray{T,3},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,3} where {FF,OF,I,T}
    params = _as_componentvector(params)
    config_norm = normalize_config(config)
    solve_type = config_norm.solver
    interp_type = config_norm.interpolator

    input_dims, num_nodes, time_len = size(input)
    new_params = expand_component_params(params, get_param_names(bucket), bucket.htypes)

    num_states = length(get_state_names(bucket))
    initstates = get(kwargs, :initstates, zeros(T, num_states * num_nodes))
    timeidx = isempty(config_norm.timeidx) ? collect(1:time_len) : config_norm.timeidx

    itpfuncs = hydrointerp(interp_type, reshape(input, input_dims * num_nodes, :), timeidx)

    solved_states = hydrosolve(
        solve_type,
        (u, p, t) -> bucket.ode_func(
            reshape(itpfuncs(t), input_dims, num_nodes),
            reshape(u, num_nodes, num_states)',
            p,
        ),
        new_params, initstates, timeidx, config_norm
    )

    solved_states_reshape = permutedims(
        reshape(solved_states, num_nodes, num_states, time_len),
        (2, 1, 3)
    )

    output = bucket.flux_func(input, solved_states_reshape, new_params)
    cat(solved_states_reshape, stack(output, dims=1), dims=1)
end

# 3D without states (multi-node)
function (bucket::HydroBucket{false,FF,OF,Vector{Int},I})(
    input::AbstractArray{T,3},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,3} where {FF,OF,I,T}
    params = _as_componentvector(params)
    new_params = expand_component_params(params, get_param_names(bucket), bucket.htypes)
    stack(bucket.flux_func(input, nothing, new_params), dims=1)
end

# Error: single-node bucket receiving 3D input
function (bucket::HydroBucket{S,FF,OF,Nothing,I})(
    input::AbstractArray{T,3},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
) where {S,FF,OF,I,T}
    error("HydroBucket without htypes only accepts 2D input (variables × time).\n" *
          "For multi-node computation, provide htypes.\n" *
          "Got input shape: $(size(input))")
end

# Error: multi-node bucket receiving 2D input
function (bucket::HydroBucket{S,FF,OF,Vector{Int},I})(
    input::AbstractArray{T,2},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
) where {S,FF,OF,I,T}
    error("HydroBucket with htypes only accepts 3D input (variables × nodes × time).\n" *
          "For single-node computation, omit htypes.\n" *
          "Got input shape: $(size(input))")
end

# ============================================================================
# Deprecated aliases for backward compatibility
# ============================================================================

const HydroMultiBucket = HydroBucket

macro hydromultibucket(args...)
    @warn "@hydromultibucket is deprecated, use @hydrobucket with htypes instead" maxlog=1
    esc(:(@hydrobucket $(args...)))
end

# Export interfaces
export HydroBucket, HydroMultiBucket, @hydrobucket, @hydromultibucket
