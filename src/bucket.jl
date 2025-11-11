"""
Bucket module - defines hydrological bucket components, supporting both symbolic and functional construction approaches.
"""

"""
    HydroBucket{S, MS, FF, OF, HT, NT} <: AbstractHydroBucket

A container component that organizes and executes a collection of hydrological flux components (fluxes and dfluxes).

It automatically compiles the components into efficient, callable functions for both single-node and multi-node simulations.

$(FIELDS)

# Type Parameters
- `S`: Whether contains state variables
- `MS`: Whether in multi-node mode
- `FF`: flux function type
- `OF`: ODE function type
- `HT`: HRU type vector type
- `NT`: Metadata type
"""
struct HydroBucket{S,MS,FF,OF,HT,NT} <: AbstractHydroBucket
    "bucket name"
    name::Symbol
    "Generated function for calculating all hydrological fluxes"
    flux_func::FF
    "Generated function for ODE calculations, or nothing if no ODE calculations needed"
    ode_func::OF
    "nodes type"
    hru_types::HT
    "Metadata about the bucket"
    infos::NT
    
    # Symbolic constructor
    function HydroBucket(;
        name::Optional{Symbol}=nothing,
        fluxes::Vector{<:AbstractFlux},
        dfluxes::Vector{<:AbstractStateFlux}=StateFlux[],
        hru_types::Vector{Int}=Int[]
    )
        inputs, outputs, states = get_var_names(vcat(fluxes, dfluxes))
        params = reduce(union, get_param_names.(vcat(fluxes, dfluxes)); init=Symbol[])
        nns = reduce(union, get_nn_names.(fluxes); init=Symbol[])
        
        infos = HydroInfos(
            inputs=inputs, states=states, outputs=outputs,
            params=params, nns=nns
        )
        
        flux_func, ode_func = build_bucket_func(fluxes, dfluxes, infos, length(hru_types) > 1)
        bucket_name = isnothing(name) ? Symbol("##bucket#", hash(infos)) : name
        hasstates, ismul = !isempty(states), length(hru_types) > 1
        
        return new{hasstates,ismul,typeof(flux_func),typeof(ode_func),typeof(hru_types),typeof(infos)}(
            bucket_name, flux_func, ode_func, hru_types, infos
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
        hru_types::Vector{Int}=Int[]
    )
        infos = HydroInfos(
            inputs=inputs, states=states, outputs=outputs,
            params=params, nns=Symbol[]
        )
        
        hasstates, ismul = !isempty(states), length(hru_types) > 1
        
        return new{hasstates,ismul,typeof(flux_func),typeof(ode_func),typeof(hru_types),typeof(infos)}(
            name, flux_func, ode_func, hru_types, infos
        )
    end
end

"""
    @hydrobucket [name] begin ... end

Macro to conveniently create a HydroBucket.

# Usage
The macro takes an optional name and a begin...end block containing the component definitions.

```jldoctest
julia> @hydrobucket :my_bucket begin
           fluxes = begin
               flux1
               flux2
           end
           dfluxes = begin
               state_flux1
           end
           hru_types = [1, 1, 2, 2]  # Optional
       end
```
"""
macro hydrobucket(args...)
    name = length(args) == 1 ? nothing : args[1]
    expr = length(args) == 1 ? args[1] : args[2]
    
    @assert Meta.isexpr(expr, :block) "Expected a begin...end block after bucket name"
    
    fluxes_expr, dfluxes_expr, hru_types_val = nothing, nothing, Int[]
    
    for assign in filter(x -> !(x isa LineNumberNode), expr.args)
        @assert Meta.isexpr(assign, :(=)) "Expected assignments in the form 'fluxes = begin...end'"
        lhs, rhs = assign.args
        
        if lhs == :fluxes
            @assert Meta.isexpr(rhs, :block) "Expected 'fluxes' to be defined in a begin...end block"
            fluxes_expr = Expr(:vect, filter(x -> !(x isa LineNumberNode), rhs.args)...)
        elseif lhs == :dfluxes
            @assert Meta.isexpr(rhs, :block) "Expected 'dfluxes' to be defined in a begin...end block"
            dfluxes_expr = Expr(:vect, filter(x -> !(x isa LineNumberNode), rhs.args)...)
        elseif lhs == :hru_types
            hru_types_val = rhs
        else
            error("Unknown assignment: $(lhs). Expected 'fluxes', 'dfluxes', or 'hru_types'")
        end
    end
    
    kwargs = [:($:(name=$name)), :($:(fluxes=$fluxes_expr)), :($:(hru_types=$hru_types_val))]
    if !isnothing(dfluxes_expr)
        push!(kwargs, :($:(dfluxes=$dfluxes_expr)))
    end
    
    return esc(:(HydroBucket(; $(kwargs...))))
end

"""
    (bucket::HydroBucket)(input, params, config; kwargs...)

Execute the HydroBucket model. This is the functor implementation for HydroBucket.

The method dispatches based on input data dimensions (2D for single-node, 3D for multi-node) 
and whether the bucket contains state variables. It solves the ODEs for state variables (if any) 
and then computes the output fluxes.

Common kwargs include `initstates`.
"""
# With states + single node
function (bucket::HydroBucket{true,false})(
    input::AbstractArray{T,2},
    params::ComponentVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,2} where {T}
    config_norm = normalize_config(config)
    solve_type = config_norm.solver
    interp_type = config_norm.interpolator
    
    param_vec, params_axes = Vector(params), getaxes(params)
    initstates = get(kwargs, :initstates, zeros(T, length(get_state_names(bucket))))
    timeidx = isempty(config_norm.timeidx) ? collect(1:size(input, 2)) : config_norm.timeidx
    # Create interpolator
    itpfuncs = hydrointerp(interp_type, input, timeidx)
    
    # Solve state variables
    solved_states = hydrosolve(
        solve_type,
        (u, p, t) -> bucket.ode_func(itpfuncs(t), u, ComponentVector(p, params_axes)),
        param_vec, Vector(initstates), timeidx, config_norm
    )
    
    # Compute flux outputs
    flux_output = bucket.flux_func(input, solved_states, params)
    flux_output_arr = stack(flux_output, dims=1)
    
    # Merge states and outputs
    vcat(solved_states, flux_output_arr)
end

# With states + multi-node
function (bucket::HydroBucket{true,true})(
    input::AbstractArray{T,3},
    params::ComponentVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,3} where {T}
    config_norm = normalize_config(config)
    solve_type = config_norm.solver
    interp_type = config_norm.interpolator
    device = config_norm.device
    
    input_dims, num_nodes, time_len = size(input)
    new_params = expand_component_params(params, get_param_names(bucket), bucket.hru_types)
    params_vec, params_axes = Vector(new_params) |> device, getaxes(new_params)
    
    num_states = length(get_state_names(bucket))
    initstates = get(kwargs, :initstates, zeros(T, num_states * num_nodes))
    timeidx = isempty(config_norm.timeidx) ? collect(1:time_len) : config_norm.timeidx
    
    # Create interpolator (flatten 3D input to 2D)
    itpfuncs = hydrointerp(interp_type, reshape(input, input_dims * num_nodes, :), timeidx)
    
    # Solve state variables
    solved_states = hydrosolve(
        solve_type,
        (u, p, t) -> bucket.ode_func(
            reshape(itpfuncs(t), input_dims, num_nodes),
            reshape(u, num_nodes, num_states)',
            ComponentVector(p, params_axes),
        ),
        params_vec, initstates, timeidx, config_norm
    )
    
    # Reshape state array to (num_states, num_nodes, time_len)
    solved_states_reshape = permutedims(
        reshape(solved_states, num_nodes, num_states, time_len),
        (2, 1, 3)
    )
    
    # Compute flux outputs
    output = bucket.flux_func(input, solved_states_reshape, new_params)
    
    # Merge states and outputs
    cat(solved_states_reshape, stack(output, dims=1), dims=1)
end

# Without states + single node
function (bucket::HydroBucket{false,false})(
    input::AbstractArray{T,2},
    params::ComponentVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,2} where {T}
    stack(bucket.flux_func(input, nothing, params), dims=1)
end

# Without states + multi-node
function (bucket::HydroBucket{false,true})(
    input::AbstractArray{T,3},
    params::ComponentVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,3} where {T}
    new_params = expand_component_params(params, get_param_names(bucket), bucket.hru_types)
    stack(bucket.flux_func(input, nothing, new_params), dims=1)
end

# Export interfaces
export HydroBucket, @hydrobucket
