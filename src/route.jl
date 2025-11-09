"""
Route module - defines hydrological routing components, supporting both symbolic and functional construction approaches.
"""

"""
    HydroRoute{FF, OF, AF, HT, NT} <: AbstractHydroRoute

Represents a hydrological routing component that combines flux calculations with an aggregation/distribution step.

It is designed for routing flows between connected nodes or grid cells in a spatial model. 
It compiles flux equations and an aggregation function into efficient callable functions.

$(FIELDS)

# Type Parameters
- `FF`: flux function type
- `OF`: ODE function type
- `AF`: aggregation function type
- `HT`: HRU type vector type
- `NT`: Metadata type
"""
struct HydroRoute{FF,OF,AF,HT,NT} <: AbstractHydroRoute
    "route name"
    name::Symbol
    "Generated function for calculating all hydrological fluxes"
    flux_func::FF
    "Generated function for ODE calculations, or nothing if no ODE calculations needed"
    ode_func::OF
    "Outflow projection/aggregation function"
    aggr_func::AF
    "nodes type"
    hru_types::HT
    "Metadata: contains keys for input, output, param, state, and nn"
    infos::NT
    
    # Symbolic constructor
    function HydroRoute(;
        rfluxes::Vector{<:AbstractHydroFlux},
        dfluxes::Vector{<:AbstractStateFlux},
        hru_types::Vector{Int},
        aggr_func::AF,
        name::Optional{Symbol}=nothing,
    ) where AF
        inputs, outputs, states = get_var_names(vcat(rfluxes, dfluxes))
        params = reduce(union, get_param_names.(vcat(rfluxes, dfluxes)); init=Symbol[])
        nns = reduce(union, get_nn_names.(rfluxes); init=Symbol[])
        
        infos = HydroInfos(
            inputs=inputs, states=states, outputs=outputs,
            params=params, nns=nns
        )
        route_name = isnothing(name) ? Symbol("##route#", hash(infos)) : name
        flux_func, ode_func = build_route_func(rfluxes, dfluxes, infos)
        
        return new{typeof(flux_func),typeof(ode_func),AF,typeof(hru_types),typeof(infos)}(
            route_name, flux_func, ode_func, aggr_func, hru_types, infos
        )
    end
    
    # Functional constructor
    function HydroRoute(
        flux_func::Function,
        ode_func::Function,
        aggr_func::AF;
        name::Symbol,
        inputs::Vector{Symbol},
        outputs::Vector{Symbol},
        states::Vector{Symbol},
        params::Vector{Symbol}=Symbol[],
        hru_types::Vector{Int}
    ) where AF
        infos = HydroInfos(
            inputs=inputs, states=states, outputs=outputs,
            params=params, nns=Symbol[]
        )
        
        return new{typeof(flux_func),typeof(ode_func),AF,typeof(hru_types),typeof(infos)}(
            name, flux_func, ode_func, aggr_func, hru_types, infos
        )
    end
end

"""
    @hydroroute [name] begin ... end

Macro to simplify the construction of a HydroRoute component.

# Usage
The macro takes an optional name and a begin...end block containing definitions for fluxes, 
state derivatives, node types, and the aggregation function.

```jldoctest
julia> @hydroroute :my_route begin
           fluxes = begin
               @hydroflux outflow ~ k * storage
           end
           dfluxes = begin
               @stateflux storage ~ inflow - outflow
           end
           hru_types = [1, 2, 3]
           aggr_func = identity  # or a matrix, or a custom function
       end
```
"""
macro hydroroute(args...)
    name = length(args) == 1 ? nothing : args[1]
    expr = length(args) == 1 ? args[1] : args[2]
    
    @assert Meta.isexpr(expr, :block) "Expected a begin...end block after route name"
    
    fluxes_expr, dfluxes_expr, aggr_func_expr, hru_types_expr = nothing, nothing, nothing, nothing
    
    for assign in filter(x -> !(x isa LineNumberNode), expr.args)
        @assert Meta.isexpr(assign, :(=)) "Expected assignments in the form 'fluxes = begin...end', etc."
        lhs, rhs = assign.args
        
        if lhs == :fluxes
            @assert Meta.isexpr(rhs, :block) "Expected 'fluxes' to be defined in a begin...end block"
            fluxes_expr = Expr(:vect, filter(x -> !(x isa LineNumberNode), rhs.args)...)
        elseif lhs == :dfluxes
            @assert Meta.isexpr(rhs, :block) "Expected 'dfluxes' to be defined in a begin...end block"
            dfluxes_expr = Expr(:vect, filter(x -> !(x isa LineNumberNode), rhs.args)...)
        elseif lhs == :aggr_func
            aggr_func_expr = rhs
        elseif lhs == :hru_types
            hru_types_expr = rhs
        else
            error("Unknown assignment: $(lhs). Expected 'fluxes', 'dfluxes', 'hru_types', or 'aggr_func'")
        end
    end
    
    err_msg = "'fluxes', 'dfluxes', 'hru_types', and 'aggr_func' must all be specified"
    @assert !isnothing(fluxes_expr) && !isnothing(dfluxes_expr) && 
            !isnothing(aggr_func_expr) && !isnothing(hru_types_expr) err_msg
    
    return esc(quote
        HydroRoute(
            rfluxes=$fluxes_expr,
            dfluxes=$dfluxes_expr,
            aggr_func=$aggr_func_expr,
            hru_types=$hru_types_expr,
            name=$(name)
        )
    end)
end

"""
    build_aggr_func(graph::DiGraph)

Build an aggregation function from a directed graph. 

The function performs routing by multiplying the outflow vector with the graph's adjacency matrix.

# Arguments
- `graph`: Directed graph representing connections between nodes

# Returns
- Aggregation function that accepts an outflow vector and returns a routed inflow vector
"""
function build_aggr_func(graph::DiGraph)
    adjacency = adjacency_matrix(graph)'
    aggr_func = (outflow) -> adjacency * outflow
    return aggr_func
end

"""
    build_aggr_func(flwdir::AbstractMatrix, positions::AbstractVector)

Build a grid-based routing aggregation function from a flow direction matrix (e.g., D8) 
and a vector of grid cell positions.

# Arguments
- `flwdir`: Flow direction matrix (D8 encoding)
- `positions`: Vector of grid cell positions

# Returns
- Aggregation function that accepts an outflow vector and returns a routed inflow vector

# Note
This implementation avoids in-place modifications and is fully compatible with Zygote automatic differentiation.
"""
function build_aggr_func(flwdir::AbstractMatrix, positions::AbstractVector)
    d8_codes = [1, 2, 4, 8, 16, 32, 64, 128]
    d8_nn_pads = [
        (1, 1, 2, 0), (2, 0, 2, 0), (2, 0, 1, 1), (2, 0, 0, 2),
        (1, 1, 0, 2), (0, 2, 0, 2), (0, 2, 1, 1), (0, 2, 2, 0)
    ]
    
    # Precompute CartesianIndex for performance
    cartesian_positions = [CartesianIndex(p...) for p in positions]
    
    function grid_routing(input::AbstractVector)
        rows, cols = size(flwdir)
        T = eltype(input)
        
        # Create input array (avoiding in-place modification)
        input_arr = zeros(T, rows, cols)
        input_arr[cartesian_positions] = input
        
        # Compute routing for each D8 direction (using pure functional operations)
        input_routed = sum(d8_codes, d8_nn_pads) do code, pad
            pad_zeros(input_arr .* (flwdir .== code), pad)
        end
        
        # Clip borders and extract results
        clip_arr = input_routed[2:rows+1, 2:cols+1]
        return clip_arr[cartesian_positions]
    end
    
    return grid_routing
end

"""
    (route::HydroRoute)(input, params, config; kwargs...)

Run the simulation for the HydroRoute component. This is the functor implementation.

It expects 3D input `(variables, nodes, timesteps)`. The method integrates the state variables 
over time using a specified ODE solver, applying the `aggr_func` at each step to handle flow 
between nodes. It then computes the final output fluxes.

Common kwargs include `initstates`.
"""
function (route::HydroRoute)(
    input::AbstractArray{T,3},
    params::ComponentVector,
    config::ConfigType=default_config();
    kwargs...
) where {T}
    config_norm = normalize_config(config)
    solve_type = config_norm.solver
    interp_type = config_norm.interpolator
    device = config_norm.device
    
    input_dims, num_nodes, time_len = size(input)
    new_pas = expand_component_params(params, get_param_names(route), route.hru_types) |> device
    params_vec, params_axes = Vector(new_pas) |> device, getaxes(new_pas)
    
    num_states = length(get_state_names(route))
    initstates = get(kwargs, :initstates, zeros(T, num_states * num_nodes))
    timeidx = isempty(config_norm.timeidx) ? collect(1:time_len) : config_norm.timeidx
    
    # Create interpolator
    itpfunc = hydrointerp(interp_type, reshape(input, input_dims * num_nodes, time_len), timeidx)
    
    # Solve state variables (avoiding in-place modifications)
    solved_states = hydrosolve(
        solve_type,
        (u, p, t) -> begin
            # Compute current time outflow and state changes
            tmp_outflow, tmp_states = route.ode_func(
                reshape(itpfunc(t), input_dims, num_nodes),
                u,
                ComponentVector(p, params_axes)
            )
            # Apply aggregation function to get inflow
            tmp_inflow_arr = route.aggr_func(tmp_outflow)
            # Return total state change (avoiding in-place .+ operation)
            tmp_states .+ tmp_inflow_arr
        end,
        params_vec, initstates, timeidx, config_norm
    )
    
    # Reshape state array
    solved_states_reshape = reshape(solved_states, num_states, num_nodes, time_len)
    
    # Compute flux outputs
    output = route.flux_func(input, solved_states_reshape, new_pas)
    
    # Merge states and outputs
    cat(solved_states_reshape, stack(output, dims=1), dims=1)
end

# Export interfaces
export HydroRoute, @hydroroute, build_aggr_func
