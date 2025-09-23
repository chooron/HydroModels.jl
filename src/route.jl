"""
    HydroRoute{FF, OF, AF, HT, NT} <: AbstractHydroRoute

Represents a hydrological routing component that combines flux calculations with an aggregation/distribution step.

It is designed for routing flows between connected nodes or grid cells in a spatial model. It compiles flux equations and an aggregation function into efficient callable functions.

$(FIELDS)
"""
struct HydroRoute{FF,OF,AF,HT,NT} <: AbstractHydroRoute
    "route names"
    name::Symbol
    "Generated function for calculating all hydrological fluxes."
    flux_func::FF
    "Generated function for ordinary differential equations (ODE) calculations, or nothing if no ODE calculations are needed."
    ode_func::OF
    "Outflow projection function"
    aggr_func::AF
    "nodes type"
    hru_types::HT
    "Metadata: contains keys for input, output, param, state, and nn"
    infos::NT

    function HydroRoute(;
        rfluxes::Vector{<:AbstractHydroFlux},
        dfluxes::Vector{<:AbstractStateFlux},
        hru_types::Vector{Int},
        aggr_func::AF,
        name::Optional{Symbol}=nothing,
    ) where AF
        inputs, outputs, states = get_var_names(vcat(rfluxes, dfluxes))
        params = reduce(union, get_param_names.(vcat(rfluxes, dfluxes)))
        nns = reduce(union, get_nn_names.(rfluxes))
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
end

"""
    @hydroroute [name] begin ... end

A macro to simplify the construction of a `HydroRoute` component.

# Usage
The macro takes an optional name and a `begin...end` block containing the definitions for fluxes, state derivatives, node types, and the aggregation function.

```julia
@hydroroute :my_route begin
    fluxes = begin
        @hydroflux outflow ~ k * storage
    end
    dfluxes = begin
        @stateflux storage ~ inflow - outflow
    end
    hru_types = [1, 2, 3]
    aggr_func = identity # or a matrix, or a custom function
end
```
"""
macro hydroroute(args...)
    name = length(args) == 1 ? nothing : args[1]
    expr = length(args) == 1 ? args[1] : args[2]

    @assert Meta.isexpr(expr, :block) "Expected a begin...end block after route name"
    fluxes_expr, dfluxes_expr, aggr_func_expr, hru_types_expr = nothing, nothing, nothing, nothing
    for assign in filter(x -> !(x isa LineNumberNode), expr.args)
        @assert Meta.isexpr(assign, :(=)) "Expected assignments in the form 'fluxes = begin...end', 'dfluxes = begin...end', and 'aggr_func = f(x)'"
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
            error("Unknown assignment: $lhs. Expected 'fluxes', 'dfluxes', 'hru_types', or 'aggr_func'")
        end
    end
    err_msg = "'fluxes', 'dfluxes', 'hru_types', and 'aggr_func' must all be specified"
    @assert !isnothing(fluxes_expr) && !isnothing(dfluxes_expr) && !isnothing(aggr_func_expr) && !isnothing(hru_types_expr) err_msg
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
$(TYPEDSIGNATURES)

Builds an aggregation function from a `DiGraph`. The function performs routing by multiplying the outflow vector with the graph's adjacency matrix.
"""
function build_aggr_func(graph::DiGraph)
    adjacency = adjacency_matrix(graph)'
    aggr_func = (outflow) -> adjacency * outflow
    return aggr_func
end

"""
$(TYPEDSIGNATURES)

Builds an aggregation function for grid-based routing from a flow direction matrix (e.g., D8) and a vector of grid cell positions.
"""
# function build_aggr_func(flwdir::AbstractMatrix, positions::AbstractVector)
#     d8_codes = [1, 2, 4, 8, 16, 32, 64, 128]
#     d8_nn_pads = [(1, 1, 2, 0), (2, 0, 2, 0), (2, 0, 1, 1), (2, 0, 0, 2), (1, 1, 0, 2), (0, 2, 0, 2), (0, 2, 1, 1), (0, 2, 2, 0)]

#     # input dims: node_num * ts_len
#     function grid_routing(input::AbstractVector, positions::AbstractVector, flwdir::AbstractMatrix)
#         # Convert input to sparse matrix
#         input_arr = Array(sparse([pos[1] for pos in positions], [pos[2] for pos in positions], input, size(flwdir)[1], size(flwdir)[2]))
#         # Calculate weighted summation
#         input_routed = sum(collect([pad_zeros(input_arr .* (flwdir .== code), arg) for (code, arg) in zip(d8_codes, d8_nn_pads)]))
#         # Clip input matrix border
#         clip_arr = input_routed[2:size(input_arr)[1]+1, 2:size(input_arr)[2]+1]
#         # Convert input matrix to vector
#         collect([clip_arr[pos[1], pos[2]] for pos in positions])
#     end
#     #* build the outflow projection function
#     aggr_func = (outflow) -> grid_routing(outflow, positions, flwdir)
#     return aggr_func
# end
function build_aggr_func(flwdir::AbstractMatrix, positions::AbstractVector)
    d8_codes = [1, 2, 4, 8, 16, 32, 64, 128]
    d8_nn_pads = [
        (1, 1, 2, 0), (2, 0, 2, 0), (2, 0, 1, 1), (2, 0, 0, 2),
        (1, 1, 0, 2), (0, 2, 0, 2), (0, 2, 1, 1), (0, 2, 2, 0)
    ]
    function grid_routing(input::AbstractVector, positions::AbstractVector, flwdir::AbstractMatrix)
        rows, cols = size(flwdir)
        input_arr = zeros(eltype(input), rows, cols)
        cartesian_positions = (p -> CartesianIndex(p...)).(positions)
        input_arr[cartesian_positions] = input
        input_routed = sum(pad_zeros(input_arr .* (flwdir .== code), arg) for (code, arg) in zip(d8_codes, d8_nn_pads))
        clip_arr = input_routed[2:rows+1, 2:cols+1]
        return clip_arr[cartesian_positions]
    end

    aggr_func = (outflow) -> grid_routing(outflow, positions, flwdir)
    return aggr_func
end

function build_route_func(fluxes::Vector{<:AbstractHydroFlux}, dfluxes::Vector{<:AbstractStateFlux}, infos::HydroModelCore.HydroInfos)
    nn_fluxes = filter(f -> f isa AbstractNeuralFlux, fluxes)

    define_calls_1 = [
        generate_var_assignments(vars=infos.inputs, target=:inputs, dims=1)...,
        generate_var_assignments(vars=infos.states, target=:states, dims=1)...,
        generate_param_assignments(params=infos.params)...,
        generate_nn_assignments(nnfluxes=nn_fluxes)...
    ]
    define_calls_2 = [
        generate_var_assignments(vars=infos.inputs, target=:inputs, dims=2)...,
        generate_var_assignments(vars=infos.states, target=:states, dims=2)...,
        generate_param_assignments(params=infos.params)...,
        generate_nn_assignments(nnfluxes=nn_fluxes)...
    ]

    multi_flux_func_expr = :(function (inputs, states, pas)
        Base.@_inline_meta
        $(define_calls_2...)
        $(generate_compute_calls(fluxes=fluxes, dims=1)...)
        return [$(infos.outputs...)]
    end)
    generated_multi_flux_func = @RuntimeGeneratedFunction(multi_flux_func_expr)

    if !isempty(infos.states)
        multi_diff_func_expr = :(function (inputs, states, pas)
            Base.@_inline_meta
            $(define_calls_1...)
            $(generate_compute_calls(fluxes=fluxes, dims=1)...)
            [$(infos.outputs...), vcat($(map(expr -> :(@. $(simplify_expr(toexpr(expr)))), reduce(vcat, get_exprs.(dfluxes)))...))]
        end)
        generated_multi_diff_func = @RuntimeGeneratedFunction(multi_diff_func_expr)
        return generated_multi_flux_func, generated_multi_diff_func
    else
        return generated_multi_flux_func, (_) -> nothing
    end
end


"""
    (route::HydroRoute)(input, params; kwargs...)

Runs the simulation for the `HydroRoute` component. This is the functor implementation.

It expects 3D input `(variables, nodes, timesteps)`. The method integrates the state variables over time using a specified ODE solver, applying the `aggr_func` at each step to handle flow between nodes. It then computes the final output fluxes.

Common `kwargs` include `solver`, `interp`, `timeidx`, and `initstates`.
"""
function (route::HydroRoute)(
    input::AbstractArray{T,3},
    params::ComponentVector,
    config::NamedTuple=DEFAULT_CONFIG;
    kwargs...
) where {T}
    solve_type = get(config, :solver, MutableSolver)
    interp_type = get(config, :interpolator, Val(DirectInterpolation))
    device = get(config, :device, identity)

    input_dims, num_nodes, time_len = size(input)
    new_pas = expand_component_params(params, get_param_names(route), route.hru_types) |> device
    params_vec, params_axes = Vector(new_pas) |> device, getaxes(new_pas)
    initstates = get(kwargs, :initstates, zeros(eltype(params), length(get_state_names(route)) * size(input, 2)))
    timeidx = get(kwargs, :timeidx, collect(1:time_len))

    itpfunc = hydrointerp(interp_type, reshape(input, input_dims * num_nodes, time_len), timeidx)
    solved_states = hydrosolve(solve_type,
        (u, p, t) -> begin
            tmp_outflow, tmp_states = route.ode_func(
                reshape(itpfunc(t), input_dims, num_nodes), u,
                ComponentVector(p, params_axes)
            )
            tmp_inflow_arr = route.aggr_func(tmp_outflow)
            tmp_states .+ tmp_inflow_arr
        end,
        params_vec, initstates, timeidx, config
    )

    solved_states_reshape = reshape(solved_states, length(get_state_names(route)), num_nodes, time_len)
    output = route.flux_func(input, solved_states_reshape, new_pas)
    cat(solved_states_reshape, stack(output, dims=1), dims=1)
end