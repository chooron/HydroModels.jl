"""
    HydroRoute{N} <: AbstractHydroRoute

Represents a hydrological routing component using customizable flux equations.

# Fields
- `rfluxes::Vector{<:AbstractHydroFlux}`: Flux functions defining routing processes (e.g., outflow calculation).
- `multi_flux_func::Function`: Compiled function to evaluate all non-state-derivative fluxes.
- `multi_ode_func::Function`: Compiled function for state derivatives used by ODE solvers.
- `aggr_func::Function`: Function to aggregate or distribute flows between connected elements (often identity or matrix multiplication).
- `infos::NamedTuple`: Metadata (inputs, outputs, states, params, nns).

# Constructor
```julia
HydroRoute(; rfluxes, dfluxes, aggr_func, name=nothing)
```
- `rfluxes`: Vector of routing fluxes (e.g., `HydroFlux`, `NeuralFlux`).
- `dfluxes`: Vector of state derivative fluxes (`StateFlux`).
- `aggr_func`: Aggregation/distribution function.
- `name`: Optional symbol identifier.

# Notes
- Builds specialized functions (`multi_flux_func`, `multi_ode_func`) for efficient computation.
- Designed for use within a larger model (e.g., `HydroModel`).
- Handles variable management and metadata tracking internally.
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
        inputs, outputs, states = get_var_names(vcat(fluxes, dfluxes))
        params = reduce(union, get_param_names.(vcat(fluxes, dfluxes)))
        nns = reduce(union, get_nn_names.(fluxes))
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
    @hydroroute name begin ... end

Macro to simplify the construction of a `HydroRoute` object.

# Usage
```julia
@hydroroute :my_route begin
    fluxes = [
        @hydroflux outflow ~ k * storage
    ]
    dfluxes = [
        @stateflux storage ~ inflow - outflow
    ]
    aggr_func = identity # or other function
end
```
Defines a `HydroRoute` with the specified name (`:my_route`), fluxes, state derivatives (`dfluxes`), and aggregation function (`aggr_func`) within the `begin...end` block.

# Arguments
- `name`: (Optional) Symbol for the route name.
- `block`: A `begin...end` block containing assignments for `fluxes`, `dfluxes`, and `aggr_func`.

# Returns
- A `HydroRoute` instance.
"""
macro hydroroute(args...)
    name = length(args) == 1 ? nothing : args[1]
    expr = length(args) == 1 ? args[1] : args[2]

    @assert Meta.isexpr(expr, :block) "Expected a begin...end block after route name"
    fluxes_expr, dfluxes_expr, aggr_func_expr = nothing, nothing, nothing
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
        else
            error("Unknown assignment: $lhs. Expected 'fluxes', 'dfluxes', or 'aggr_func'")
        end
    end
    err_msg = "'fluxes', 'dfluxes', and 'aggr_func' must all be specified"
    @assert !isnothing(fluxes_expr) && !isnothing(dfluxes_expr) && !isnothing(aggr_func_expr) err_msg
    return esc(quote
        HydroRoute(
            rfluxes=$fluxes_expr,
            dfluxes=$dfluxes_expr,
            aggr_func=$aggr_func_expr,
            name=$(name)
        )
    end)
end

"""
    Build the aggregation function from a given graph
"""
function build_aggr_func(graph::DiGraph)
    adjacency = adjacency_matrix(graph)'
    aggr_func = (outflow) -> adjacency * outflow
    return aggr_func
end

"""
    Build the aggregation function from a given flow direction matrix and positions
"""
function build_aggr_func(flwdir::AbstractMatrix, positions::AbstractVector)
    d8_codes = [1, 2, 4, 8, 16, 32, 64, 128]
    d8_nn_pads = [(1, 1, 2, 0), (2, 0, 2, 0), (2, 0, 1, 1), (2, 0, 0, 2), (1, 1, 0, 2), (0, 2, 0, 2), (0, 2, 1, 1), (0, 2, 2, 0)]

    # input dims: node_num * ts_len
    function grid_routing(input::AbstractVector, positions::AbstractVector, flwdir::AbstractMatrix)
        # Convert input to sparse matrix
        input_arr = Array(sparse([pos[1] for pos in positions], [pos[2] for pos in positions], input, size(flwdir)[1], size(flwdir)[2]))
        # Calculate weighted summation
        input_routed = sum(collect([pad_zeros(input_arr .* (flwdir .== code), arg) for (code, arg) in zip(d8_codes, d8_nn_pads)]))
        # Clip input matrix border
        clip_arr = input_routed[2:size(input_arr)[1]+1, 2:size(input_arr)[2]+1]
        # Convert input matrix to vector
        collect([clip_arr[pos[1], pos[2]] for pos in positions])
    end
    #* build the outflow projection function
    aggr_func = (outflow) -> grid_routing(outflow, positions, flwdir)
    return aggr_func
end

""" Helper to build multi-dimensional bucket and routing functions. """
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
        $(generate_compute_calls(fluxes, dims=2)...)
        return [$(outputs...), vcat($(map(expr -> :(@. $(simplify_expr(toexpr(expr)))), reduce(vcat, get_exprs.(dfluxes)))...))]
    end)
    generated_multi_flux_func = @RuntimeGeneratedFunction(multi_flux_func_expr)

    if !isempty(names.state_names)
        multi_diff_func_expr = :(function (inputs, states, pas)
            Base.@_inline_meta
            $(define_calls_1...)
            $(generate_compute_calls(fluxes, dims=1)...)
            $(generate_states_expression(dfluxes=dfluxes, broadcast=true))
        end)
        generated_multi_diff_func = @RuntimeGeneratedFunction(multi_diff_func_expr)
        return generated_multi_flux_func, generated_multi_diff_func
    else
        return generated_multi_flux_func, (_) -> nothing
    end
end


"""
    (route::HydroRoute)(input::AbstractArray{T,3}, params::ComponentVector; kwargs...)

Runs the simulation for the `HydroRoute` component.

# Arguments
- `input::AbstractArray{T,3}`: Input data array with shape (variables, nodes, timesteps).
- `params::ComponentVector`: Parameters for the routing fluxes.

# Keyword Arguments
- `initstates`: Initial states for the component. Defaults to zeros.
- `ptyidx`: Indices mapping parameters to nodes. Defaults to `1:num_nodes`.
- `styidx`: Indices mapping states to nodes. Defaults to `1:num_nodes`.
- `interp`: Interpolation method for time-varying inputs (e.g., `DirectInterpolation`).
- `solver`: ODE solver instance (e.g., `ManualSolver`).
- `timeidx`: Time indices for the simulation. Defaults to `1:time_len`.
- `device`: Computing device (e.g., `identity` for CPU).

# Returns
- `AbstractArray`: Output array containing state evolution and calculated fluxes, shape (output_variables, nodes, timesteps).

# Example
```julia
# Assuming 'route' is a HydroRoute instance
output = route(input_data, parameters; solver=Tsit5())
```

# Notes
- Expects 3D input (variables, nodes, timesteps).
- Integrates state variables over time using the specified solver.
- Applies the `multi_flux_func` to calculate outputs based on inputs and states.
"""
function (route::HydroRoute)(
    input::AbstractArray{T,3},
    params::ComponentVector;
    solver::S=ManualSolver(mutable=true),
    interp::I=DirectInterpolation,
    timeidx::AbstractVector=collect(1:size(input, 3)),
    initstates::AbstractArray=zeros(eltype(params), length(get_state_names(route)), size(input, 2)),
    device=solver.dev,
    kwargs...
) where {T,S,I}
    input_dims, num_nodes, time_len = size(input)
    initstates_expand = expand_component_initstates(initstates, styidx) |> device
    new_pas = expand_component_params(params, get_param_names(route), ptyidx) |> device
    params_vec, params_axes = Vector(new_pas) |> device, getaxes(new_pas)

    itpfunc = interp(reshape(input, input_dims * num_nodes, time_len), timeidx)
    solved_states = solver(
        (u, p, t) -> begin
            tmp_outflow, tmp_states = route.ode_func(
                reshape(itpfunc(t), input_dims, num_nodes), u,
                ComponentVector(p, params_axes)
            )
            tmp_inflow_arr = route.aggr_func(tmp_outflow)
            tmp_states .+ tmp_inflow_arr
        end,
        params_vec, initstates_expand, timeidx
    )

    solved_states_reshape = reshape(solved_states, length(get_state_names(route)), num_nodes, time_len)
    output = route.flux_func(input, solved_states_reshape, new_pas)
    cat(solved_states_reshape, stack(output, dims=1), dims=1)
end