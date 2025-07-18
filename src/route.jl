"""
    HydroRoute{N} <: AbstractHydroRoute

Represents a hydrological routing component using customizable flux equations.

# Fields
- `rfluxes::Vector{<:AbstractFlux}`: Flux functions defining routing processes (e.g., outflow calculation).
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
struct HydroRoute{N} <: AbstractHydroRoute
    "Routing function"
    rfluxes::Vector{<:AbstractFlux}
    "State derivative function"
    dfluxes::Vector{<:AbstractStateFlux}
    "Generated function for calculating all hydrological fluxes."
    flux_func::Function
    "Generated function for ordinary differential equations (ODE) calculations, or nothing if no ODE calculations are needed."
    ode_func::Function
    "Outflow projection function"
    aggr_func::Function
    "Metadata: contains keys for input, output, param, state, and nn"
    infos::NamedTuple

    function HydroRoute(;
        rfluxes::Vector{<:AbstractFlux},
        dfluxes::Vector{<:AbstractStateFlux},
        aggr_func::Function,
        name::Optional{Symbol}=nothing,
    )
        #* Extract all variable names of funcs and dfuncs
        inputs, outputs, states = get_vars(rfluxes, dfluxes)
        params = reduce(union, get_params.(rfluxes))
        nns = reduce(union, get_nns.(rfluxes))
        #* Setup the name information of the hydrobucket
        infos = (; inputs, outputs, states, params, nns)
        #* define the route name
        route_name = isnothing(name) ? Symbol("##route#", hash(infos)) : name
        #* build the route function
        flux_func, ode_func = build_route_func(rfluxes, dfluxes, infos)
        return new{route_name}(rfluxes, dfluxes, flux_func, ode_func, aggr_func, infos)
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
    @assert (!isnothing(fluxes_expr) && !isnothing(dfluxes_expr) && !isnothing(aggr_func_expr),
        "'fluxes', 'dfluxes', and 'aggr_func' must all be specified")
    return esc(quote
        fluxes = $fluxes_expr
        dfluxes = $dfluxes_expr
        aggr_func = $aggr_func_expr
        HydroRoute(
            rfluxes=fluxes,
            dfluxes=dfluxes,
            aggr_func=aggr_func,
            name=$(name)
        )
    end)
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
    kwargs...
) where {T<:Number}
    input_dims, num_nodes, time_len = size(input)

    #* get kwargs
    ptyidx = get(kwargs, :ptyidx, collect(1:num_nodes))
    styidx = get(kwargs, :styidx, collect(1:num_nodes))
    interp = get(kwargs, :interp, DirectInterpolation)
    solver = get(kwargs, :solver, ManualSolver(mutable=true))
    timeidx = get(kwargs, :timeidx, collect(1:time_len))
    device = get(kwargs, :device, solver.dev)

    #* prepare initstates
    initstates = get(kwargs, :initstates, zeros(eltype(params), length(get_state_names(route)), num_nodes)) |> device
    initstates_ = initstates isa ComponentVector ? initstates[get_state_names(route)] : initstates
    initstates_mat = expand_component_initstates(initstates_, styidx) |> device

    #* prepare states parameters and nns
    new_pas = expand_component_params(params, get_param_names(route), ptyidx) |> device
    params_vec, params_axes = Vector(new_pas) |> device, getaxes(new_pas)

    #* prepare input function
    itpfunc = interp(reshape(input, input_dims * num_nodes, time_len), timeidx)
    solved_states = solver(
        (u, p, t) -> begin
            tmp_outflow, tmp_states = route.ode_func(
                itpfunc(t), u,
                ComponentVector(p, params_axes)
            )
            tmp_inflow_arr = stack(route.aggr_func.(tmp_outflow), dims=1)
            tmp_output_state = stack(tmp_states, dims=1)
            tmp_output_state .+ tmp_inflow_arr
        end,
        params_vec, initstates_mat, timeidx
    )
    #* run flux functions
    output = route.flux_func(input, solved_states, new_pas)
    cat(solved_states, stack(output, dims=1), dims=1)
end