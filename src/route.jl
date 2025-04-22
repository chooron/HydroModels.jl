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
    "Generated function for calculating all hydrological fluxes."
    multi_flux_func::Function
    "Generated function for ordinary differential equations (ODE) calculations, or nothing if no ODE calculations are needed."
    multi_ode_func::Function
    "Outflow projection function"
    aggr_func::Function
    "Metadata: contains keys for input, output, param, state, and nn"
    infos::NamedTuple

    function HydroRoute(;
        rfluxes::Vector{<:AbstractFlux},
        dfluxes::Vector{<:AbstractStateFlux},
        aggr_func::Function,
        name::Union{Symbol,Nothing}=nothing,
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
        multi_flux_func, multi_ode_func = build_route_func(rfluxes, dfluxes, infos)
        return new{route_name}(rfluxes, multi_flux_func, multi_ode_func, aggr_func, infos)
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
    @assert !isnothing(fluxes_expr) && !isnothing(dfluxes_expr) && !isnothing(aggr_func_expr) "'fluxes', 'dfluxes', and 'aggr_func' must all be specified"
    return esc(quote
        let
            fluxes = $fluxes_expr
            dfluxes = $dfluxes_expr
            aggr_func = $aggr_func_expr
            HydroRoute(
                rfluxes=fluxes,
                dfluxes=dfluxes,
                aggr_func=aggr_func,
                name=$(name)
            )
        end
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
    new_pas = expand_component_params(params, ptyidx) |> device
    params_vec, params_axes = Vector(new_pas) |> device, getaxes(new_pas)

    #* prepare input function
    itpfuncs = interp.(eachslice(input, dims=1), Ref(timeidx))
    solved_states = solver(
        (u, p, t) -> begin
            tmp_states, tmp_outflow = route.multi_ode_func(
                ntuple(i -> itpfuncs[i](t), input_dims),
                eachslice(u, dims=1),
                ComponentVector(p, params_axes)
            )
            tmp_states_arr = reduce(hcat, tmp_states)
            tmp_inflow_arr = reduce(hcat, route.aggr_func.(tmp_outflow))
            # todo 这里元编程表达一直存在问题
            tmp_states_arr .+ tmp_inflow_arr |> permutedims
        end,
        params_vec, initstates_mat, timeidx
    )
    #* run other functions
    output = route.multi_flux_func(eachslice(input, dims=1), eachslice(solved_states, dims=1), new_pas)
    cat(solved_states, stack(output, dims=1), dims=1)
end

"""
    RapidRoute{N} <: AbstractRoute

Implements the RAPID (Routing Application for Parallel computatIon of Discharge) routing model.

# Fields
- `adjacency::AbstractMatrix`: Sparse matrix representing the river network connectivity (upstream to downstream).
- `infos::NamedTuple`: Metadata including standard variable names (`inputs`, `outputs`, `states`, `params`). States are typically `[:discharge]`. Params are `[:rapid_k, :rapid_x]`.

# Constructor
```julia
RapidRoute(; network::AbstractGraph, name::Union{Symbol, Nothing}=nothing)
```
- `network`: A `Graphs.jl` compatible directed graph representing the river network.
- `name`: Optional symbol identifier.

# Notes
- Uses the Muskingum-Cunge method formulated for efficient matrix operations.
- Suitable for large-scale river networks.
- Parameters `k` (wave celerity) and `x` (diffusion) are expected per node via the `params` argument during simulation call.
"""
struct RapidRoute{N} <: AbstractRoute
    "Routing adjacency matrix"
    adjacency::AbstractMatrix
    "Metadata: contains keys for input, output, param, state, and nn"
    infos::NamedTuple

    function RapidRoute(
        fluxes::Pair{Vector{Num},Vector{Num}};
        network::DiGraph,
        name::Union{Symbol,Nothing}=nothing,
    )
        @parameters rapid_k rapid_x
        #* Extract all variable names of funcs and dfuncs
        inputs, outputs = fluxes[1], fluxes[2]
        @assert length(inputs) == length(outputs) == 1 "The length of inputs and outputs must be the 1, but got inputs: $(length(inputs)) and outputs: $(length(outputs))"
        #* Setup the name information of the hydrobucket
        @parameters rapid_k rapid_x
        infos = (; inputs=inputs, outputs=outputs, states=Num[], params=[rapid_k, rapid_x])
        route_name = isnothing(name) ? Symbol("##route#", hash(infos)) : name
        #* generate adjacency matrix from network
        adjacency = adjacency_matrix(network)'
        return new{route_name}(adjacency, infos)
    end
end

"""
    (route::RapidRoute)(input::Array, params::ComponentVector; kwargs...)

Executes the RAPID routing simulation.

# Arguments
- `input::Array`: Lateral inflow data, typically shape (1, nodes, timesteps) where the variable is inflow.
- `params::ComponentVector`: Parameters containing `:rapid_k` and `:rapid_x` for each node.

# Keyword Arguments
- `delta_t::Float64`: Simulation time step in seconds (default: 1.0).
- `device`: Computing device (e.g., `identity` for CPU, `gpu` function). Default `identity`.
- `ptyidx`: Indices mapping parameters to nodes. Default `1:num_nodes`.
- `interp`: Interpolation method for time-varying inputs. Default `DirectInterpolation`.
- `solver`: ODE solver instance. Default `ManualSolver(mutable=true)`.
- `timeidx`: Time indices for the simulation. Default `1:time_len`.
- `initstates`: Initial discharge values. Default zeros.

# Returns
- `AbstractArray{T,3}`: Output array of discharge, shape (1, nodes, timesteps).

# Example
```julia
# Assuming 'rapid_route' is a RapidRoute instance
discharge = rapid_route(lateral_inflow, parameters; delta_t=3600.0)
```

# Notes
- Solves the Muskingum-Cunge equations using a matrix-based approach.
- Calculates Muskingum coefficients (c0, c1, c2) internally based on `k`, `x`, and `delta_t`.
- Requires `input` to be 3D (variables, nodes, time).
"""
function (route::RapidRoute)(input::Array, params::ComponentVector; kwargs...)
    #* get the parameter types and state types
    ptyidx = get(kwargs, :ptyidx, 1:size(input, 2))
    device = get(kwargs, :device, identity)
    delta_t = get(kwargs, :delta_t, 1.0)
    #* get the interpolation type and solver type
    interp = get(kwargs, :interp, DirectInterpolation)
    solver = get(kwargs, :solver, ManualSolver(mutable=true))
    #* get the time index
    timeidx = get(kwargs, :timeidx, collect(1:size(input, 3)))
    #* get the initial states
    initstates = zeros(eltype(params), size(input, 2)) |> device

    #* var num * node num * ts len
    itpfuncs = interp(input[1, :, :], timeidx)

    #* prepare the parameters for the routing function
    expand_params = expand_component_params(params, ptyidx)[:params] |> device
    k_ps, x_ps = expand_params[:rapid_k], expand_params[:rapid_x]
    c0 = @. ((delta_t / k_ps) - (2 * x_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    c1 = @. ((delta_t / k_ps) + (2 * x_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    c2 = @. ((2 * (1 - x_ps)) - (delta_t / k_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    new_params = ComponentVector(c0=c0, c1=c1, c2=c2) |> device
    new_params_vec, new_params_axes = Vector(new_params) |> device, getaxes(new_params)
    A = (p) -> Matrix(I, size(route.adjacency)...) .- diagm(p.c0) * route.adjacency

    function du_func(u, p, t)
        q_out_t1 = u
        q_gen = itpfuncs(t)
        ps = ComponentVector(p, new_params_axes)
        #* Ax = b, x is the q_out(t+1)
        rflux_b = ps.c0 .* q_gen .+ ps.c1 .* (route.adjacency * q_out_t1 .+ q_gen) .+ ps.c2 .* q_out_t1
        #* solve the linear equation (simple solve by matrix inversion)
        A(ps) \ (rflux_b .- A(ps) * u)
    end

    #* solve the ode
    sol_arr = solver(du_func, new_params_vec, initstates, timeidx, convert_to_array=true)
    return reshape(sol_arr, 1, size(sol_arr)...)
end
