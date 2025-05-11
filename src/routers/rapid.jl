module RapidRoute

using ..HydroModels
using ..HydroModels: AbstractRoute
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
    expand_params = expand_component_params(params, get_param_names(route), ptyidx)[:params] |> device
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
    sol_arr = solver(du_func, new_params_vec, initstates, timeidx)
    return reshape(sol_arr, 1, size(sol_arr)...)
end

end