"""
	HydroRoute(; rfunc::AbstractHydroFlux, rstate::Num, aggr_func::AbstractHydroFlux, name::Union{Symbol,Nothing}=nothing)

Represents a routing component for simulating water movement through a hydrological network.

# Arguments
- `rfunc::AbstractHydroFlux`: Flow calculation function that determines outflow from each node
- `rstate::Num`: State variable representing water storage in routing system
- `aggr_func::AbstractHydroFlux`: Flow projection function that distributes water to downstream nodes
- `name::Union{Symbol,Nothing}=nothing`: Optional identifier for the routing component

# Fields
- `rfunc::AbstractHydroFlux`: Flow calculation function
- `aggr_func::AbstractHydroFlux`: Flow projection function
- `meta::HydroMeta`: Component metadata including:
  - `inputs`: Required input variables
  - `outputs`: Generated output variables
  - `states`: State variables (from `rstate`)
  - `params`: Required parameters
  - `nns`: Neural network components (if any)

# Description
HydroRoute implements water routing through a network of connected nodes, managing both 
local flow calculations and inter-node water transfer.

## Component Structure
1. Flow Calculation (`rfunc`):
   - Computes outflow from each node
   - Uses local state and input variables
   - Can be parameter-based or neural network-based

2. Flow Projection (`aggr_func`):
   - Determines water distribution between nodes
   - Handles network connectivity
   - Maintains mass conservation

3. State Management:
   - Tracks water storage in each node
   - Updates based on inflow/outflow balance
   - Handles temporal evolution

## Implementation Types
1. Parameter-Based:
   ```julia
   # Using traditional routing equations
   route = HydroRoute(
       rfunc=HydroFlux([:storage] => [:outflow], :(k * storage)),
       rstate=:storage,
       aggr_func=HydroFlux([:outflow] => [:inflow], :outflow)
   )
   ```

2. Neural Network-Based:
   ```julia
   # Using learned flow relationships
   route = HydroRoute(
       rfunc=NeuralFlux([:storage] => [:outflow]),
       rstate=:storage,
       aggr_func=HydroFlux([:outflow] => [:inflow], :outflow)
   )
   ```

## Usage Notes
- State variables are automatically managed
- Mass conservation is enforced
- Supports both lumped and distributed modeling
- Compatible with various network topologies
- Can be combined with other hydrological components
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
        infos = (;inputs, outputs, states, params, nns)
        #* define the route name
        route_name = isnothing(name) ? Symbol("##route#", hash(infos)) : name
        #* build the route function
        multi_flux_func, multi_ode_func = build_route_func(rfluxes, dfluxes, infos)
        return new{route_name}(rfluxes, multi_flux_func, multi_ode_func, aggr_func, infos)
    end
end

"""
    @hydroroute name begin
        fluxes = begin
            ...
        end
        dfluxes = begin
            ...
        end
        aggr_func = f(x)
    end

A macro to construct a HydroRoute object.

# Arguments
- `name: The name of the route, a Symbol
- `expr`: A code block containing the following component definitions:
- `fluxes`: An array of HydroFlux or NeuralFlux objects, defining the flow calculations for the route
- `dfluxes`: An array of StateFlux objects, defining the changes in state variables
- `aggr_func`: A projection function of type Function, used to ensure state variables remain within valid ranges

# Example
```julia
route = @hydroroute :route1 begin
    fluxes = [
        @hydroflux a ~ k1 * b - k2 * c
        @hydroflux d ~ k1 * b - k2 * c
    ]
    dfluxes = [
        @stateflux c ~ b - a - d
    ]
    aggr_func = x -> max(0, x)
end
```
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
    (route::HydroRoute)(
        input::AbstractArray{T,3}, 
        params::ComponentVector;
        initstates::ComponentVector=ComponentVector(),
        config::Dict=Dict()
    ) where {T<:Number}

Execute routing simulation for a river network.

# Arguments
- `input::AbstractArray{T,3}`: Input data with dimensions (variables, nodes, timesteps)
- `params::ComponentVector`: Model parameters organized by component
- `initstates::ComponentVector=ComponentVector()`: Initial states for all nodes
- `config::Dict=Dict()`: Configuration options including:
  - `solver::AbstractHydroSolver`: ODE solver (e.g., `ManualSolver`, `Tsit5`)
  - `interp::Function`: Input interpolation method (e.g., `DirectInterpolation`)
  - `ptyidx::AbstractVector{Int}`: Parameter type indices for nodes
  - `styidx::AbstractVector{Int}`: State type indices for nodes

# Returns
- `AbstractArray{T,3}`: Output array with dimensions (variables, nodes, timesteps) containing:
  - State variables (e.g., storage)
  - Flow variables (e.g., outflow)
  - Additional derived variables

# Description
Executes a complete routing simulation by:
1. Processing inputs and parameters
2. Solving flow equations
3. Managing state evolution
4. Computing derived variables

## Simulation Steps
1. Input Processing:
   - Validate dimensions and data types
   - Apply interpolation if needed
   - Organize parameters by node

2. State Evolution:
   - Initialize state variables
   - Apply routing equations
   - Update states through time

3. Flow Computation:
   - Calculate node outflows
   - Project flows downstream
   - Maintain mass balance

## Example Usage
```julia
# Setup and run routing simulation
output = route(
    input,                          # (vars, nodes, time)
    params;                         # component parameters
    initstates=init_states,        # initial conditions
    config=Dict(
        :solver => ManualSolver{true}(),
        :interp => DirectInterpolation,
        :ptyidx => 1:n_nodes
    )
)

# Access results
storage = output[1:n_states, :, :]         # state variables
outflow = output[n_states+1:end, :, :]     # flow variables
```

# Notes
- Input dimensions must match network configuration
- Parameters must be properly structured
- States are automatically managed
- Mass conservation is enforced
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
    input_reshape = reshape(input, input_dims * num_nodes, time_len)
    itpfuncs = interp(input_reshape, timeidx)
    solved_states = solver(
        (u, p, t) -> begin
            tmp_input = reshape(itpfuncs(t), input_dims, num_nodes)
            tmp_states, tmp_outflow = route.multi_ode_func(eachslice(tmp_input, dims=1), eachslice(u, dims=1), ComponentVector(p, params_axes))
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

Implements the RAPID (Routing Application for Parallel computatIon of Discharge) routing scheme 
for large-scale river networks.

# Fields
- `adjacency::SparseMatrixCSC`: Network connectivity matrix
- `meta::HydroMeta`: Component metadata including:
  - `inputs`: Required input variables (e.g., runoff)
  - `outputs`: Generated output variables (e.g., discharge)
  - `params`: Required parameters (k, x)
  - `states`: State variables (storage)

# Constructor
    RapidRoute(;
        network::DiGraph,
        name::Union{Symbol,Nothing}=nothing
    )

# Arguments
- `network::DiGraph`: Directed graph representing river network topology
- `name::Union{Symbol,Nothing}=nothing`: Optional identifier for the routing component

# Description
RapidRoute implements the Muskingum-Cunge routing method using a state-space formulation,
designed for efficient parallel computation in large river networks.

## Mathematical Formulation
The routing scheme uses a state-space approach:
```math
Q(t+Δt) = C₁Q(t) + C₂Q(t-Δt) + C₃q(t)
```
where:
- Q: Channel discharge
- q: Lateral inflow
- C₁, C₂, C₃: Muskingum coefficients
- Δt: Time step

## Parameters
- `k`: Wave celerity parameter [T]
- `x`: Diffusion parameter [0-0.5]

## Implementation Example
```julia
# Create RAPID routing component
rapid = RapidRoute(
    network=river_network,    # DiGraph object
    name=:rapid_routing      # Optional name
)

# Run routing simulation
output = rapid(
    input,                   # Lateral inflows (vars, nodes, time)
    params;                  # Parameters (k, x for each node)
    kwargs=Dict(
        :delta_t => 3600.0,  # Time step in seconds
        :device => identity  # Computing device (CPU/GPU)
    )
)
```

# Notes
- Efficient for large-scale applications
- Parallel computation ready
- Mass conservation guaranteed
- Numerically stable formulation
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
        infos = (;inputs=inputs, outputs=outputs, states=Num[], params=[rapid_k, rapid_x])
        route_name = isnothing(name) ? Symbol("##route#", hash(infos)) : name
        #* generate adjacency matrix from network
        adjacency = adjacency_matrix(network)'
        return new{route_name}(adjacency, infos)
    end
end

"""
    (route::RapidRoute)(
        input::AbstractArray{T,3}, 
        params::ComponentVector; 
        kwargs...
    ) where {T<:Number}

Execute RAPID routing simulation for a river network.

# Arguments
- `input::AbstractArray{T,3}`: Lateral inflow data with dimensions (variables, nodes, timesteps)
- `params::ComponentVector`: Model parameters with required fields:
  - `rapid_k`: Wave celerity parameter for each node
  - `rapid_x`: Diffusion parameter for each node
- `kwargs`: Configuration options including:
  - `delta_t::Float64=1.0`: Time step in seconds
  - `device::Function=identity`: Computing device (e.g., `identity` for CPU, `cu` for GPU)
  - `interp::Function=DirectInterpolation`: Input interpolation method
  - `solver::AbstractHydroSolver=ManualSolver{true}()`: ODE solver
  - `ptyidx::AbstractVector{Int}`: Parameter type indices for nodes
  - `timeidx::AbstractVector{Int}`: Time indices for simulation

# Returns
- `AbstractArray{T,3}`: Output array with dimensions (variables, nodes, timesteps) containing:
  - Channel discharge for each node
  - Additional derived variables if specified

# Description
Executes a RAPID routing simulation using the Muskingum-Cunge method:

1. Parameter Processing:
   - Computes Muskingum coefficients (C₀, C₁, C₂)
   - Handles parameter distribution across nodes
   - Prepares sparse matrix operations

2. State Evolution:
   - Updates discharge states using state-space formulation
   - Maintains numerical stability
   - Ensures mass conservation

3. Device Management:
   - Supports CPU and GPU computation
   - Efficient memory handling
   - Automatic device placement

## Example Usage
```julia
# Basic simulation
output = route(
    input,                # Lateral inflows
    params;              # k and x parameters
    delta_t=3600.0,     # 1-hour timestep
    device=identity     # Run on CPU
)

# Advanced configuration
output = route(
    input,
    params;
    delta_t=1800.0,                    # 30-min timestep
    device=gpu_device(),                         # Run on GPU
    interp=DirectInterpolation,        # Linear interpolation
    solver=ManualSolver{true}(),       # Manual timestepping
    ptyidx=1:n_nodes,                 # All nodes
    timeidx=1:n_timesteps             # Full simulation
)
```

# Notes
- Input dimensions must match network size
- Parameters must be properly structured
- Mass conservation is guaranteed
- GPU acceleration supported
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
