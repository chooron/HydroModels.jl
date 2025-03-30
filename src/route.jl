"""
	HydroRoute(; rfunc::AbstractHydroFlux, rstate::Num, proj_func::AbstractHydroFlux, name::Union{Symbol,Nothing}=nothing)

Represents a routing structure for hydrological modeling.

# Arguments
- `rfunc::AbstractHydroFlux`: The routing function used for flow calculations.
- `rstate::Num`: The state variable for routing.
- `proj_func::AbstractHydroFlux`: Function for projecting outflow to downstream nodes.
- `name::Union{Symbol,Nothing}=nothing`: Optional name for the routing instance. If not provided, will be automatically generated.

# Fields
- `rfunc::AbstractHydroFlux`: The routing function used for flow calculations.
- `proj_func::Function`: Function for projecting outflow to downstream nodes.
- `meta::HydroMeta`: Contains metadata about the routing instance, including input, output, state, parameter and neural network names.

# Description
HydroRoute is a structure that represents a routing system in a hydrological model.
It uses a specified routing function (`rfunc`) to calculate flow between nodes and a 
projection function (`proj_func`) to determine how water moves between connected nodes.

The routing process involves:
1. Calculating outflow from each node using the routing function
2. Projecting outflows to downstream nodes using the projection function
3. Updating node states based on inflow, outflow and locally generated runoff

The structure supports both traditional parameter-based routing functions and neural network
based routing functions through the AbstractHydroFlux interface.

The metadata (`meta`) is automatically constructed from the provided functions and contains:
- Input names (excluding the routing state variable)
- Output names
- Parameter names
- State name (from `rstate`)
- Neural network names (if any)

This structure serves as the base for more specific routing implementations like GridRoute
and VectorRoute.
"""
struct HydroRoute <: AbstractHydroRoute
    "Name of the routing instance"
    name::Symbol
    "Routing function"
    rfluxes::Vector{<:AbstractFlux}
    "Generated function for calculating all hydrological fluxes."
    multi_flux_func::Function
    "Generated function for ordinary differential equations (ODE) calculations, or nothing if no ODE calculations are needed."
    multi_ode_func::Function
    "Outflow projection function"
    proj_func::Function
    "Metadata: contains keys for input, output, param, state, and nn"
    infos::NamedTuple

    function HydroRoute(;
        rfluxes::Vector{<:AbstractFlux},
        dfluxes::Vector{<:AbstractStateFlux},
        proj_func::Function,
        name::Union{Symbol,Nothing}=nothing,
    )
        #* Extract all variable names of funcs and dfuncs
        input_names, output_names, state_names = get_var_names(rfluxes, dfluxes)
        param_names = reduce(union, get_param_names.(rfluxes))
        nn_names = reduce(union, get_nn_names.(rfluxes))
        #* Setup the name information of the hydrobucket
        infos = (;inputs=input_names, outputs=output_names, states=state_names, params=param_names, nns=nn_names)
        #* define the route name
        route_name = isnothing(name) ? Symbol("##route#", hash(infos)) : name
        #* build the route function
        multi_flux_func, multi_ode_func = build_route_func(rfluxes, dfluxes, infos)
        return new(route_name, rfluxes, multi_flux_func, multi_ode_func, proj_func, infos)
    end
end

"""
	GridRoute(;
		rfluxes::Vector{<:AbstractHydroFlux},
        dfluxes::Vector{<:AbstractStateFlux},
        proj_func::Function,
		name::Union{Symbol,Nothing}=nothing
	)

Create a HydroRoute instance for grid-based river routing.

# Arguments
- `rfunc::AbstractHydroFlux`: Routing function that calculates outflow from each node
- `rstate::Num`: Symbolic variable representing the routing state
- `flwdir::AbstractMatrix`: Flow direction matrix using D8 encoding (1-128)
- `positions::AbstractVector`: Vector of (row,col) positions for each node in the grid
- `aggtype::Symbol=:matrix`: Aggregation type for flow routing:
	- `:matrix`: Uses graph-based adjacency matrix
	- `:network`: Uses network-based adjacency matrix
- `name::Union{Symbol,Nothing}=nothing`: Optional name for the routing instance

# Returns
`HydroRoute`: A configured routing instance for grid-based networks

# Description
Creates a routing structure for grid-based river networks using either a matrix-based
approach (matrix) or network-based approach (network) for flow accumulation. The function
verifies that node positions match node IDs and constructs appropriate projection
functions based on the chosen aggregation type.
"""
function GridRoute(;
    rfluxes::Vector{<:AbstractFlux},
    dfluxes::Vector{<:AbstractStateFlux},
    flwdir::AbstractMatrix,
    positions::AbstractVector,
    aggtype::Symbol=:matrix,
    name::Union{Symbol,Nothing}=nothing,
)
    if aggtype == :matrix
        d8_codes = [1, 2, 4, 8, 16, 32, 64, 128]
        d8_nn_pads = [(1, 1, 2, 0), (2, 0, 2, 0), (2, 0, 1, 1), (2, 0, 0, 2), (1, 1, 0, 2), (0, 2, 0, 2), (0, 2, 1, 1), (0, 2, 2, 0)]

        #* input dims: node_num * ts_len
        function grid_routing(input::AbstractVector, positions::AbstractVector, flwdir::AbstractMatrix)
            #* Convert input to sparse matrix
            input_arr = Array(sparse([pos[1] for pos in positions], [pos[2] for pos in positions], input, size(flwdir)[1], size(flwdir)[2]))
            #* Calculate weighted summation
            input_routed = sum(collect([pad_zeros(input_arr .* (flwdir .== code), arg) for (code, arg) in zip(d8_codes, d8_nn_pads)]))
            #* Clip input matrix border
            clip_arr = input_routed[2:size(input_arr)[1]+1, 2:size(input_arr)[2]+1]
            #* Convert input matrix to vector
            collect([clip_arr[pos[1], pos[2]] for pos in positions])
        end
        #* build the outflow projection function
        proj_func = (outflow) -> grid_routing(outflow, positions, flwdir)

    elseif aggtype == :network
        network = build_grid_digraph(flwdir, positions)
        #* build the outflow projection function
        adjacency = adjacency_matrix(network)'
        proj_func = (outflow) -> adjacency * outflow
    else
        @error "the $aggtype is not support"
    end

    return HydroRoute(; rfluxes, dfluxes, proj_func, name)
end

"""
	VectorRoute(;
		rfluxes::Vector{<:AbstractHydroFlux},
		rvars::Vector{Num},
        rstates::Vector{Num},
        proj_func::Function,
		network::DiGraph,
		name::Union{Symbol,Nothing}=nothing
	)

Create a HydroRoute instance for vector-based river routing.

# Arguments
- `rfunc::AbstractHydroFlux`: Routing function that calculates outflow from each node
- `rstate::Num`: Symbolic variable representing the routing state
- `network::DiGraph`: Directed graph representing the river network connectivity
- `name::Union{Symbol,Nothing}=nothing`: Optional name for the routing instance

# Returns
`HydroRoute`: A configured routing instance for vector-based networks

# Description
Creates a routing structure for vector-based river networks using a graph-based approach
for flow accumulation. The function verifies that the number of nodes matches the number
of node IDs and constructs a projection function based on the network's adjacency matrix.
The adjacency matrix is used to route flow between connected nodes in the network.
"""
function VectorRoute(;
    rfluxes::Vector{<:AbstractFlux},
    dfluxes::Vector{<:AbstractStateFlux},
    network::DiGraph,
    name::Union{Symbol,Nothing}=nothing,
)
    adjacency = adjacency_matrix(network)'
    proj_func = (outflow) -> adjacency * outflow
    return HydroRoute(; rfluxes, dfluxes, proj_func, name)
end

"""
	(route::HydroRoute)(input::Array, pas::ComponentVector; config::NamedTuple=NamedTuple(), kwargs...)

Run the routing model for given input and parameters.

# Arguments
- `input`: Input data array with dimensions (variables × nodes × time)
- `pas`: ComponentVector containing parameters, initial states and neural network parameters (if applicable)
- `config`: Configuration options including:
  - `ptypes`: Parameter types to use for each node (default: all parameter types)
  - `interp`: Interpolation method for input data (default: LinearInterpolation)
  - `solver`: AbstractHydroSolver to use for ODE solving (default: ODESolver())
  - `timeidx`: Vector of time indices (default: 1:size(input,3))

# Returns
Array with dimensions (states+outputs × nodes × time) containing:
- Routing states for each node over time
- Routed outflow for each node over time

# Description
This function executes the routing model by:
1. Validating input dimensions and parameter/state configurations
2. Setting up interpolation and solver configurations
3. Building parameter functions for either neural network or regular routing
4. Solving the routing equations using the specified solver
5. Computing outflows based on states and parameters

The routing can use either neural network based routing functions (AbstractNeuralFlux) or
regular routing functions, with parameters extracted accordingly.
"""
function (route::HydroRoute)(input::AbstractArray{T,3}, params::ComponentVector; kwargs...,) where {T}
    input_dims, num_nodes, time_len = size(input)

    #* get kwargs
    ptyidx = get(kwargs, :ptyidx, 1:num_nodes)
    styidx = get(kwargs, :styidx, 1:num_nodes)
    interp = get(kwargs, :interp, LinearInterpolation)
    solver = get(kwargs, :solver, ManualSolver{true}())
    timeidx = get(kwargs, :timeidx, collect(1:time_len))
    initstates = get(kwargs, :initstates, zeros(eltype(params), length(get_state_names(route)), num_nodes))

    #* prepare states parameters and nns
    new_pas = expand_component_params(params, ptyidx)
    initstates_mat = expand_component_initstates(initstates, styidx)

    #* prepare input function
    input_reshape = reshape(input, input_dims * num_nodes, time_len)
    itpfuncs = interp(input_reshape, timeidx)
    solved_states = solver(
        (u, p, t) -> begin
            tmp_input = reshape(itpfuncs(t), input_dims, num_nodes)
            tmp_states, tmp_outflow = route.multi_ode_func(eachslice(tmp_input, dims=1), eachslice(u, dims=1), p)
            tmp_states_arr = reduce(hcat, tmp_states)
            tmp_inflow_arr = reduce(hcat, route.proj_func.(tmp_outflow))
            # todo 这里元编程表达一直存在问题
            tmp_states_arr .+ tmp_inflow_arr |> permutedims
        end,
        new_pas, initstates_mat, timeidx
    )
    #* run other functions
    output = route.multi_flux_func(eachslice(input, dims=1), eachslice(solved_states, dims=1), new_pas)
    output_arr = if length(output) > 1
        permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), output), (3, 1, 2))
    else
        reshape(output[1], 1, num_nodes, time_len)
    end
    cat(solved_states, output_arr, dims=1)
end

"""
	RapidRoute<: AbstractRoute

A structure representing a vector-based routing scheme for hydrological modeling.

# Fields
- `rfunc::AbstractVector{<:AbstractRouteFlux}`: Vector of routing flux functions for each node.
- `network::DiGraph`: A directed graph representing the routing network topology.
- `infos::NamedTuple`: Contains information about the VectorRoute instance, including input, output, state, and parameter names.

# Constructor
	RapidRoute(
		name::Symbol;
		rfunc::AbstractRouteFlux,
		network::DiGraph
	)

Constructs a `RapidRoute` object with the given name, routing flux function, and network structure.

# Arguments
- `name::Symbol`: A symbol representing the name of the routing scheme.
- `rfunc::AbstractRouteFlux`: The routing flux function to be applied at each node.
- `network::DiGraph`: A directed graph representing the routing network topology.

The constructor extracts variable names, parameter names, and neural network names from the provided
routing flux function to set up the internal information structure of the `RapidRoute` object.

Note: from Rapid
"""
struct RapidRoute <: AbstractRoute
    "Name of the routing instance"
    name::Symbol
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
        input_names, output_names = tosymbol.(inputs), tosymbol.(outputs)
        @assert length(inputs) == length(outputs) == 1 "The length of inputs and outputs must be the 1, but got inputs: $(length(inputs)) and outputs: $(length(outputs))"
        #* Setup the name information of the hydrobucket
        infos = (;inputs=input_names, outputs=output_names, states=Symbol[], params=[:rapid_k, :rapid_x])
        route_name = isnothing(name) ? Symbol("##route#", hash(infos)) : name
        #* generate adjacency matrix from network
        adjacency = adjacency_matrix(network)'
        return new(route_name, adjacency, infos)
    end
end

function (route::RapidRoute)(input::Array, params::ComponentVector; kwargs...)
    #* get the parameter types and state types
    ptyidx = get(kwargs, :ptyidx, 1:size(input, 2))
    delta_t = get(kwargs, :delta_t, 1.0)
    #* get the interpolation type and solver type
    interp = get(kwargs, :interp, LinearInterpolation)
    solver = get(kwargs, :solver, ManualSolver{true}())
    #* get the time index
    timeidx = get(kwargs, :timeidx, collect(1:size(input, 3)))
    #* get the initial states
    initstates = get(kwargs, :initstates, zeros(eltype(params), size(input, 2)))

    #* var num * node num * ts len
    itpfuncs = interp(input[1, :, :], timeidx)

    #* prepare the parameters for the routing function
    expand_params = expand_component_params(params, ptyidx)
    k_ps = view(expand_params, :rapid_k)
    x_ps = view(expand_params, :rapid_x)
    c0 = @. ((delta_t / k_ps) - (2 * x_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    c1 = @. ((delta_t / k_ps) + (2 * x_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    c2 = @. ((2 * (1 - x_ps)) - (delta_t / k_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    new_params = ComponentVector(c0=c0, c1=c1, c2=c2)
    new_params_vec, new_params_axes = Vector(new_params), getaxes(new_params)
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
