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
struct HydroRoute{N, M <: ComponentVector} <: AbstractHydroRoute
	"Routing function"
	rfunc::AbstractFlux
	"Outflow projection function"
	route_func::Function
	proj_func::Function
	"Metadata: contains keys for input, output, param, state, and nn"
	meta::M

	function HydroRoute(;
		rfunc::AbstractFlux,
		rstates::Vector{Num},
		proj_func::Function,
		name::Union{Symbol, Nothing} = nothing,
	)
		#* Extract all variable names of funcs and dfuncs
		inputs, outputs = get_input_vars(rfunc), get_output_vars(rfunc)
		params, nns = get_param_vars(rfunc), get_nn_vars(rfunc)
		#* Setup the name information of the hydrobucket
		meta = ComponentVector(inputs = setdiff(inputs, rstates), outputs = outputs, states = rstates, params = params, nns = nns)
		#* define the route name
		route_name = isnothing(name) ? Symbol("##route#", hash(meta)) : name
		#* build the route function
		route_func = build_route_func(rfunc, rstates)

		return new{route_name, typeof(meta)}(
			rfunc,
			route_func,
			proj_func,
			meta,
		)
	end
end

"""
	GridRoute(;
		rfunc::AbstractHydroFlux,
		rstate::Num,
		flwdir::AbstractMatrix,
		positions::AbstractVector,
		aggtype::Symbol=:type1,
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
	rfunc::AbstractHydroFlux,
	rstates::Vector{Num},
	flwdir::AbstractMatrix,
	positions::AbstractVector,
	aggtype::Symbol = :matrix,
	name::Union{Symbol, Nothing} = nothing,
)
	if aggtype == :matrix
		d8_codes = [1, 2, 4, 8, 16, 32, 64, 128]
		d8_nn_pads = [(1, 1, 2, 0), (2, 0, 2, 0), (2, 0, 1, 1), (2, 0, 0, 2), (1, 1, 0, 2), (0, 2, 0, 2), (0, 2, 1, 1), (0, 2, 2, 0)]

		#* input dims: node_num * ts_len
		function grid_routing(input::AbstractVector, positions::AbstractVector, flwdir::AbstractMatrix)
			#* 转换为input的稀疏矩阵
			input_arr = Array(sparse([pos[1] for pos in positions], [pos[2] for pos in positions], input, size(flwdir)[1], size(flwdir)[2]))
			#* 计算权重求和结果
			input_routed = sum(collect([pad_zeros(input_arr .* (flwdir .== code), arg) for (code, arg) in zip(d8_codes, d8_nn_pads)]))
			#* 裁剪输入矩阵边框
			clip_arr = input_routed[2:size(input_arr)[1]+1, 2:size(input_arr)[2]+1]
			#* 将输入矩阵转换为向量
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

	return HydroRoute(; rfunc, rstates, proj_func, name)
end

"""
	VectorRoute(;
		rfunc::AbstractHydroFlux,
		rstate::Num,
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
	rfunc::AbstractHydroFlux,
	rstates::Vector{Num},
	network::DiGraph,
	name::Union{Symbol, Nothing} = nothing,
)
	#* generate adjacency matrix from network
	adjacency = adjacency_matrix(network)'
	#* build the outflow projection function
	proj_func = (outflow) -> adjacency * outflow
	return HydroRoute(; rfunc, rstates, proj_func, name)
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
function (route::HydroRoute{N, M})(
	input::AbstractArray{T, 3},
	pas::ComponentVector;
	config::NamedTuple = NamedTuple(),
	kwargs...,
) where {N, M, T}
	input_dims, num_nodes, time_len = size(input)
	#* get the parameter types and state types
	ptyidx = get(config, :ptyidx, 1:size(input, 2))
	styidx = get(config, :styidx, 1:size(input, 2))
	#* get the interpolation type and solver type
	interp = get(config, :interp, LinearInterpolation)
	solver = get(config, :solver, ManualSolver{true}())
	#* get the time index
	timeidx = get(config, :timeidx, collect(1:size(input, 3)))

	#* prepare parameter and nn parameter
	params_len = length(get_param_names(route))
	states_len = length(get_state_names(route))
	#* convert to matrix (params_len/states_len, params_types/state_types)
	initstates_mat = reshape(Vector(view(pas, :initstates)), :, states_len)'
	params_mat = reshape(Vector(view(pas, :params)), :, params_len)'
	extract_initstates_mat = view(initstates_mat, :, styidx)
	extract_params_mat = view(params_mat, :, ptyidx)
	nn_params_vec = isempty(get_nn_vars(route)) ? Vector{eltype(pas)}[] : Vector(view(pas, :nns))
	vcat_pas = vcat(nn_params_vec, vec(extract_params_mat))
	nn_idx_bounds = 1:length(nn_params_vec)
	params_idx_bound = length(nn_params_vec)+1:length(vcat_pas)

	#* prepare input function
	input_reshape = reshape(input, input_dims * num_nodes, time_len)
	itpfunc_list = interp.(eachslice(input_reshape, dims = 1), Ref(timeidx))
	ode_input_func(t) = map(itpfunc -> itpfunc(t), itpfunc_list)

	#* define the ODE function
	function du_func(u, p, t)
		@views ps, nn_ps = reshape(view(p, params_idx_bound), params_len, num_nodes), view(p, nn_idx_bounds)
		route_output = route.route_func.(ode_input_func(t), eachslice(u, dims = 2), eachslice(ps, dims = 2), Ref(nn_ps))
		route_output_mat = reduce(hcat, route_output)
		q_gen, q_out = view(route_output_mat, 1, :), view(route_output_mat, 2, :)
		q_in = route.proj_func(q_out)
		q_in .+ q_gen .- q_out
	end

	#* Call the solve_prob method to solve the state of bucket at the specified timeidx
	solved_states = solver(du_func, vcat_pas, extract_initstates_mat, timeidx)
	#* run other functions
	route_output_vec = map(1:size(input, 3)) do i
		input_ = @view input[:, :, i]
		states_ = @view solved_states[:, :, i]
		reduce(hcat, route.route_func.(eachslice(input_, dims = 2), eachslice(states_, dims = 2), eachslice(extract_params_mat, dims = 2), Ref(nn_params_vec)))
	end
	route_output_arr = reduce((m1, m2) -> cat(m1, m2, dims = 3), route_output_vec)
	cat(solved_states, route_output_arr, dims = 1)
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
struct RapidRoute{N, M<:ComponentVector} <: AbstractRoute
	"Routing adjacency matrix"
	adjacency::AbstractMatrix
	"Metadata: contains keys for input, output, param, state, and nn"
	meta::M

	function RapidRoute(
		fluxes::Pair{Vector{Num}, Vector{Num}};
		network::DiGraph,
		name::Union{Symbol, Nothing} = nothing,
	)
		#* Extract all variable names of funcs and dfuncs
		inputs, outputs = fluxes[1], fluxes[2]
		@assert length(inputs) == length(outputs) == 1 "The length of inputs and outputs must be the 1, but got inputs: $(length(inputs)) and outputs: $(length(outputs))"
		#* Extract all parameters names of funcs and dfuncs
		@parameters rapid_k rapid_x
		#* Setup the name information of the hydrobucket
		meta = ComponentVector(inputs = inputs, outputs = outputs, states = Num[], params = [rapid_k, rapid_x])
		route_name = isnothing(name) ? Symbol("##route#", hash(meta)) : name
		#* generate adjacency matrix from network
		adjacency = adjacency_matrix(network)'
		return new{route_name, typeof(meta)}(
			adjacency,
			meta,
		)
	end
end

function (route::RapidRoute)(
	input::Array,
	pas::ComponentVector;
	config::NamedTuple = NamedTuple(),
	kwargs...,
)
	input_dims, num_nodes, time_len = size(input)
	#* get the parameter types and state types
	ptyidx = get(config, :ptyidx, 1:size(input, 2))
	delta_t = get(config, :delta_t, 1.0)
	#* get the interpolation type and solver type
	interp = get(config, :interp, LinearInterpolation)
	solver = get(config, :solver, ManualSolver{true}())
	#* get the time index
	timeidx = get(config, :timeidx, collect(1:size(input, 3)))

	#* var num * node num * ts len
	itp_funcs = interp.(eachslice(input[1, :, :], dims = 1), Ref(timeidx))

	#* prepare the parameters for the routing function
	params = view(pas, :params)
	k_ps = view(view(params, :rapid_k), ptyidx)
	x_ps = view(view(params, :rapid_x), ptyidx)
	c0 = @. ((delta_t / k_ps) - (2 * x_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
	c1 = @. ((delta_t / k_ps) + (2 * x_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
	c2 = @. ((2 * (1 - x_ps)) - (delta_t / k_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
	A = (p) -> Matrix(I, size(route.adjacency)...) .- diagm(p.c0) * route.adjacency

	function du_func(u, p, t)
		q_out_t1 = u
		q_gen = [itp_func(t) for itp_func in itp_funcs]
		#* Ax = b, x is the q_out(t+1)
		rflux_b = p.c0 .* q_gen .+ p.c1 .* (route.adjacency * q_out_t1 .+ q_gen) .+ p.c2 .* q_out_t1
		#* solve the linear equation (simple solve by matrix inversion)
		A(p) \ (rflux_b .- A(p) * u)
	end

	#* solve the ode
	sol_arr = solver(du_func, ComponentVector(c0 = c0, c1 = c1, c2 = c2), zeros(size(input)[2]), timeidx, convert_to_array = true)
	return reshape(sol_arr, 1, size(sol_arr)...)
end
