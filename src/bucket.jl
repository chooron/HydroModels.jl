"""
	HydroBucket(; funcs::Vector{<:AbstractHydroFlux}, dfuncs::Vector{<:AbstractStateFlux}=StateFlux[], name::Union{Symbol,Nothing}=nothing, sort_funcs::Bool=false)

Represents a hydrological bucket model component.

# Arguments
- `funcs::Vector{<:AbstractHydroFlux}`: A vector of flux functions that describe the hydrological processes.
- `dfuncs::Vector{<:AbstractStateFlux}`: A vector of state derivative functions (default is an empty vector of StateFlux).
- `name::Union{Symbol,Nothing}`: Optional name for the bucket. If not provided, a name will be automatically generated from state variable names.
- `sort_funcs::Bool`: Whether to sort the flux functions (default is false).

# Fields
- `funcs::Vector{<:AbstractHydroFlux}`: Vector of flux functions describing hydrological processes.
- `dfuncs::Vector{<:AbstractStateFlux}`: Vector of state derivative functions for ODE calculations.
- `flux_func::Function`: Combined function for calculating all hydrological fluxes.
- `ode_func::Union{Nothing,Function}`: Function for ordinary differential equations (ODE) calculations, or nothing if no ODE calculations are needed.
- `meta::HydroMeta`: Contains metadata about the bucket, including input, output, state, parameter, and neural network names.

# Description
HydroBucket is a structure that encapsulates the behavior of a hydrological bucket model. 
It combines multiple flux functions and state derivative functions to model water movement 
and storage within a hydrological unit.

The structure automatically extracts relevant information from the provided functions to 
populate the metadata, which includes names of:
- Inputs: Variables that drive the model
- Outputs: Variables produced by the model
- States: Internal model states that evolve over time
- Parameters: Model parameters that control behavior
- Neural Networks: Any neural network components (if applicable)

The `flux_func` and `ode_func` are constructed based on the provided `funcs` and `dfuncs`, 
enabling efficient calculation of fluxes and state changes over time.

This structure is particularly useful for building complex hydrological models by combining 
multiple HydroBucket instances to represent different components of a water system.

"""
struct HydroBucket{N, M<:ComponentVector, S} <: AbstractBucket
	"Vector of flux functions describing hydrological processes."
	fluxes::Vector{<:AbstractHydroFlux}
	"Vector of state derivative functions for ODE calculations."
	dfluxes::Vector{<:AbstractStateFlux}
	"Generated function for calculating all hydrological fluxes."
	flux_func::Function
	"Generated function for ordinary differential equations (ODE) calculations, or nothing if no ODE calculations are needed."
	ode_func::Union{Nothing, Function}
	"Metadata about the bucket, including input, output, state, parameter, and neural network names."
	meta::M

	function HydroBucket(;
		fluxes::Vector{<:AbstractHydroFlux},
		dfluxes::Vector{<:AbstractStateFlux} = StateFlux[],
		name::Union{Symbol, Nothing} = nothing,
		sort_fluxes::Bool = false,
	)
        #* sort the fluxes if needed
        fluxes = sort_fluxes ? sort_fluxes(fluxes) : fluxes
		#* Extract all variable names of fluxes and dfluxes
        bucket_inputs, bucket_outputs, bucket_states = get_all_vars(vcat(fluxes, dfluxes))
		bucket_params = reduce(union, get_param_vars.(vcat(fluxes, dfluxes)))
		bucket_nns_ntp = reduce(merge, map(flux -> NamedTuple(get_nn_vars(flux)), fluxes))
        #* Setup the meta data of the bucket
		meta = ComponentVector(inputs = bucket_inputs, outputs = bucket_outputs, states = bucket_states, params = bucket_params, nns = bucket_nns_ntp)
		#* Construct a function for ordinary differential calculation based on dfunc and funcs
		flux_func, ode_func = build_ele_func(fluxes, dfluxes, meta)
		bucket_name = isnothing(name) ? Symbol("##bucket#", hash(meta)) : name
		return new{bucket_name, typeof(meta), !isempty(bucket_states)}(
			fluxes,
			dfluxes,
			flux_func,
			ode_func,
			meta,
		)
	end
end

"""
	(ele::HydroBucket)(input::Matrix, pas::ComponentVector; config::NamedTuple=NamedTuple(), kwargs...)
	(ele::HydroBucket)(input::Array, pas::ComponentVector; config::NamedTuple=NamedTuple(), kwargs...)

Run the HydroBucket model for given input and parameters.

# Arguments
- `input`: Input data in one of these formats:
  - Matrix: dimensions are variables × time
  - Array: dimensions are variables × nodes × time 
- `pas`: ComponentVector containing model parameters and initial states
- `config`: Optional configuration with these fields:
  - `solver`: Solver to use for ODEs (default: ManualSolver{true}())
  - `interp`: Interpolation method (default: LinearInterpolation)
  - `timeidx`: Time indices (default: 1:size(input, last_dim))
  - `ptyidx`: Parameter type indices for multi-node runs
  - `styidx`: State type indices for multi-node runs

# Returns
Matrix or Array containing model outputs:
- For single node input: Matrix of size (states+outputs) × time
- For multi-node input: Array of size (states+outputs) × nodes × time

# Details
The function handles both single node and multi-node model runs:

For single node runs:
- Takes input time series for one location
- Uses provided parameters and initial states
- Solves ODEs if model has state variables
- Calculates fluxes using model's flux function
- Returns combined states and fluxes

For multi-node runs:
- Processes multiple locations simultaneously
- Can use shared or independent parameters
- Handles state propagation for each node
- Returns results for all nodes

The input dimensions must match the number of input variables defined in the model.
Required parameters and initial states must be present in the pas argument.
"""
function (ele::HydroBucket{N, M, true})(
	input::AbstractArray{T, 2},
	pas::ComponentVector;
	config::NamedTuple = NamedTuple(),
	kwargs...,
) where {N, M, T}
	#* get kwargs
	solver = get(config, :solver, ManualSolver{true}())
	interp = get(config, :interp, DataInterpolations.LinearInterpolation)
	timeidx = get(config, :timeidx, collect(1:size(input, 2)))

	initstates_vec = Vector(view(pas, :initstates))
	model_params = Vector(view(pas, :params))
	nn_params = isempty(get_nn_vars(ele)) ? Vector{eltype(pas)}[] : Vector(view(pas, :nns)) 

	#* prepare input interpolation
	itpfuncs = interp(input, timeidx)

	#* prepare parameter and nn parameter
	pas_vec = vcat(model_params, nn_params)
	params_idx_bound = 1:length(model_params)
	nn_idx_bounds = length(model_params)+1:length(pas_vec)

	#* define the ODE function
	function du_func(u, p, t)
		@views ps, nn_ps = p[params_idx_bound], p[nn_idx_bounds]
		ele.ode_func(itpfuncs(t), u, ps, nn_ps)
	end

	#* solve the problem by call the solver
	solved_states = solver(du_func, pas_vec, initstates_vec, timeidx)

	#* calculate output, slice input on time dim, then calculate each output
	tmp_flux_func(i, s) = ele.flux_func(i, s, model_params, nn_params)
	flux_output = tmp_flux_func.(eachslice(input, dims = 2), eachslice(solved_states, dims = 2))
	#* convert vector{vector} to matrix
	reduce(hcat, flux_output)
end

function (ele::HydroBucket{N, M, false})(
	input::AbstractArray{T, 2},
	pas::ComponentVector;
	kwargs...,
) where {N, M, T}
	model_params = Vector(view(pas, :params))
	nn_params = isempty(get_nn_vars(ele)) ? Vector{eltype(pas)}[] : Vector(view(pas, :nns)) 
	#* calculate output, slice input on time dim, then calculate each output
	tmp_flux_func(i) = ele.flux_func(i, nothing, model_params, nn_params)
	flux_output = tmp_flux_func.(eachslice(input, dims = 2))
	#* convert vector{vector} to matrix
	reduce(hcat, flux_output)
end

function (ele::HydroBucket{N, M, true})(
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

	#* prepare states parameters and nns
	params_len, states_len = length(get_param_names(ele)), length(get_state_names(ele))
	initstates_mat = reshape(Vector(view(pas, :initstates)), :, states_len)'
	params_mat = reshape(Vector(view(pas, :params)), :, params_len)'
	extract_initstates_mat, extract_params_mat = view(initstates_mat, :, styidx), view(params_mat, :, ptyidx)
	extract_params_vec = vec(extract_params_mat)
	nn_params = isempty(get_nn_vars(ele)) ? Vector{eltype(pas)}[] : Vector(view(pas, :nns))
	vcat_pas = vcat(extract_params_vec, nn_params)
	params_idx_bound = 1:length(extract_params_vec)
	nn_idx_bounds = length(extract_params_vec)+1:length(vcat_pas)

	#* prepare input function
	input_reshape = reshape(input, input_dims * num_nodes, time_len)
	itpfuncs  = interp(input_reshape, timeidx)

	#* define the ODE function
	function du_func(u, p, t)
		@views ps, nn_ps = reshape(view(p, params_idx_bound), params_len, num_nodes), view(p, nn_idx_bounds)
		tmp_input = reshape(itpfuncs(t), input_dims, num_nodes)
		ele.ode_func.(eachslice(tmp_input, dims = 2), eachslice(u, dims = 2), eachslice(ps, dims = 2), Ref(nn_ps))
	end

	#* Call the solve_prob method to solve the state of bucket at the specified timeidx
	solved_states = solver(du_func, vcat_pas, extract_initstates_mat, timeidx)

	#* run other functions
	tmp_flux_func(i, s, p) = ele.flux_func(i, s, p, nn_params)
	ele_output_vec = map(1:size(input, 3)) do i
		input_ = view(input, :, :, i)
		states_ = view(solved_states, :, :, i)
		reduce(hcat, ele.flux_func.(eachslice(input_, dims = 2), eachslice(states_, dims = 2), eachslice(extract_params_mat, dims = 2), Ref(nn_params)))
	end
	reduce((m1, m2) -> cat(m1, m2, dims = 3), ele_output_vec)
end

function (ele::HydroBucket{N, M, false})(
	input::AbstractArray{T, 3},
	pas::ComponentVector;
	config::NamedTuple = NamedTuple(),
	kwargs...,
) where {N, M, T}
	ptyidx = get(config, :ptyidx, 1:size(input, 2))

	#* prepare parameter and nn parameter
	params_len = length(get_param_names(ele))
	#* convert to matrix (params_len, params_types)
	params_mat = reshape(Vector(view(pas, :params)), :, params_len)'
	extract_params_mat = view(params_mat, :, ptyidx)

	nn_params_vec = isempty(get_nn_vars(ele)) ? Vector{eltype(pas)}[] : Vector(view(pas, :nns)) 
	#* run other functions
	tmp_flux_func(i, p) = ele.flux_func(i, nothing, p, nn_params_vec)
	ele_output_vec = map(1:size(input, 3)) do i
		input_ = view(input, :, :, i)
		reduce(hcat, tmp_flux_func.(eachslice(input_, dims = 2), eachslice(extract_params_mat, dims = 2)))
	end
	reduce((m1, m2) -> cat(m1, m2, dims = 3), ele_output_vec)
end
