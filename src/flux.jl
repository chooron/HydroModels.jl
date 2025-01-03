"""
	HydroFlux

Represents a simple flux component in a hydrological model.

# Arguments
- `inputs::Vector{Num}`: A vector of input variables.
- `outputs::Vector{Num}`: A vector of output variables. 
- `params::Vector{Num}`: A vector of parameter variables.
- `exprs::Vector{Num}`: A vector of expressions describing the formulas for output variables.
- `name::Union{Symbol,Nothing}=nothing`: Optional name for the flux component.

# Fields
- `name::Symbol`: Name of the flux component
- `exprs::Vector{Num}`: Vector of expressions describing the formulas for output variables
- `func::Function`: Compiled function that calculates the flux
- `meta::ComponentVector`: Metadata containing input, output, and parameter variables

# Constructors
	HydroFlux(inputs::Vector{Num}, outputs::Vector{Num}, params::Vector{Num}; 
			  exprs::Vector{Num}, name::Union{Symbol,Nothing}=nothing)
	HydroFlux(fluxes::Pair{Vector{Num},Vector{Num}}, params::Vector{Num}=Num[]; 
			  exprs::Vector{Num}, name::Union{Symbol,Nothing}=nothing)

# Description
`HydroFlux` is a structure that encapsulates a simple flux calculation in a hydrological model.
It can be constructed either by providing explicit inputs, outputs, parameters, and expressions,
or by specifying input/output fluxes as a pair along with parameters and expressions.

The structure automatically compiles the provided expressions into an efficient calculation 
function, which can be used to compute flux values given input and parameter values.

If no name is provided, a unique name is generated based on the hash of the expressions.
The number of expressions must match the number of output variables.

This structure is particularly useful for representing straightforward hydrological processes
where the relationship between inputs and outputs can be expressed as simple mathematical formulas.
"""
struct HydroFlux{N, M <: ComponentArray} <: AbstractHydroFlux
	"Vector of expressions describing the formulas for output variables"
	exprs::Vector{Num}
	"Compiled function that calculates the flux"
	func::Function
	"Metadata about the flux, including input, output, and parameter names"
	meta::M

	function HydroFlux(
		inputs::Vector{T},
		outputs::Vector{T},
		params::Vector{T};
		exprs::Vector{T},
		name::Union{Symbol, Nothing} = nothing,
	) where {T <: Num}
		#* construct meta
		meta = ComponentVector(inputs = inputs, outputs = outputs, params = params)
		@assert length(exprs) == length(outputs) "The number of expressions and outputs must match, but got expressions: $(length(exprs)) and outputs: $(length(outputs))"

		#* build flux function
		flux_func = build_flux_func(inputs, outputs, params, exprs)

        #* use hash of exprs to name the flux
		flux_name = isnothing(name) ? Symbol("##hydro_flux#", hash(meta)) : name
		return new{flux_name, typeof(meta)}(exprs, flux_func, meta)
	end

	#* construct hydro flux with input fluxes and output fluxes
	function HydroFlux(
		fluxes::Pair{Vector{Num}, Vector{Num}},
		params::Vector{Num} = Num[];
		exprs::Vector{Num},
		name::Union{Symbol, Nothing} = nothing,
	)
		return HydroFlux(fluxes[1], fluxes[2], params, exprs = exprs, name = name)
	end
end

"""
	(flux::AbstractHydroFlux)(input::Union{Vector,Matrix,Array}, pas::ComponentVector; config::NamedTuple=NamedTuple(), kwargs...)

Apply the simple flux model to input data of various dimensions.

# Arguments
- `input`: Input data, which can be:
  - `Vector`: A vector of input values for a single time step.
  - `Matrix`: A matrix of input values, where each column represents a different time step.
  - `Array`: A 3D array of input values, with dimensions (var_names, node_names, ts_len).
- `pas::ComponentVector`: A component vector containing parameter values.
- `config::NamedTuple`: Configuration options for the flux calculation:
  - `ptyidx`: Indices specifying which parameter types to use (only for 3D array input)
- `kwargs...`: Additional keyword arguments (unused in this function), provided for compatibility with the component callable function API

# Returns
- For vector input: The result of applying the flux function to the input and parameters.
- For matrix input: A matrix where each column is the result of applying the flux function to the corresponding input column.
- For 3D array input: A 3D array of flux outputs, with dimensions (output_var_names, node_names, ts_len).
"""
function (flux::AbstractHydroFlux)(input::AbstractVector, pas::ComponentVector; kwargs...)
	params_vec = Vector(view(pas, :params))
	flux.func(input, params_vec)
end

function (flux::AbstractHydroFlux)(input::AbstractArray{T, 2}, pas::ComponentVector; kwargs...) where {T}
	params_vec = Vector(view(pas, :params))
	output_arr = reduce(hcat, flux.func.(eachslice(input, dims = 2), Ref(params_vec)))
	output_arr
end

function (flux::AbstractHydroFlux)(input::AbstractArray{T, 3}, pas::ComponentVector; config::NamedTuple = NamedTuple(), kwargs...) where {T}
	ptyidx = get(config, :ptyidx, 1:size(input, 2))
	#* prepare parameter and nn parameter
	params_len = length(get_param_names(flux))

	#* convert to matrix (params_len, params_types)
	params_mat = reshape(Vector(view(pas, :params)), :, params_len)'
	extract_params_mat = view(params_mat, :, ptyidx)

	#* array dims: (var_names * node_names * ts_len)
	output_vec = map(1:size(input, 3)) do i
		input_ = @view input[:, :, i]
		reduce(hcat, flux.func.(eachslice(input_, dims = 2), eachslice(extract_params_mat, dims = 2)))
	end
	#* if there is only one time step, return the output as a matrix
	if length(output_vec) == 1
		tmp_output_arr = output_vec[1]
		return reshape(tmp_output_arr, size(tmp_output_arr)..., 1)
	else
		return reduce((m1, m2) -> cat(m1, m2, dims = 3), output_vec)
	end
end


"""
	NeuralFlux <: AbstractNeuralFlux

Represents a neural network-based flux component in a hydrological model.

# Fields
- `expr::Symbolics.Arr{Num,1}`: Expressions describing the formulas for output variables.
- `func::Function`: A compiled function that calculates the flux using the neural network.
- `meta::ComponentVector`: Contains metadata about the flux, including input, output, and neural network parameter names.
- `nninfos::NamedTuple`: Contains information about the neural network's input and output structure.

# Constructors
	NeuralFlux(inputs::Vector{Num}, outputs::Vector{Num}, chain::LuxCore.AbstractLuxLayer)

# Description
`NeuralFlux` is a structure that encapsulates a neural network-based flux calculation in a hydrological model. 
It combines symbolic mathematics with neural networks to create a flexible and powerful representation of complex hydrological processes.

The structure automatically handles the integration of the neural network into the symbolic framework,
allowing for seamless incorporation of machine learning models into traditional hydrological equations.

The constructor requires:
- A vector of input variables
- A vector of output variables  
- A Lux neural network chain with a specified name (must be a Symbol)

The structure initializes the neural network parameters and creates the necessary symbolic variables and expressions
to integrate the neural network into the hydrological model framework.

This structure is particularly useful for representing complex, non-linear relationships in hydrological systems
where traditional equations may be insufficient or unknown.
"""
struct NeuralFlux{N, M <: ComponentVector} <: AbstractNeuralFlux
	"Array of expressions describing the formulas for output variables"
	exprs::Vector
	"Compiled function that calculates the flux using the neural network"
	func::Function
	"Metadata about the flux, including input, output, and neural network parameter names"
	meta::M
	"Information about the neural network's input and output structure"
	nninfos::NamedTuple

	function NeuralFlux(
		inputs::Vector{T},
		outputs::Vector{T},
		chain::LuxCore.AbstractLuxLayer;
		name::Union{Symbol, Nothing} = nothing,
	) where {T <: Num}
		#* Check chain name
		@assert chain.name isa Symbol "Neural network chain should have a name with Symbol type"
		chain_name = chain.name

		#* Initialize parameters
		init_params = ComponentVector(Lux.initialparameters(StableRNG(42), chain))
		params_axes = getaxes(init_params)

		#* Define parameter variables 
		chain_params = first(@parameters $chain_name[1:length(init_params)] [guess = Vector(init_params), description = "$chain_name Neural network parameters"])
		lazy_params = Symbolics.array_term((x, axes) -> ComponentVector(x, axes), chain_params, params_axes, size = size(chain_params))

		#* Define neural network input/output variables
		nn_input_name = Symbol(chain_name, :_input)
		nn_output_name = Symbol(chain_name, :_output)
		nn_input = first(@variables $(nn_input_name)[1:length(inputs)] [description = "$chain_name Neural network input"])
		nn_output = first(@variables $(nn_output_name)[1:length(outputs)] [description = "$chain_name Neural network output"])

		#* Build expressions and function
		flux_expr = LuxCore.stateless_apply(chain, nn_input, lazy_params)
		nn_func = (x, p) -> LuxCore.stateless_apply(chain, x, ComponentVector(p, params_axes))

		#* Create metadata
		meta = ComponentVector(inputs = inputs, outputs = outputs, nns = NamedTuple{Tuple([chain_name])}([chain_params]))
		nninfos = (inputs = nn_input, outputs = nn_output, nns = chain_params)

		flux_name = isnothing(name) ? Symbol("##neural_flux#", meta) : name

		new{flux_name, typeof(meta)}(
			[flux_expr], nn_func,
			meta, nninfos,
		)
	end

	#* construct neural flux with input fluxes and output fluxes
	function NeuralFlux(fluxes::Pair{Vector{Num}, Vector{Num}}, chain, name::Union{Symbol, Nothing} = nothing)
		return NeuralFlux(fluxes[1], fluxes[2], chain, name = name)
	end
end

"""
	(flux::AbstractFlux)(input::Union{Vector,Matrix,Array}, pas::ComponentVector; ptypes::AbstractVector{Symbol}=Symbol[], kwargs...)

Apply the flux model (simple or neural) to input data of various dimensions.

# Arguments
- `input`: Input data, which can be:
  - `Vector`: A vector of input values for a single time step.
  - `Matrix`: A matrix of input values, where each column represents a different time step.
  - `Array`: A 3D array of input values, with dimensions (var_names, node_names, ts_len).
- `pas::ComponentVector`: A component vector containing parameter values.
- `ptypes::AbstractVector{Symbol}`: A vector of symbols representing parameter categories (only used for `Array` input).
- `kwargs...`: Additional keyword arguments (unused in this function), provided for compatibility with the component callable function API

# Returns
- For vector input: The result of applying the flux function to the input and parameters.
- For matrix input: A matrix where each column is the result of applying the flux function to the corresponding input column.
- For 3D array input: A 3D array of flux outputs, with dimensions (output_var_names, node_names, ts_len).

# Note
For neural flux models, the parameters are accessed from `pas[:nn]` instead of `pas[:params]`.
"""
function (flux::AbstractNeuralFlux)(input::AbstractVector, pas::ComponentVector; kwargs...)
	nn_params_vec = pas[:nns][get_nn_names(flux)[1]]
	flux.func(input, nn_params_vec)
end

function (flux::AbstractNeuralFlux)(input::AbstractArray{T, 2}, pas::ComponentVector; kwargs...) where {T}
	nn_params_vec = pas[:nns][get_nn_names(flux)[1]]
	output_arr = flux.func(input, nn_params_vec)
	output_arr
end

function (flux::AbstractNeuralFlux)(input::AbstractArray{T, 3}, pas::ComponentVector; kwargs...) where {T}
	nn_params = pas[:nns][get_nn_names(flux)[1]]
	#* array dims: (ts_len * node_names * var_names)
	flux_output_vec = [flux.func(input[:, i, :], nn_params) for i in 1:size(input)[2]]
	flux_output_arr = reduce((m1, m2) -> cat(m1, m2, dims = 3), flux_output_vec)
	permutedims(flux_output_arr, (1, 3, 2))
end

"""
    StateFlux <: AbstractStateFlux

Represents a state flux component in a hydrological model.

# Fields
- `name::Symbol`: Name of the state flux component
- `expr::Num`: Expression describing the formula for the state variable
- `meta::ComponentVector`: Metadata containing input, state, and parameter variables

# Constructors
    StateFlux(inputs::Vector{Num}, state::Num, params::Vector{Num}=Num[]; 
              expr::Num, name::Union{Symbol,Nothing}=nothing)
    StateFlux(fluxes::Pair{Vector{Num},Vector{Num}}, state::Num)
    StateFlux(states::Pair{Num,Num})

# Description
StateFlux is a structure that represents a state flux in a hydrological model. It encapsulates 
the relationship between input fluxes and a state variable. The structure stores the mathematical
expression that defines how the state changes based on inputs and parameters.

The structure can be constructed in three ways:
1. By explicitly providing inputs, state variable, parameters, and the state expression
2. By providing input/output fluxes as a pair along with the state variable - automatically constructs
   the state expression as the difference between sum of input fluxes and sum of output fluxes
3. By providing a pair of state variables - constructs a simple state transition expression

If no name is provided, a unique name is generated based on the state variable.

This structure is particularly useful in building complex hydrological models where state 
variables evolve over time based on various input fluxes and parameters.
"""
struct StateFlux{N, M <: ComponentVector} <: AbstractStateFlux
    "flux expressions to descripe the formula of the state variable"
    exprs::Vector{Num}
    "bucket information: keys contains: input, output, param, state"
    meta::M

    function StateFlux(
        inputs::Vector{T},
        state::T,
        params::Vector{T} = T[];
        expr::T,
        name::Union{Symbol, Nothing} = nothing,
    ) where {T <: Num}
        #* Convert to a symbol based on the variable
        meta = ComponentVector(inputs = inputs, states = [state], params = params)
        flux_name = isnothing(name) ? Symbol("##state_flux#", meta) : name
        return new{flux_name, typeof(meta)}(
            [expr],
            meta,
        )
    end
    #* construct state flux with input fluxes and output fluxes
    StateFlux(fluxes::Pair{Vector{Num}, Vector{Num}}, state::Num) = StateFlux(vcat(fluxes[1], fluxes[2]), state, expr = sum(fluxes[1]) - sum(fluxes[2]))
    #* construct state flux with state variables
    StateFlux(states::Pair{Num, Num}) = StateFlux([states[2]], states[1], expr = states[2] - states[1])
end

"""
    (flux::AbstractStateFlux)(input::Union{Vector,Matrix,Array}, pas::ComponentVector; config::NamedTuple=NamedTuple(), kwargs...)

Apply the state flux model to input data of various dimensions.

# Arguments
- `input`: Input data, which can be:
  - `Vector`: A vector of input values for a single time step.
  - `Matrix`: A matrix of input values, where each column represents a different time step.
  - `Array`: A 3D array of input values, with dimensions (var_names, node_names, ts_len).
- `pas::ComponentVector`: A component vector containing parameter values.
- `config::NamedTuple`: Configuration options for the flux calculation:
  - `ptyidx`: Indices specifying which parameter types to use (only for 3D array input)
- `kwargs...`: Additional keyword arguments (unused in this function), provided for compatibility with the component callable function API.

# Returns
This function does not actually return a value, as state flux models cannot be run directly.

# Notes
- State flux models cannot be run directly and will throw an error if attempted.
- To use state flux models, they should be incorporated into a HydroFlux or other composite flux model.
"""
(::AbstractStateFlux)(::AbstractVector, ::ComponentVector; kwargs...) = @error "State Flux cannot run directly, please using HydroFlux to run"
(::AbstractStateFlux)(::AbstractArray{T, 2}, ::ComponentVector; kwargs...) where {T} = @error "State Flux cannot run directly, please using HydroFlux to run"
(::AbstractStateFlux)(::AbstractArray{T, 3}, ::ComponentVector; kwargs...) where {T} = @error "State Flux cannot run directly, please using HydroFlux to run"

