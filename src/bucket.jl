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
struct HydroBucket{S} <: AbstractBucket
    "Name of the bucket"
    name::Symbol
    "Vector of flux functions describing hydrological processes."
    fluxes::Vector{<:AbstractHydroFlux}
    "Vector of state derivative functions for ODE calculations."
    dfluxes::Vector{<:AbstractStateFlux}
    "Generated function for calculating all hydrological fluxes."
    flux_func::Function
    "Generated function for calculating all hydrological fluxes."
    multi_flux_func::Function
    "Generated function for ordinary differential equations (ODE) calculations, or nothing if no ODE calculations are needed."
    ode_func::Union{Nothing,Function}
    "Generated function for ordinary differential equations (ODE) calculations, or nothing if no ODE calculations are needed."
    multi_ode_func::Union{Nothing,Function}
    "Metadata about the bucket, including input, output, state, parameter, and neural network names."
    meta::ComponentVector

    function HydroBucket(;
        name::Union{Symbol,Nothing}=nothing,
        fluxes::Vector{<:AbstractHydroFlux},
        dfluxes::Vector{<:AbstractStateFlux}=StateFlux[],
        sort_fluxes::Bool=false,
    )
        #* sort the fluxes if needed
        fluxes = sort_fluxes ? sort_fluxes(fluxes) : fluxes
        #* Extract all variable names of fluxes and dfluxes
        bucket_inputs, bucket_outputs, bucket_states = get_all_vars(vcat(fluxes, dfluxes))
        bucket_params = reduce(union, get_param_vars.(vcat(fluxes, dfluxes)))
        bucket_nns_ntp = reduce(merge, map(flux -> NamedTuple(get_nn_vars(flux)), fluxes))
        #* Setup the meta data of the bucket
        meta = ComponentVector(inputs=bucket_inputs, outputs=bucket_outputs, states=bucket_states, params=bucket_params, nns=bucket_nns_ntp)
        #* Construct a function for ordinary differential calculation based on dfunc and funcs
        flux_func, multi_flux_func, ode_func, multi_ode_func = build_ele_func(fluxes, dfluxes, meta)
        bucket_name = isnothing(name) ? Symbol("##bucket#", hash(meta)) : name
        return new{!isempty(bucket_states)}(bucket_name, fluxes, dfluxes, flux_func, multi_flux_func, ode_func, multi_ode_func, meta)
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
function (ele::HydroBucket{true})(
    input::AbstractArray{T,2},
    pas::ComponentVector;
    config::NamedTuple=NamedTuple(),
    kwargs...,
) where {T}
    #* get kwargs
    solver = get(config, :solver, ManualSolver{true}())
    interp = get(config, :interp, DataInterpolations.LinearInterpolation)
    timeidx = get(config, :timeidx, collect(1:size(input, 2)))
    #* solve ode functions
    itpfuncs = interp(input, timeidx)
    solved_states = solver(
        (u, p, t) -> ele.ode_func(itpfuncs(t), u, p),
        pas, Vector(view(pas, :initstates)), timeidx
    )
    #* concatenate states and fluxes 
    flux_output = ele.flux_func(eachslice(input, dims=1), eachslice(solved_states, dims=1), pas)
    vcat(solved_states, permutedims(reduce(hcat, flux_output)))
end

(ele::HydroBucket{false})(input::AbstractArray{T,2}, pas::ComponentVector; kwargs...) where {T} = begin
    flux_output = ele.flux_func(input, nothing, pas)
    permutedims(reduce(hcat, flux_output))
end

function (ele::HydroBucket{true})(
    input::AbstractArray{T,3},
    pas::ComponentVector;
    config::NamedTuple=NamedTuple(),
    kwargs...,
) where {T}
    input_dims, num_nodes, time_len = size(input)

    #* get kwargs
    ptyidx = get(config, :ptyidx, 1:size(input, 2))
    styidx = get(config, :styidx, 1:size(input, 2))
    interp = get(config, :interp, LinearInterpolation)
    solver = get(config, :solver, ManualSolver{true}())
    timeidx = get(config, :timeidx, collect(1:size(input, 3)))

    #* prepare states parameters and nns
    params = view(pas, :params)
    nn_params = isempty(get_nn_vars(ele)) ? ones(eltype(pas), num_nodes) : view(pas, :nns)
    expand_params = ComponentVector(NamedTuple{Tuple(get_param_names(ele))}([params[p][ptyidx] for p in get_param_names(ele)]))
    new_pas = ComponentVector(params=expand_params, nns=nn_params)
    initstates_mat = view(reshape(Vector(view(pas, :initstates)), num_nodes, :)', :, styidx)

    #* prepare input function
    input_reshape = reshape(input, input_dims * num_nodes, time_len)
    itpfuncs = interp(input_reshape, timeidx)
    solved_states = solver(
        (u, p, t) -> begin
            tmp_input = reshape(itpfuncs(t), input_dims, num_nodes)
            reduce(hcat, ele.multi_ode_func(eachslice(tmp_input, dims=1), eachslice(u, dims=1), p)) |> permutedims
        end,
        new_pas, initstates_mat, timeidx
    )
    #* run other functions
    output = ele.multi_flux_func(eachslice(input, dims=1), eachslice(solved_states, dims=1), new_pas)
    output_arr = length(output) > 1 ? permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), output), (3, 1, 2)) : reshape(output[1], 1, num_nodes, time_len)
    cat(solved_states, output_arr, dims=1)
end

function (ele::HydroBucket{false})(
    input::AbstractArray{T,3},
    pas::ComponentVector;
    config::NamedTuple=NamedTuple(),
    kwargs...,
) where {T}
    ptyidx = get(config, :ptyidx, 1:size(input, 2))
    params = view(pas, :params)
    nn_params = isempty(get_nn_vars(ele)) ? ones(eltype(pas), size(input, 2)) : view(pas, :nns)
    expand_params = ComponentVector(NamedTuple{Tuple(get_param_names(ele))}([params[p][ptyidx] for p in get_param_names(ele)]))
    new_pas = ComponentVector(params=expand_params, nns=nn_params)
    #* run other functions
    output = ele.flux_func(eachslice(input, dims=1), nothing, new_pas)
    permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), output), (3, 1, 2))
end
