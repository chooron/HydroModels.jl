"""
    HydroBucket(; 
        funcs::Vector{<:AbstractHydroFlux}, 
        dfuncs::Vector{<:AbstractStateFlux}=StateFlux[], 
        name::Union{Symbol,Nothing}=nothing, 
        sort_funcs::Bool=false
    )

Represents a hydrological bucket model component that handles both single-node and multi-node computations.

# Arguments
- `funcs::Vector{<:AbstractHydroFlux}`: Vector of flux functions describing hydrological processes
- `dfuncs::Vector{<:AbstractStateFlux}`: Vector of state derivative functions for ODE calculations (default: empty)
- `name::Union{Symbol,Nothing}`: Optional bucket identifier. Defaults to auto-generated name from state variables
- `sort_funcs::Bool}`: Whether to topologically sort flux functions based on dependencies (default: false)

# Fields
- `fluxes::Vector{<:AbstractHydroFlux}`: Vector of flux functions describing hydrological processes
- `dfluxes::Vector{<:AbstractStateFlux}`: Vector of state derivative functions for ODE calculations
- `flux_func::Function`: Generated function for single-node flux calculations
- `multi_flux_func::Function`: Generated function for multi-node parallel flux calculations
- `ode_func::Union{Nothing,Function}`: Generated function for single-node ODE calculations
- `multi_ode_func::Union{Nothing,Function}`: Generated function for multi-node parallel ODE calculations
- `meta::ComponentVector`: Metadata containing model structure information

# Description
HydroBucket is a type-stable implementation of a hydrological bucket model that supports both 
single-node and distributed (multi-node) computations. It automatically generates optimized 
functions for flux and ODE calculations based on the provided process functions.

## Model Structure
The bucket model consists of:
- Process functions (`fluxes`): Define water movement between storages
- State derivatives (`dfluxes`): Define rate of change for state variables
- Generated functions: Optimized implementations for both single and multi-node calculations

## Metadata Components
The `meta` field tracks:
- `inputs`: External forcing variables (e.g., precipitation, temperature)
- `outputs`: Model-generated variables (e.g., runoff, evaporation)
- `states`: Internal storage variables (e.g., soil moisture, groundwater)
- `params`: Model parameters controlling process behavior
- `nn_vars`: Neural network components (if any) for hybrid modeling

## Performance Features
- Type-stable computations for both single and multi-node cases
- Efficient broadcasting operations for vectorized calculations
- Automatic function generation with optimized broadcasting
- Support for ComponentArray parameters for structured data handling

## Usage Notes
1. For single-node simulations: Use `flux_func` and `ode_func`
2. For multi-node simulations: Use `multi_flux_func` and `multi_ode_func`
3. Parameters should be provided as ComponentVector for type stability
4. Broadcasting operations are automatically handled for multi-node cases

See also: [`AbstractHydroFlux`](@ref), [`AbstractStateFlux`](@ref), [`ComponentVector`](@ref)
"""
struct HydroBucket{S} <: AbstractBucket
    "Name of the bucket"
    name::Symbol
    "Vector of flux functions describing hydrological processes."
    fluxes::Vector{<:AbstractHydroFlux}
    "Vector of state derivative functions for ODE calculations."
    dfluxes::Vector{<:AbstractStateFlux}
    "Generated function for calculating all hydrological fluxes. (Supports single-node data)"
    flux_func::Function
    "Generated function for calculating all hydrological fluxes. (Supports multi-nodes data)"
    multi_flux_func::Function
    "Generated function for ordinary differential equations (ODE) calculations, or nothing if no ODE calculations are needed. (Supports single-node data)"
    ode_func::Union{Nothing,Function}
    "Generated function for ordinary differential equations (ODE) calculations, or nothing if no ODE calculations are needed. (Supports multi-nodes data)"
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
    input::AbstractArray{T,2}, pas::ComponentVector;
    config::NamedTuple=NamedTuple(), kwargs...
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
    flux_output = ele.flux_func(eachslice(input, dims=1), nothing, pas)
    permutedims(reduce(hcat, flux_output))
end

function (ele::HydroBucket{true})(
    input::AbstractArray{T,3}, pas::ComponentVector;
    config::NamedTuple=NamedTuple(), kwargs...
) where {T}
    input_dims, num_nodes, time_len = size(input)

    #* get kwargs
    ptyidx = get(config, :ptyidx, 1:size(input, 2))
    styidx = get(config, :styidx, 1:size(input, 2))
    interp = get(config, :interp, LinearInterpolation)
    solver = get(config, :solver, ManualSolver{true}())
    timeidx = get(config, :timeidx, collect(1:size(input, 3)))

    #* prepare states parameters and nns
    new_pas = expand_component_params(pas, ptyidx)
    initstates_mat = expand_component_initstates(pas, styidx)

    #* prepare input function
    itpfuncs = interp(reshape(input, input_dims * num_nodes, time_len), timeidx)
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
    input::AbstractArray{T,3}, pas::ComponentVector;
    config::NamedTuple=NamedTuple(), kwargs...,
) where {T}
    ptyidx = get(config, :ptyidx, 1:size(input, 2))
    new_pas = expand_component_params(pas, ptyidx)
    #* run other functions
    output = ele.flux_func(eachslice(input, dims=1), nothing, new_pas)
    permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), output), (3, 1, 2))
end
