"""
    HydroBucket{N,S} <: AbstractBucket

Represents a hydrological bucket model component that handles both single-node and multi-node 
computations. The type parameters are:
- `N`: Component name (encoded at type level for better type stability)
- `S`: Boolean indicating if the bucket has state variables

# Arguments
- `fluxes::Vector{<:AbstractHydroFlux}`: Vector of flux functions
- `dfluxes::Vector{<:AbstractStateFlux}=StateFlux[]`: Vector of state derivative functions
- `name::Union{Symbol,Nothing}=nothing`: Optional bucket identifier
- `sort_fluxes::Bool=false`: Whether to sort fluxes by dependencies

# Fields
- `fluxes::Vector{<:AbstractHydroFlux}`: Vector of flux functions
- `dfluxes::Vector{<:AbstractStateFlux}`: Vector of state derivative functions
- `flux_funcs::Vector{<:Function}`: Generated functions for flux calculations
- `ode_funcs::Union{Nothing,Vector}`: Generated functions for ODE calculations
- `infos::NamedTuple`: Component metadata

# Description
HydroBucket is a type-stable implementation of a hydrological bucket model that supports both 
single-node and distributed (multi-node) computations. It automatically generates optimized 
functions for flux and ODE calculations based on the provided process functions.

## Model Structure
The bucket model consists of:
- Process functions (`fluxes`): Define water movement between storages
- State derivatives (`dfluxes`): Define rate of change for state variables
- Generated functions: Optimized implementations for calculations

## Metadata Components
The `infos` field tracks:
- `inputs`: External forcing variables (e.g., precipitation)
- `outputs`: Model-generated variables (e.g., runoff)
- `states`: Internal storage variables (e.g., soil moisture)
- `params`: Model parameters
- `nns`: Neural network components (if any)

## Usage Notes
1. For single-node simulations: Use `flux_funcs[1]` and `ode_funcs[1]`
2. For multi-node simulations: Use `flux_funcs[2]` and `ode_funcs[2]`
3. Parameters should be provided as ComponentVector for type stability
4. Broadcasting operations are automatically handled for multi-node cases
"""
struct HydroBucket{N,S} <: AbstractBucket
    "Vector of flux functions describing hydrological processes."
    fluxes::Vector{<:AbstractHydroFlux}
    "Vector of state derivative functions for ODE calculations."
    dfluxes::Vector{<:AbstractStateFlux}
    "Generated function for calculating all hydrological fluxes. (Supports single-node data, multi-nodes data)"
    flux_funcs::Vector{<:Function}
    "Generated function for ordinary differential equations (ODE) calculations, or nothing if no ODE calculations are needed. (Supports single-node data)"
    ode_funcs::Union{Nothing,Vector}
    "Metadata about the bucket, including input, output, state, parameter, and neural network names."
    infos::NamedTuple

    function HydroBucket(;
        name::Union{Symbol,Nothing}=nothing,
        fluxes::Vector{<:AbstractHydroFlux},
        dfluxes::Vector{<:AbstractStateFlux}=StateFlux[],
        sort_fluxes::Bool=false,
    )
        #* sort the fluxes if needed
        fluxes = sort_fluxes ? sort_fluxes(fluxes) : fluxes
        #* Extract all variable names of fluxes and dfluxes
        inputs, outputs, states = get_vars(fluxes, dfluxes)
        params = reduce(union, get_params.(vcat(fluxes, dfluxes)))
        nns = reduce(union, get_nns.(fluxes))
        infos = (; inputs=inputs, outputs=outputs, states=states, params=params, nns=nns)
        #* Construct a function for ordinary differential calculation based on dfunc and funcs
        flux_funcs, ode_funcs = build_ele_func(fluxes, dfluxes, infos)
        bucket_name = isnothing(name) ? Symbol("##bucket#", hash(infos)) : name
        return new{bucket_name,!isempty(states)}(fluxes, dfluxes, flux_funcs, ode_funcs, infos)
    end
end

"""
    (bucket::HydroBucket{N,S})(input::AbstractArray{T,2}, params::ComponentVector; kwargs...)
    (bucket::HydroBucket{N,S})(input::AbstractArray{T,3}, params::ComponentVector; kwargs...)

Run the HydroBucket model for given input and parameters.

# Arguments
- `input`: Input data array with dimensions:
  - For 2D input: variables × time (single node)
  - For 3D input: variables × nodes × time (multi-node)
- `params`: ComponentVector containing model parameters and initial states
- `kwargs`: Optional keyword arguments:
  - `solver`: ODE solver to use (default: ManualSolver{true}())
  - `interp`: Interpolation method (default: LinearInterpolation)
  - `timeidx`: Time indices (default: 1:size(input, last_dim))
  - `ptyidx`: Parameter indices for multi-node runs
  - `styidx`: State indices for multi-node runs
  - `initstates`: Initial states for ODE solving

# Returns
Array containing model outputs:
- For 2D input: Array of size (states+outputs) × time
- For 3D input: Array of size (states+outputs) × nodes × time

# Description
This function executes the bucket model by:
1. Setting up interpolation for input time series
2. Preparing parameters and initial states
3. Solving ODEs if the model has state variables (S=true)
4. Computing fluxes using the model's flux functions
5. Combining states and fluxes into the output array

## Single Node Operation (2D input)
- Processes one location's time series
- Uses provided parameters and initial states
- Solves ODEs if model has states
- Returns combined states and fluxes

## Multi-Node Operation (3D input)
- Processes multiple locations simultaneously
- Supports both shared and distributed parameters
- Handles state evolution for each node
- Returns results for all nodes

Note: Input dimensions must match the number of input variables defined in the model's
metadata. All required parameters and initial states must be present in the params argument.
"""
function (ele::HydroBucket{N,true})(input::AbstractArray{T,2}, params::ComponentVector; kwargs...) where {T,N}
    #* get kwargs
    solver = get(kwargs, :solver, ManualSolver(mutable=true))
    interp = get(kwargs, :interp, DirectInterpolation)
    timeidx = get(kwargs, :timeidx, collect(1:size(input, 2)))
    #* prepare initstates
    initstates = get(kwargs, :initstates, zeros(eltype(params), length(get_state_names(ele))))
    initstates_ = initstates isa ComponentVector ? Vector(initstates[get_state_names(ele)]) : initstates
    #* get params axes
    param_vec, params_axes = Vector(params), getaxes(params)
    #* solve ode functions
    itpfuncs = interp(input, timeidx)
    solved_states = solver(
        (u, p, t) -> ele.ode_funcs[1](itpfuncs(t), u, ComponentVector(p, params_axes)),
        param_vec, initstates_, timeidx
    )
    #* concatenate states and fluxes 
    flux_output = ele.flux_funcs[1](eachslice(input, dims=1), eachslice(solved_states, dims=1), params)
    vcat(solved_states, stack(flux_output, dims=1))
end

function (ele::HydroBucket{N,true})(input::AbstractArray{T,3}, params::ComponentVector; kwargs...) where {T,N}
    input_dims, num_nodes, time_len = size(input)

    #* get kwargs
    ptyidx = get(kwargs, :ptyidx, collect(1:num_nodes))
    styidx = get(kwargs, :styidx, collect(1:num_nodes))
    device = get(kwargs, :device, identity)
    solver = get(kwargs, :solver, ManualSolver(mutable=true))
    interp = get(kwargs, :interp, DirectInterpolation)
    timeidx = get(kwargs, :timeidx, collect(1:size(input, 3)))

    #* prepare initstates
    initstates = get(kwargs, :initstates, zeros(eltype(params), length(get_state_names(ele)), num_nodes)) |> device
    initstates_ = initstates isa ComponentVector ? initstates[get_state_names(ele)] : initstates
    initstates_mat = expand_component_initstates(initstates_, styidx) |> device

    #* prepare states parameters and nns
    new_params = expand_component_params(params, ptyidx) |> device
    params_vec, params_axes = Vector(new_params) |> device, getaxes(new_params)

    #* prepare input function
    itpfuncs = interp.(eachslice(input, dims=1), Ref(timeidx))

    solved_states = solver(
        (u, p, t) -> ele.ode_funcs[2](
            ntuple(i -> itpfuncs[i](t), length(itpfuncs)),
            eachslice(u, dims=1),
            ComponentVector(p, params_axes)
        ),
        params_vec, initstates_mat, timeidx
    )

    #* run other functions
    output = ele.flux_funcs[2](eachslice(input, dims=1), eachslice(solved_states, dims=1), new_params)
    cat(solved_states, stack(output, dims=1), dims=1)
end

(ele::HydroBucket{N,false})(input::AbstractArray{T,2}, params::ComponentVector; kwargs...) where {T,N} = begin
    stack(ele.flux_funcs[1](eachslice(input, dims=1), nothing, params), dims=1)
end

function (ele::HydroBucket{N,false})(input::AbstractArray{T,3}, params::ComponentVector; kwargs...) where {T,N}
    ptyidx = get(kwargs, :ptyidx, 1:size(input, 2))
    output = ele.flux_funcs[2](eachslice(input, dims=1), nothing, expand_component_params(params, ptyidx))
    stack(output, dims=1)
end