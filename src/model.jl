"""
    HydroModel{N} <: AbstractModel

Represents a complete hydrological model that integrates multiple components for simulating 
water processes. The type parameter `N` encodes the model name at type level for better 
type stability.

# Arguments
- `name::Symbol`: Model identifier
- `components::Vector{<:AbstractComponent}`: Vector of model components
- `solver::Union{Nothing,ODESolver}=nothing`: Optional ODE solver for state equations
- `interp::Union{Nothing,Function}=nothing`: Optional interpolation method for inputs

# Fields
- `components::Vector{<:AbstractComponent}`: Hydrological computation elements
- `varindices::Vector{Vector{Int}}`: Input variable indices for each component
- `infos::NamedTuple`: Model metadata including:
  - `name`: Model identifier
  - `inputs`: External input variables
  - `outputs`: Model output variables
  - `states`: State variables
  - `allvars`: All model variables
  - `nns`: Neural network components

# Description
HydroModel provides a framework for building and running hydrological simulations by 
connecting multiple components in a type-stable and efficient manner.

## Model Structure
The model consists of:
- Components: Individual process modules (e.g., buckets, fluxes, routing)
- Connections: Automatic variable mapping between components
- States: Time-evolving system variables
- Parameters: Model coefficients and settings

## Component Integration
Components are automatically connected based on their input/output interfaces:
1. Variable matching between components
2. Execution order determination
3. Data flow management
4. State variable handling

## Usage Notes
1. Model Construction:
   ```julia
   model = HydroModel(:lumped_model,
       components=[
           HydroBucket(:bucket1, fluxes=[rainfall_flux, evap_flux]),
           UnitHydrograph(:routing1, response=[0.3, 0.5, 0.2])
       ]
   )
   ```

2. Model Simulation:
   ```julia
   # Single-node run
   output = model(input, params)  # input: (variables, timesteps)

   # Multi-node run
   output = model(input, params)  # input: (variables, nodes, timesteps)
   ```

3. Configuration Options:
   - `solver`: ODE solver for state equations (e.g., Tsit5())
   - `interp`: Input interpolation method
   - Component-specific options via kwargs

# Notes
- Components are executed in sequence based on dependencies
- State variables are integrated automatically when using an ODE solver
- Input interpolation is applied when timesteps don't align
- Type stability is maintained throughout the computation
"""
struct HydroModel{N} <: AbstractModel
    "hydrological computation elements"
    components::Vector{<:AbstractComponent}
    "meta data of hydrological model"
    infos::NamedTuple
    "input variables index for each components"
    _varindices::AbstractVector{AbstractVector{Int}}
    "output variables index for sort output variables"
    _outputindices::AbstractVector{Int}

    function HydroModel(;
        components::Vector{C},
        name::Union{Symbol,Nothing}=nothing,
        sort_components::Bool=false,
    ) where {C<:AbstractComponent}
        components = sort_components ? sort_components(components) : components
        input_names, output_names, state_names = get_var_names(components)
        param_names = reduce(union, get_param_names.(components))
        nn_names = reduce(union, get_nn_names.(components))
        input_idx, output_idx = _prepare_indices(components, input_names, vcat(state_names, output_names))
        infos = (; inputs=input_names, outputs=output_names, states=state_names, params=param_names, nns=nn_names)
        model_name = isnothing(name) ? Symbol("##model#", hash(infos)) : name
        new{model_name}(components, infos, input_idx, output_idx)
    end
end


"""
    _prepare_indices(components::Vector{<:AbstractComponent}, 
                    input_names::Vector{Symbol}, 
                    vcat_names::Vector{Symbol}) -> Tuple{Vector{Vector{Int}}, Vector{Int}}

Prepare input and output variable indices for connecting components in a hydrological model.

# Arguments
- `components::Vector{<:AbstractComponent}`: Vector of model components to process
- `input_names::Vector{Symbol}`: Initial list of input variable names
- `vcat_names::Vector{Symbol}`: Concatenated list of all variable names (inputs, states, outputs)

# Returns
- `Tuple{Vector{Vector{Int}}, Vector{Int}}`: A tuple containing:
  - First element: Vector of input indices for each component
  - Second element: Vector of output indices for all components

# Description
This internal function manages the variable connections between components by:
1. Mapping each component's input variables to global input indices
2. Tracking state and output variables as they become available
3. Building the complete variable dependency chain

## Process Flow
1. For each component:
   - Map its input variables to current input name indices
   - Add its state and output variables to available inputs
   - Map its outputs to the complete variable list

## Implementation Details
- Uses `findfirst` for efficient index mapping
- Maintains order of variable declarations
- Handles both direct and derived variables
- Ensures proper variable propagation

# Notes
- This is an internal function used by HydroModel constructor
- Variable names must be unique across all components
- The order of components affects the variable propagation chain
"""
function _prepare_indices(components::Vector{<:AbstractComponent}, input_names::Vector{Symbol}, vcat_names::Vector{Symbol})
    input_idx, output_idx = Vector{Int}[], Vector{Int}()
    for component in components
        #* extract input index
        tmp_input_idx = map((nm) -> findfirst(varnm -> varnm == nm, input_names), get_input_names(component))
        push!(input_idx, tmp_input_idx)
        #* extract output index
        tmp_cpt_vcat_names = vcat(get_state_names(component), get_output_names(component))
        input_names = vcat(input_names, tmp_cpt_vcat_names)
    end
    for name in vcat_names
        push!(output_idx, findfirst(varnm -> varnm == name, input_names))
    end
    return input_idx, output_idx
end

"""
    (model::HydroModel)(input::AbstractArray, params::ComponentVector; kwargs...)

Run a hydrological model simulation with the given input data and parameters.

# Arguments
- `input::AbstractArray`: Input data with dimensions:
  - `Matrix`: Shape (variables, timesteps) for single-node simulation
  - `Array{3}`: Shape (variables, nodes, timesteps) for distributed simulation
- `params::ComponentVector`: Model parameters for all components
- `initstates::ComponentVector=ComponentVector()`: Initial states for components
- `config::Union{Dict,Vector{<:Dict}}=Dict()`: Configuration for components
- `kwargs...`: Additional keyword arguments passed to components

# Returns
- `Matrix`: For 2D input, returns matrix of shape (output_variables, timesteps)
- `Array{3}`: For 3D input, returns array of shape (output_variables, nodes, timesteps)

# Description
This function runs the hydrological model simulation by:
1. Processing input data through each component sequentially
2. Managing state variables and their evolution
3. Collecting and organizing output variables
4. Handling both single-node and distributed simulations

## Simulation Process
1. Initialize:
   - Set up component configurations
   - Prepare initial states (default or user-provided)
   - Organize input data structure

2. Component Execution:
   - Pass relevant inputs to each component
   - Update states based on component calculations
   - Collect outputs and prepare for next component

3. Output Processing:
   - Gather all component outputs
   - Extract requested output variables
   - Maintain proper dimensionality

## Configuration Options
- Component-specific settings via `config` dict
- Initial state values via `initstates`
- Solver settings for ODE components
- Interpolation methods for input data

# Example
```julia
# Single-node simulation
model = HydroModel(:lumped_model, components=[bucket, routing])
output = model(
    input,                     # shape: (variables, timesteps)
    params,                    # component parameters
    initstates=init_states,    # initial conditions
    config=Dict(              # component settings
        :solver => ManualSolver{true}(),
        :interp => LinearInterpolation
    )
)

# Distributed simulation
output = model(
    input,                     # shape: (variables, nodes, timesteps)
    params,                    # parameters per node
    config=[                   # settings per component
        Dict(                 
            :solver => ManualSolver{true}(),
            :interp => LinearInterpolation
        ),
        Dict(                  
            :solver => ManualSolver{true}(),
            :interp => LinearInterpolation
        ),
        ... # number is consistent with components
    ]
)
```

# Notes
- Input variables must match component requirements
- Parameters must be properly structured in ComponentVector
- State initialization is automatic if not provided
- Configuration can be shared or component-specific
"""
function (model::HydroModel{N})(
    input::AbstractArray{T,2},
    params::ComponentVector;
    kwargs...,
) where {N,T<:Number}
    initstates = get(kwargs, :initstates, get_default_states(model, eltype(input)))
    config = get(kwargs, :config, NamedTuple())
    comp_configs = config isa NamedTuple ? fill(config, length(model.components)) : config
    @assert length(comp_configs) == length(model.components) "component configs length must be equal to components length"
    @assert size(input, 1) == length(get_input_names(model)) "input variables length must be equal to input variables length"
    outputs = input
    for (idx_, comp_, config_) in zip(model._varindices, model.components, comp_configs)
        tmp_outputs = comp_(outputs[idx_, :], params; initstates=initstates[get_state_names(comp_)], config_...)
        outputs = cat(outputs, tmp_outputs, dims=1)
    end
    return outputs[model._outputindices, :]
end

function (model::HydroModel{N})(
    input::AbstractArray{T,3},
    params::ComponentVector;
    kwargs...,
) where {N,T<:Number}
    config = get(kwargs, :config, NamedTuple())
    initstates = get(kwargs, :initstates, get_default_states(model, size(input, 2), eltype(input)))
    comp_configs = config isa NamedTuple ? fill(config, length(model.components)) : config
    @assert length(comp_configs) == length(model.components) "component configs length must be equal to components length"
    @assert size(input, 1) == length(get_input_names(model)) "input variables length must be equal to input variables length"
    outputs = input
    for (idx_, comp_, config_) in zip(model._varindices, model.components, comp_configs)
        tmp_outputs = comp_(outputs[idx_, :, :], params; initstates=initstates[get_state_names(comp_)], config_...)
        outputs = cat(outputs, tmp_outputs, dims=1)
    end
    return outputs[model._outputindices, :, :]
end