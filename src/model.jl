"""
Model module - defines complete hydrological model, combining multiple components.
"""

"""
    HydroModel{CS, NT, VI, OI, DS} <: AbstractModel

Represents a complete hydrological model composed of a sequence of components (e.g., HydroBucket, HydroRoute).

It automatically manages the data flow between components, routing the outputs of one component as inputs to the next.

$(FIELDS)

# Type Parameters
- `CS`: Components tuple type
- `NT`: Metadata type
- `VI`: Variable indices type
- `OI`: Output indices type
- `DS`: Default states type
"""
struct HydroModel{CS,NT,VI,OI,DS} <: AbstractModel
    "hydrological model name"
    name::Symbol
    "hydrological computation elements"
    components::CS
    "meta data of hydrological model"
    infos::NT
    "input variables index for each components"
    _varindices::VI
    "output variables index for sort output variables"
    _outputindices::OI
    "default initial states"
    _defaultstates::DS
    
    function HydroModel(;
        components::Tuple,
        name::Optional{Symbol}=nothing,
    )
        inputs, outputs, states = get_var_names(components)
        params = reduce(union, get_param_names.(components); init=Symbol[])
        nns = reduce(union, get_nn_names.(components); init=Symbol[])
        
        input_idx, output_idx = _prepare_indices(components, inputs, vcat(states, outputs))
        
        infos = HydroInfos(
            inputs=inputs, states=states, outputs=outputs,
            params=params, nns=nns
        )
        
        model_name = isnothing(name) ? Symbol("##model#", hash(infos)) : name
        states_axes = Axis(NamedTuple{Tuple(states)}(1:length(states)))
        default_states = ComponentVector(zeros(length(states)), states_axes)
        
        new{typeof(components),typeof(infos),typeof(input_idx),typeof(output_idx),typeof(default_states)}(
            model_name, components, infos, input_idx, output_idx, default_states
        )
    end
end

"""
    _prepare_indices(components, input_names, var_names)

Prepare input and output indices for the model.

# Arguments
- `components`: Tuple of components
- `input_names`: Input variable names
- `var_names`: All variable names (states + outputs)

# Returns
- `input_idx`: Vector of input variable indices for each component
- `output_idx`: Vector of final output variable indices
"""
function _prepare_indices(
    components::CT,
    input_names::Vector{Symbol},
    var_names::Vector{Symbol}
) where {CT}
    input_idx = Vector{Int}[]
    output_idx = Int[]
    
    # Accumulate available variable names
    available_vars = copy(input_names)
    
    for component in components
        # Find input indices for current component
        comp_input_names = get_input_names(component)
        tmp_input_idx = map(comp_input_names) do nm
            findfirst(==(nm), available_vars)
        end
        
        if nothing in tmp_input_idx
            missing_vars = comp_input_names[tmp_input_idx .=== nothing]
            @warn "Input variables $missing_vars not found in available_vars"
        end
        
        push!(input_idx, filter(!isnothing, tmp_input_idx))
        
        # Add component outputs and states to available variables
        comp_output_state = vcat(get_state_names(component), get_output_names(component))
        available_vars = vcat(available_vars, comp_output_state)
    end
    
    # Find indices for final outputs
    for name in var_names
        idx = findfirst(==(name), available_vars)
        if !isnothing(idx)
            push!(output_idx, idx)
        end
    end
    
    return input_idx, output_idx
end

"""
    @hydromodel [name] begin ... end

Macro to conveniently create a HydroModel from a sequence of components.

# Usage
The macro takes an optional name and a begin...end block containing an ordered list of component instances.

```jldoctest
julia> # Assuming bucket1 and route1 are pre-defined components
julia> @hydromodel :my_full_model begin
           bucket1
           route1
       end
```
"""
macro hydromodel(args...)
    name = length(args) == 1 ? nothing : args[1]
    expr = length(args) == 1 ? args[1] : args[2]
    
    @assert Meta.isexpr(expr, :block) "Expected a begin...end block after model name"
    
    components = filter(x -> !(x isa LineNumberNode), expr.args)
    
    return esc(:(HydroModel(; name=$(name), components=tuple($(components...)))))
end

"""
    (model::HydroModel)(input, params, config; kwargs...)

Execute the full hydrological model simulation. This is the functor implementation for HydroModel.

It sequentially runs each component in the model, automatically handling the routing of data between them. 
It supports both 2D (single-node) and 3D (multi-node) inputs.

Common kwargs include `initstates` and `config` (for component-specific settings).

# Arguments
- `input`: Input array with dimensions (variables, timesteps) or (variables, nodes, timesteps)
- `params`: ComponentVector of parameters
- `config`: Configuration object or configuration tuple (for multi-component)
- `kwargs`: Additional keyword arguments

# Returns
- Output array containing all state and output variables
"""
function (model::HydroModel)(
    input::AbstractArray{T,D},
    params::ComponentVector,
    config::Union{ConfigType,Tuple}=default_config();
    kwargs...
) where {T,D}
    # Normalize configuration
    comp_configs = if config isa Tuple || config isa NamedTuple{(:components,)}
        # Multi-component configuration
        length(config) == length(model.components) || 
            error("Component configs length must equal components length")
        config
    else
        # Single configuration, apply to all components
        ntuple(_ -> normalize_config(config), length(model.components))
    end
    
    @assert D in (2, 3) "Input array dimension must be 2 or 3"
    @assert size(input, 1) == length(get_input_names(model)) "Input variables length must equal model input variables length"
    
    # Get initial states
    initstates = get(kwargs, :initstates, model._defaultstates)
    
    # Validate initstates contains required state names
    required_states = get_state_names(model)
    if !isempty(required_states)
        for state in required_states
            if !hasproperty(initstates, state)
                throw(ArgumentError("Missing required state: $state in initstates"))
            end
        end
    end
    
    # Initialize output as input
    outputs = input
    
    # Execute each component sequentially
    for (idx_, comp_, config_) in zip(model._varindices, model.components, comp_configs)
        # Use view to avoid copying, improving performance
        tmp_input = @view outputs[idx_, ntuple(_ -> Colon(), D - 1)...]
        
        # Get component initial states
        comp_state_names = get_state_names(comp_)
        comp_initstates = if !isempty(comp_state_names)
            initstates[comp_state_names]
        else
            nothing
        end
        
        # Execute component
        tmp_output = if !isnothing(comp_initstates)
            comp_(tmp_input, params, config_; initstates=comp_initstates)
        else
            comp_(tmp_input, params, config_)
        end
        
        # Concatenate output (avoiding in-place modification)
        outputs = vcat(outputs, tmp_output)
    end
    
    # Return final output
    return @view outputs[model._outputindices, ntuple(_ -> Colon(), D - 1)...]
end

# Specialized version for 2D input
function (model::HydroModel)(
    input::AbstractArray{T,2},
    params::ComponentVector,
    config::Union{ConfigType,Tuple}=default_config();
    kwargs...
) where {T}
    # Normalize configuration
    comp_configs = if config isa Tuple
        length(config) == length(model.components) || 
            error("Component configs length must equal components length")
        map(normalize_config, config)
    else
        ntuple(_ -> normalize_config(config), length(model.components))
    end
    
    @assert size(input, 1) == length(get_input_names(model)) "Input variables length mismatch"
    
    initstates = get(kwargs, :initstates, model._defaultstates)
    outputs = input
    
    for (idx_, comp_, config_) in zip(model._varindices, model.components, comp_configs)
        tmp_input = @view outputs[idx_, :]
        
        comp_state_names = get_state_names(comp_)
        comp_initstates = !isempty(comp_state_names) ? initstates[comp_state_names] : nothing
        
        tmp_output = if !isnothing(comp_initstates)
            comp_(tmp_input, params, config_; initstates=comp_initstates)
        else
            comp_(tmp_input, params, config_)
        end
        
        outputs = vcat(outputs, tmp_output)
    end
    
    return @view outputs[model._outputindices, :]
end

# Specialized version for 3D input
function (model::HydroModel)(
    input::AbstractArray{T,3},
    params::ComponentVector,
    config::Union{ConfigType,Tuple}=default_config();
    kwargs...
) where {T}
    # Normalize configuration
    comp_configs = if config isa Tuple
        length(config) == length(model.components) || 
            error("Component configs length must equal components length")
        map(normalize_config, config)
    else
        ntuple(_ -> normalize_config(config), length(model.components))
    end
    
    @assert size(input, 1) == length(get_input_names(model)) "Input variables length mismatch"
    
    initstates = get(kwargs, :initstates, model._defaultstates)
    outputs = input
    
    for (idx_, comp_, config_) in zip(model._varindices, model.components, comp_configs)
        tmp_input = @view outputs[idx_, :, :]
        
        comp_state_names = get_state_names(comp_)
        comp_initstates = !isempty(comp_state_names) ? initstates[comp_state_names] : nothing
        
        tmp_output = if !isnothing(comp_initstates)
            comp_(tmp_input, params, config_; initstates=comp_initstates)
        else
            comp_(tmp_input, params, config_)
        end
        
        outputs = vcat(outputs, tmp_output)
    end
    
    return @view outputs[model._outputindices, :, :]
end

# Export interfaces
export HydroModel, @hydromodel
