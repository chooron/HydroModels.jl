"""    
    check(; component, input, pas, initstates, timeidx)
    check(; component, input, pas, initstates, timeidx, ptypes, stypes)

Validate component inputs, parameters, and initial states, parameter types and state types (for 3D input).

# Arguments
- `component::AbstractComponent`: Component to validate
- `input::AbstractArray{<:Number,2}`: Input data with shape (variables, timesteps)
- `pas::ComponentVector`: Parameter collection
- `initstates::ComponentVector`: Initial states collection
- `timeidx::AbstractVector`: Time indices

- `ptypes::AbstractVector{Int}`: Parameter type indices
- `stypes::AbstractVector{Int}`: State type indices
"""
function check(component::AbstractComponent, input::AbstractArray{<:Number,2}, pas::ComponentVector, initstates::ComponentVector, timeidx::AbstractVector)
    check_input(component, input, timeidx)
    check_params(component, pas)
    check_initstates(component, initstates)
    check_nns(component, pas)
end

function check(component::AbstractComponent, input::AbstractArray{<:Number,3}, pas::ComponentVector, initstates::ComponentVector, timeidx::AbstractVector, ptypes::AbstractVector{Int}, stypes::AbstractVector{Int})
    check_input(component, input, timeidx)
    check_ptypes(component, input, ptypes)
    check_stypes(component, input, stypes)
    check_params(component, pas)
    check_initstates(component, initstates)
    check_nns(component, pas)
end

"""    
    check_input(component, input, timeidx)

Validate input data dimensions against component requirements.

# Arguments
- `component::AbstractComponent`: Component to validate
- `input::AbstractArray`: Input data (2D or 3D array)
- `timeidx::AbstractVector`: Time indices
"""
function check_input(component::AbstractComponent, input::AbstractArray, timeidx::AbstractVector)
    component_name = get_name(component)
    input_names = get_input_names(component)
    expected_vars = length(input_names)
    actual_vars = size(input, 1)
    
    @assert actual_vars == expected_vars "Input variables in component '$(component_name)' do not match required dimensions. Expected $(expected_vars) variables ($(input_names)), got $(actual_vars) variables"
    
    time_dim = ndims(input) == 2 ? 2 : 3
    expected_timesteps = length(timeidx)
    actual_timesteps = size(input, time_dim)
    
    @assert actual_timesteps == expected_timesteps "Time steps in component '$(component_name)' do not match required length. Expected $(expected_timesteps) steps, got $(actual_timesteps) steps"
end

"""    
    check_ptypes(component, input, ptypes)

Validate parameter type count matches the node count in input data.

# Arguments
- `component::AbstractComponent`: Component to validate
- `input::AbstractArray{<:Number,3}`: 3D input data
- `ptypes::AbstractVector{Int}`: Parameter type indices
"""
function check_ptypes(component::AbstractComponent, input::AbstractArray{<:Number,3}, ptypes::AbstractVector{Int})
    component_name = get_name(component)
    expected_count = size(input, 2)
    actual_count = length(ptypes)
    
    @assert actual_count == expected_count "Parameter type count mismatch in component '$(component_name)'. Expected $(expected_count) parameter types, got $(actual_count)."
end

"""    
    check_stypes(component, input, stypes)

Validate state type count matches the node count in input data.

# Arguments
- `component::AbstractComponent`: Component to validate
- `input::AbstractArray{<:Number,3}`: 3D input data
- `stypes::AbstractVector{Int}`: State type indices
"""
function check_stypes(component::AbstractComponent, input::AbstractArray{<:Number,3}, stypes::AbstractVector{Int})
    component_name = get_name(component)
    expected_count = size(input, 2)
    actual_count = length(stypes)
    
    @assert actual_count == expected_count "State type count mismatch in component '$(component_name)'. Expected $(expected_count) state types, got $(actual_count)."
end

"""    
    check_params(component, pas)

Validate that parameter collection contains all required parameters.

# Arguments
- `component::AbstractComponent`: Component to validate
- `pas::ComponentVector`: Parameter collection
"""
function check_params(component::AbstractComponent, pas::ComponentVector)
    param_names = get_param_names(component)
    component_name = get_name(component)
    available_params = keys(pas[:params])
    
    for param_name in param_names
        @assert param_name in available_params "Parameter '$(param_name)' in component '$(component_name)' is required but not found in parameters. Available parameters: $(available_params)"
    end
end

"""    
    check_initstates(component, initstates)

Validate that initial state collection contains all required states.

# Arguments
- `component::AbstractComponent`: Component to validate
- `initstates::ComponentVector`: Initial state collection
"""
function check_initstates(component::AbstractComponent, initstates::ComponentVector)
    state_names = get_state_names(component)
    component_name = get_name(component)
    available_states = keys(initstates)
    
    for state_name in state_names
        @assert state_name in available_states "Initial state '$(state_name)' in component '$(component_name)' is required but not found in initial states. Available states: $(available_states)"
    end
end

"""    
    check_nns(component, pas)

Validate that parameter collection contains all required neural networks.

# Arguments
- `component::AbstractComponent`: Component to validate
- `pas::ComponentVector`: Parameter collection
"""
function check_nns(component::AbstractComponent, pas::ComponentVector)
    nn_names = get_nn_names(component)
    component_name = get_name(component)
    
    if !isempty(nn_names)
        available_nns = keys(pas[:nns])
        
        for nn_name in nn_names
            @assert nn_name in available_nns "Neural network '$(nn_name)' in component '$(component_name)' is required but not found in neural networks. Available networks: $(available_nns)"
        end
    end
end