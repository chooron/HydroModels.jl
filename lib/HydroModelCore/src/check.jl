# ================================================================================================
# 主检查函数
# ================================================================================================

"""
    check(component, input, pas, initstates, timeidx)

$(SIGNATURES)

Validate component inputs, parameters, initial states, and neural networks.

# Arguments
- `component::AbstractComponent`: Component to validate
- `input::AbstractArray{<:Number}`: Input data (2D or 3D array)
- `pas::ComponentVector`: Parameter collection
- `initstates::ComponentVector`: Initial states collection
- `timeidx::AbstractVector`: Time indices

# Details
This function performs comprehensive validation including:
- Input dimensions and variable count
- Parameter existence and types
- Initial state existence and dimensions
- Neural network existence

Dispatches to specialized checking functions based on input dimensionality.
"""
function check(component::AbstractComponent, input::AbstractArray{<:Number}, pas::ComponentVector, initstates::ComponentVector, timeidx::AbstractVector)
    dims = ndims(input)
    check_input(component, input, timeidx, Val(dims))
    check_params(component, pas)
    check_initstates(component, initstates, Val(dims))
    check_nns(component, pas)
    return nothing
end


# ================================================================================================
# 输入检查（使用 Val 分派优化）
# ================================================================================================

"""
    check_input(component, input, timeidx, ::Val{dims})

$(SIGNATURES)

Validate input data dimensions against component requirements.
Uses Val dispatch for compile-time optimization based on input dimensionality.

# Arguments
- `component::AbstractComponent`: Component to validate
- `input::AbstractArray`: Input data
- `timeidx::AbstractVector`: Time indices
- `::Val{dims}`: Dimensionality indicator (2 or 3)
"""
function check_input(component::AbstractComponent, input::AbstractArray, timeidx::AbstractVector, ::Val{dims}) where {dims}
    _check_input_vars(component, input, dims)
    _check_input_time(component, input, timeidx, dims)
    return nothing
end

"""
$(SIGNATURES)

Check if the number of input variables matches component requirements.
"""
@inline function _check_input_vars(component::AbstractComponent, input::AbstractArray, dims::Int)
    component_name = get_name(component)
    input_names = get_input_names(component)
    expected_vars = length(input_names)
    actual_vars = size(input, 1)

    if actual_vars != expected_vars
        throw(DimensionMismatch(
            "Input variables in component '$component_name' do not match required dimensions.\n" *
            "  Expected: $expected_vars variables $(input_names)\n" *
            "  Got:      $actual_vars variables\n" *
            "  Hint: Check that input array has shape ($(expected_vars), timesteps$(dims == 3 ? ", grids" : ""))"
        ))
    end
    return nothing
end

"""
$(SIGNATURES)

Check if the number of timesteps matches timeidx length.
"""
@inline function _check_input_time(component::AbstractComponent, input::AbstractArray, timeidx::AbstractVector, dims::Int)
    component_name = get_name(component)
    time_dim = dims == 2 ? 2 : 3  # For 2D: (vars, time), For 3D: (vars, grids, time)
    expected_timesteps = length(timeidx)
    actual_timesteps = size(input, time_dim)

    if actual_timesteps != expected_timesteps
        throw(DimensionMismatch(
            "Time steps in component '$component_name' do not match required length.\n" *
            "  Expected: $expected_timesteps steps\n" *
            "  Got:      $actual_timesteps steps\n" *
            "  Hint: timeidx length should match input dimension $time_dim"
        ))
    end
    return nothing
end


# ================================================================================================
# 参数检查
# ================================================================================================

"""
    check_params(component, pas)

$(SIGNATURES)

Validate that parameter collection contains all required parameters.

# Arguments
- `component::AbstractComponent`: Component to validate
- `pas::ComponentVector`: Parameter collection

# Throws
- `KeyError`: If any required parameter is missing
"""
function check_params(component::AbstractComponent, pas::ComponentVector)
    param_names = get_param_names(component)
    isempty(param_names) && return nothing

    component_name = get_name(component)
    available_params = keys(pas[:params])

    missing_params = Symbol[]
    for param_name in param_names
        if !(param_name in available_params)
            push!(missing_params, param_name)
        end
    end

    if !isempty(missing_params)
        throw(KeyError(
            "Missing parameters in component '$component_name':\n" *
            "  Required but missing: $(missing_params)\n" *
            "  Available parameters: $(collect(available_params))\n" *
            "  Hint: Initialize missing parameters before running the component"
        ))
    end
    return nothing
end

"""
$(SIGNATURES)

Check parameter types and values for validity.

# Arguments
- `component::AbstractComponent`: Component to validate
- `pas::ComponentVector`: Parameter collection
- `allow_negative::Bool=false`: Whether negative values are allowed
"""
function check_param_values(component::AbstractComponent, pas::ComponentVector; allow_negative::Bool=false)
    param_names = get_param_names(component)
    isempty(param_names) && return nothing

    component_name = get_name(component)

    for param_name in param_names
        param_value = pas[:params][param_name]

        # Check for NaN
        if any(isnan, param_value)
            @warn "Parameter '$param_name' in component '$component_name' contains NaN values"
        end

        # Check for Inf
        if any(isinf, param_value)
            @warn "Parameter '$param_name' in component '$component_name' contains Inf values"
        end

        # Check for negative values if not allowed
        if !allow_negative && any(<(0), param_value)
            @warn "Parameter '$param_name' in component '$component_name' contains negative values"
        end
    end
    return nothing
end


# ================================================================================================
# 初始状态检查（使用 Val 分派）
# ================================================================================================

"""
    check_initstates(component, initstates, ::Val{dims})

$(SIGNATURES)

Validate that initial state collection contains all required states.
Uses Val dispatch to handle different dimensionalities.

# Arguments
- `component::AbstractComponent`: Component to validate
- `initstates::ComponentVector`: Initial state collection
- `::Val{dims}`: Dimensionality indicator
"""
function check_initstates(component::AbstractComponent, initstates::ComponentVector, ::Val{dims}) where {dims}
    state_names = get_state_names(component)
    isempty(state_names) && return nothing

    component_name = get_name(component)
    available_states = keys(initstates)

    missing_states = Symbol[]
    for state_name in state_names
        if !(state_name in available_states)
            push!(missing_states, state_name)
        end
    end

    if !isempty(missing_states)
        throw(KeyError(
            "Missing initial states in component '$component_name':\n" *
            "  Required but missing: $(missing_states)\n" *
            "  Available states: $(collect(available_states))\n" *
            "  Hint: Initialize missing states before running the component"
        ))
    end

    # Check dimensions if dims is specified
    _check_state_dimensions(component, initstates, Val(dims))
    return nothing
end

"""
$(SIGNATURES)

Check that state dimensions are consistent with input dimensionality.
"""
@inline function _check_state_dimensions(component::AbstractComponent, initstates::ComponentVector, ::Val{2})
    # For 2D input, states should be scalars or 1D arrays
    state_names = get_state_names(component)
    component_name = get_name(component)

    for state_name in state_names
        state_value = initstates[state_name]
        if ndims(state_value) > 1
            @warn "Initial state '$state_name' in component '$component_name' has dimension $(ndims(state_value)), expected scalar or 1D for 2D input"
        end
    end
    return nothing
end

@inline function _check_state_dimensions(component::AbstractComponent, initstates::ComponentVector, ::Val{3})
    # For 3D input, states should be 1D or 2D arrays
    state_names = get_state_names(component)
    component_name = get_name(component)

    for state_name in state_names
        state_value = initstates[state_name]
        if ndims(state_value) > 2
            @warn "Initial state '$state_name' in component '$component_name' has dimension $(ndims(state_value)), expected 1D or 2D for 3D input"
        end
    end
    return nothing
end

"""
$(SIGNATURES)

Check initial state values for validity.
"""
function check_state_values(component::AbstractComponent, initstates::ComponentVector)
    state_names = get_state_names(component)
    isempty(state_names) && return nothing

    component_name = get_name(component)

    for state_name in state_names
        state_value = initstates[state_name]

        # Check for NaN
        if any(isnan, state_value)
            @warn "Initial state '$state_name' in component '$component_name' contains NaN values"
        end

        # Check for Inf
        if any(isinf, state_value)
            @warn "Initial state '$state_name' in component '$component_name' contains Inf values"
        end
    end
    return nothing
end


# ================================================================================================
# 神经网络检查
# ================================================================================================

"""
    check_nns(component, pas)

$(SIGNATURES)

Validate that parameter collection contains all required neural networks.

# Arguments
- `component::AbstractComponent`: Component to validate
- `pas::ComponentVector`: Parameter collection

# Throws
- `KeyError`: If any required neural network is missing
"""
function check_nns(component::AbstractComponent, pas::ComponentVector)
    nn_names = get_nn_names(component)
    isempty(nn_names) && return nothing

    component_name = get_name(component)
    available_nns = keys(pas[:nns])

    missing_nns = Symbol[]
    for nn_name in nn_names
        if !(nn_name in available_nns)
            push!(missing_nns, nn_name)
        end
    end

    if !isempty(missing_nns)
        throw(KeyError(
            "Missing neural networks in component '$component_name':\n" *
            "  Required but missing: $(missing_nns)\n" *
            "  Available networks: $(collect(available_nns))\n" *
            "  Hint: Initialize missing neural networks before running the component"
        ))
    end
    return nothing
end


# ================================================================================================
# 批量检查辅助函数
# ================================================================================================

"""
$(SIGNATURES)

Validate multiple components at once.

# Arguments
- `components::AbstractVector{<:AbstractComponent}`: Components to validate
- `input::AbstractArray`: Input data
- `pas::ComponentVector`: Parameter collection
- `initstates::ComponentVector`: Initial states collection
- `timeidx::AbstractVector`: Time indices

# Returns
- `Vector{Bool}`: Vector indicating which components passed validation
"""
function check_multiple(components::AbstractVector{<:AbstractComponent}, input::AbstractArray, pas::ComponentVector, initstates::ComponentVector, timeidx::AbstractVector)
    results = Bool[]
    for component in components
        try
            check(component, input, pas, initstates, timeidx)
            push!(results, true)
        catch e
            @warn "Component '$(get_name(component))' failed validation: $e"
            push!(results, false)
        end
    end
    return results
end

"""
$(SIGNATURES)

Check if all components in a collection are valid.
Returns true only if all components pass validation.
"""
function check_all(components::AbstractVector{<:AbstractComponent}, input::AbstractArray, pas::ComponentVector, initstates::ComponentVector, timeidx::AbstractVector)
    return all(check_multiple(components, input, pas, initstates, timeidx))
end


# ================================================================================================
# 依赖性检查
# ================================================================================================

"""
$(SIGNATURES)

Check if component dependencies are satisfied.
Verifies that all required inputs are available from previous components' outputs.

# Arguments
- `component::AbstractComponent`: Component to check
- `available_vars::AbstractVector{Symbol}`: Variables available from previous components

# Returns
- `Tuple{Bool, Vector{Symbol}}`: (is_satisfied, missing_vars)
"""
function check_dependencies(component::AbstractComponent, available_vars::AbstractVector{Symbol})
    required_inputs = get_input_names(component)
    missing_vars = setdiff(required_inputs, available_vars)
    return isempty(missing_vars), missing_vars
end

"""
$(SIGNATURES)

Check dependencies for a sequence of components.
Returns true if the component chain is valid (each component's inputs are satisfied by previous outputs).
"""
function check_component_chain(components::AbstractVector{<:AbstractComponent}, initial_vars::AbstractVector{Symbol}=Symbol[])
    available_vars = copy(initial_vars)

    for (idx, component) in enumerate(components)
        is_satisfied, missing_vars = check_dependencies(component, available_vars)

        if !is_satisfied
            component_name = get_name(component)
            throw(ErrorException(
                "Component chain validation failed at component #$idx ('$component_name'):\n" *
                "  Missing inputs: $(missing_vars)\n" *
                "  Available variables: $(available_vars)\n" *
                "  Hint: Reorder components or add missing input sources"
            ))
        end

        # Add this component's outputs to available vars
        union!(available_vars, get_output_names(component))
    end

    return true
end
