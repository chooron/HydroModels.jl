"""
BMI implementation for lumped (0-dimensional) hydrological models.

This module provides a Basic Model Interface (BMI) wrapper for HydroModels components,
enabling interoperability with other modeling frameworks.
"""

using BasicModelInterface
const BMI = BasicModelInterface
using ComponentArrays
using Dates

"""
    LumpedModelBMI{C}

BMI wrapper for lumped (single-point) hydrological models.

# Fields
- `component::C`: The HydroModels component (HydroBucket, HydroModel, etc.)
- `start_time::Date`: Simulation start time
- `end_time::Date`: Simulation end time
- `current_time::Date`: Current simulation time
- `dt::Number`: Time step size
- `time_unit::String`: Time unit (default: "days")
- `time_index::Int`: Current time step index
- `timeidx::AbstractVector{<:Date}`: Complete time vector
- `forcing::AbstractArray`: Forcing data (inputs × timesteps)
- `pas::ComponentVector`: Model parameters
- `initstates::ComponentVector`: Initial states
- `config::NamedTuple`: Model configuration
- `current_state::ComponentVector`: Current state values
- `outputs::Union{Nothing,AbstractArray}`: Stored model outputs
- `lat::Union{Nothing,Float64}`: Latitude (optional)
- `lon::Union{Nothing,Float64}`: Longitude (optional)
"""
@kwdef mutable struct LumpedModelBMI{C}
    component::C
    start_time::Date
    end_time::Date
    current_time::Date
    dt::Number
    time_unit::String = "days"
    time_index::Int = 1
    timeidx::AbstractVector{<:Date}
    forcing::AbstractArray
    pas::ComponentVector
    initstates::ComponentVector
    config::NamedTuple
    current_state::ComponentVector
    outputs::Union{Nothing,AbstractArray} = nothing
    lat::Union{Nothing,Float64} = nothing
    lon::Union{Nothing,Float64} = nothing
end

"""
    validate_bmi_inputs(component, forcing, pas, initstates, timeidx)

Validate inputs for BMI initialization.

# Arguments
- `component`: HydroModels component
- `forcing`: Forcing data array
- `pas`: Parameter ComponentVector
- `initstates`: Initial states ComponentVector
- `timeidx`: Time index vector

# Throws
- Error if validation fails
"""
function validate_bmi_inputs(
    component::AbstractComponent,
    forcing::AbstractArray,
    pas::ComponentVector,
    initstates::ComponentVector,
    timeidx::AbstractVector{<:Date}
)
    # Check forcing dimensions
    input_names = get_input_names(component)
    if size(forcing, 1) != length(input_names)
        error("Forcing array first dimension ($(size(forcing, 1))) must match number of inputs ($(length(input_names)))")
    end
    if size(forcing, 2) != length(timeidx)
        error("Forcing array second dimension ($(size(forcing, 2))) must match timeidx length ($(length(timeidx)))")
    end

    # Check parameters
    param_names = get_param_names(component)
    for pname in param_names
        if !haskey(pas, pname)
            error("Missing required parameter: $pname")
        end
    end

    # Check initial states
    state_names = get_state_names(component)
    for sname in state_names
        if !haskey(initstates, sname)
            error("Missing required initial state: $sname")
        end
    end

    # Check time vector
    if length(timeidx) < 2
        error("Time vector must have at least 2 elements")
    end
    if !issorted(timeidx)
        error("Time vector must be sorted")
    end

    return true
end

"""
    BMI.initialize(component, forcing, pas, initstates, timeidx, config; lat, lon)

Initialize a BMI model from a HydroModels component.

# Arguments
- `component`: HydroModels component (HydroBucket, HydroModel, etc.)
- `forcing`: Forcing data array (n_inputs × n_timesteps)
- `pas`: Parameter ComponentVector
- `initstates`: Initial states ComponentVector
- `timeidx`: Time index vector (Date array)
- `config`: Configuration NamedTuple

# Keyword Arguments
- `lat::Union{Nothing,Float64}`: Latitude (optional)
- `lon::Union{Nothing,Float64}`: Longitude (optional)

# Returns
- `LumpedModelBMI`: Initialized BMI model
"""
function BMI.initialize(
    component::C,
    forcing::AbstractArray,
    pas::ComponentVector,
    initstates::ComponentVector,
    timeidx::AbstractVector{<:Date},
    config::NamedTuple;
    lat::Union{Nothing,Float64}=nothing,
    lon::Union{Nothing,Float64}=nothing
) where {C<:AbstractComponent}
    # Validate inputs
    validate_bmi_inputs(component, forcing, pas, initstates, timeidx)

    # Calculate time step
    dt = length(timeidx) > 1 ? (timeidx[2] - timeidx[1]).value : 1

    return LumpedModelBMI{C}(
        component=component,
        start_time=timeidx[1],
        end_time=timeidx[end],
        current_time=timeidx[1],
        dt=dt,
        time_unit="days",
        time_index=1,
        timeidx=timeidx,
        forcing=forcing,
        pas=pas,
        initstates=initstates,
        config=config,
        current_state=deepcopy(initstates),
        outputs=nothing,
        lat=lat,
        lon=lon
    )
end

"""
    BMI.update(model::LumpedModelBMI)

Advance the model by one time step.

Updates the model state and increments the time index.
"""
function BMI.update(model::LumpedModelBMI)
    # Check if we've reached the end
    if model.time_index > length(model.timeidx)
        error("Cannot update: already at end time")
    end

    # Get forcing for current time step
    current_forcing = model.forcing[:, model.time_index:model.time_index]

    # Run model for one time step
    output = model.component(
        current_forcing,
        model.pas,
        model.config,
        initstates=model.current_state
    )

    # Extract states from output
    state_names = get_state_names(model.component)
    num_states = length(state_names)

    if num_states > 0
        # Update current state (states are first rows of output)
        for (i, sname) in enumerate(state_names)
            model.current_state[sname] = output[i, end]
        end
    end

    # Store outputs
    model.outputs = output

    # Increment time
    model.time_index += 1
    if model.time_index <= length(model.timeidx)
        model.current_time = model.timeidx[model.time_index]
    end

    return nothing
end

"""
    BMI.update_until(model::LumpedModelBMI, time::Date)

Advance the model to a specific time.

# Arguments
- `model`: BMI model
- `time`: Target time

# Returns
- Model outputs up to the target time
"""
function BMI.update_until(model::LumpedModelBMI, time::Date)
    # Validate time
    if time < model.current_time
        error("Target time must be >= current time")
    end
    if time > model.end_time
        error("Target time must be <= end time")
    end

    # Find target index
    target_index = findfirst(t -> t >= time, model.timeidx)
    if isnothing(target_index)
        error("Target time $time not found in time vector")
    end

    # Get forcing for time range
    forcing_slice = model.forcing[:, model.time_index:target_index]

    # Run model
    output = model.component(
        forcing_slice,
        model.pas,
        model.config,
        initstates=model.current_state
    )

    # Update state
    state_names = get_state_names(model.component)
    num_states = length(state_names)

    if num_states > 0
        for (i, sname) in enumerate(state_names)
            model.current_state[sname] = output[i, end]
        end
    end

    # Store outputs
    model.outputs = output

    # Update time tracking
    model.time_index = target_index
    model.current_time = model.timeidx[target_index]

    return output
end

"""
    BMI.finalize(model::LumpedModelBMI)

Complete the simulation and clean up memory.

# Returns
- Final model outputs
"""
function BMI.finalize(model::LumpedModelBMI)
    # Run to end if not already there
    if model.current_time < model.end_time
        output = BMI.update_until(model, model.end_time)
    else
        output = model.outputs
    end

    # Clear large data structures to free memory
    model.outputs = nothing
    model.forcing = nothing

    return output
end

# ============================================================================
# Information Functions
# ============================================================================

get_component_name(model::LumpedModelBMI) = get_name(model.component)

get_input_item_count(model::LumpedModelBMI) = length(get_input_names(model.component))

get_output_item_count(model::LumpedModelBMI) = length(get_output_names(model.component)) + length(get_state_names(model.component))

get_input_var_names(model::LumpedModelBMI) = string.(get_input_names(model.component))

get_output_var_names(model::LumpedModelBMI) = string.(vcat(get_state_names(model.component), get_output_names(model.component)))

# ============================================================================
# Variable Information Functions
# ============================================================================

get_var_grid(::LumpedModelBMI, ::String) = 0

get_var_type(model::LumpedModelBMI, ::String) = repr(eltype(model.forcing))

function get_var_units(model::LumpedModelBMI, name::String)
    name_sym = Symbol(name)
    infos = model.component.infos

    if haskey(infos.outputs, name_sym)
        return get(infos.outputs[name_sym], :unit, "")
    elseif haskey(infos.inputs, name_sym)
        return get(infos.inputs[name_sym], :unit, "")
    elseif haskey(infos.states, name_sym)
        return get(infos.states[name_sym], :unit, "")
    else
        error("Variable $name not found in model")
    end
end

get_var_itemsize(model::LumpedModelBMI, ::String) = sizeof(eltype(model.forcing))

get_var_nbytes(model::LumpedModelBMI, name::String) = get_var_itemsize(model, name) * length(model.timeidx)

get_var_location(model::LumpedModelBMI, name::String) = "node"

# ============================================================================
# Time Functions
# ============================================================================

get_current_time(model::LumpedModelBMI) = model.current_time

get_start_time(model::LumpedModelBMI) = model.start_time

get_end_time(model::LumpedModelBMI) = model.end_time

get_time_units(model::LumpedModelBMI) = model.time_unit

get_time_step(model::LumpedModelBMI) = model.dt

# ============================================================================
# Getter Functions
# ============================================================================

function get_value(model::LumpedModelBMI, name::String, dest::AbstractArray)
    dest .= copy(get_value_ptr(model, name))
    return dest
end

function get_value_ptr(model::LumpedModelBMI, name::String)
    # Check if we have outputs
    if isnothing(model.outputs)
        error("No outputs available. Run update() or update_until() first.")
    end

    # Get variable names
    state_names = get_state_names(model.component)
    output_names = get_output_names(model.component)
    input_names = get_input_names(model.component)

    name_sym = Symbol(name)

    # Check states
    if name_sym in state_names
        var_index = findfirst(==(name_sym), state_names)
        return model.outputs[var_index, :]
    end

    # Check outputs (come after states)
    if name_sym in output_names
        var_index = findfirst(==(name_sym), output_names)
        return model.outputs[length(state_names) + var_index, :]
    end

    # Check inputs (from forcing)
    if name_sym in input_names
        var_index = findfirst(==(name_sym), input_names)
        return model.forcing[var_index, :]
    end

    error("Variable $name not found in model")
end

function get_value_at_indices(model::LumpedModelBMI, name::String, dest::AbstractArray, indices::AbstractArray{Int})
    values = get_value_ptr(model, name)
    dest .= values[indices]
    return dest
end

# ============================================================================
# Setter Functions
# ============================================================================

function set_value(model::LumpedModelBMI, name::String, src::AbstractArray)
    name_sym = Symbol(name)

    # Get variable names
    state_names = get_state_names(model.component)
    input_names = get_input_names(model.component)

    # Can only set states or inputs
    if name_sym in state_names
        # Update current state
        model.current_state[name_sym] = src[1]  # Set to first value

        # Update outputs if they exist
        if !isnothing(model.outputs)
            var_index = findfirst(==(name_sym), state_names)
            model.outputs[var_index, :] .= src
        end
    elseif name_sym in input_names
        # Update forcing
        var_index = findfirst(==(name_sym), input_names)
        model.forcing[var_index, :] .= src
    else
        error("Cannot set variable $name (not a state or input)")
    end

    return nothing
end

function set_value_at_indices(
    model::LumpedModelBMI,
    name::String,
    indices::AbstractArray{Int},
    src::AbstractArray
)
    name_sym = Symbol(name)

    # Get variable names
    state_names = get_state_names(model.component)
    input_names = get_input_names(model.component)

    if name_sym in state_names
        if !isnothing(model.outputs)
            var_index = findfirst(==(name_sym), state_names)
            model.outputs[var_index, indices] .= src
        end
        # Update current state if last index is included
        if maximum(indices) == size(model.outputs, 2)
            model.current_state[name_sym] = src[end]
        end
    elseif name_sym in input_names
        var_index = findfirst(==(name_sym), input_names)
        model.forcing[var_index, indices] .= src
    else
        error("Cannot set variable $name at indices")
    end

    return nothing
end

# ============================================================================
# Grid Functions (for lumped models, all return fixed values)
# ============================================================================

get_grid_rank(::LumpedModelBMI, ::Int) = 2

get_grid_size(::LumpedModelBMI, ::Int) = 1

get_grid_type(::LumpedModelBMI, ::Int) = "uniform_rectilinear"

get_grid_shape(::LumpedModelBMI, ::Int, ::AbstractVector{Int}) = (1, 1)

get_grid_spacing(::LumpedModelBMI, ::Int, ::AbstractVector{Float64}) = (0.0, 0.0)

function get_grid_origin(model::LumpedModelBMI, ::Int, dest::AbstractVector{Float64})
    # For lumped models, return lat/lon if provided, otherwise (0, 0)
    if !isnothing(model.lat) && !isnothing(model.lon)
        dest[1] = model.lat
        dest[2] = model.lon
    else
        dest[1] = 0.0
        dest[2] = 0.0
    end
    return dest
end

get_grid_x(::LumpedModelBMI, ::Int, ::AbstractVector{Float64}) = 1

get_grid_y(::LumpedModelBMI, ::Int, ::AbstractVector{Float64}) = 1

get_grid_z(::LumpedModelBMI, ::Int, ::AbstractVector{Float64}) = 1

get_grid_node_count(::LumpedModelBMI, ::Int) = 1

get_grid_edge_count(::LumpedModelBMI, ::Int) = 1

get_grid_face_count(::LumpedModelBMI, ::Int) = 1

function get_grid_edge_nodes(model::LumpedModelBMI, grid::Int, edge_nodes::AbstractVector{Int})
    # For lumped model, single edge connecting the node to itself
    edge_nodes[1] = 1
    edge_nodes[2] = 1
    return edge_nodes
end

get_grid_face_edges(::LumpedModelBMI, ::Int, ::AbstractVector{Int}) = 1

get_grid_face_nodes(::LumpedModelBMI, ::Int, ::AbstractVector{Int}) = 1

get_grid_nodes_per_face(::LumpedModelBMI, ::Int, ::AbstractVector{Int}) = 1

# ============================================================================
# Helper Functions
# ============================================================================

"""
    save_state(model::LumpedModelBMI)

Save current model state for later restoration.

# Returns
- NamedTuple with current state, time index, current time, and outputs
"""
function save_state(model::LumpedModelBMI)
    return (
        current_state = deepcopy(model.current_state),
        time_index = model.time_index,
        current_time = model.current_time,
        outputs = isnothing(model.outputs) ? nothing : copy(model.outputs)
    )
end

"""
    restore_state!(model::LumpedModelBMI, state)

Restore model to a previously saved state.

# Arguments
- `model`: BMI model
- `state`: State tuple from save_state()
"""
function restore_state!(model::LumpedModelBMI, state)
    model.current_state = deepcopy(state.current_state)
    model.time_index = state.time_index
    model.current_time = state.current_time
    model.outputs = isnothing(state.outputs) ? nothing : copy(state.outputs)
    return nothing
end

"""
    step_forward!(model::LumpedModelBMI, n_steps::Int=1)

Step the model forward by n time steps.

# Arguments
- `model`: BMI model
- `n_steps`: Number of steps to advance (default: 1)
"""
function step_forward!(model::LumpedModelBMI, n_steps::Int=1)
    for _ in 1:n_steps
        if model.time_index > length(model.timeidx)
            @warn "Already at end time"
            break
        end
        BMI.update(model)
    end
    return nothing
end

"""
    get_time_step_count(model::LumpedModelBMI)

Get total number of time steps.
"""
get_time_step_count(model::LumpedModelBMI) = length(model.timeidx)

"""
    get_remaining_steps(model::LumpedModelBMI)

Get number of remaining time steps.
"""
get_remaining_steps(model::LumpedModelBMI) = length(model.timeidx) - model.time_index + 1
