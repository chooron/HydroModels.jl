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

    available_vars = copy(input_names)

    for component in components
        comp_input_names = get_input_names(component)
        tmp_input_idx = map(comp_input_names) do nm
            findfirst(==(nm), available_vars)
        end

        if nothing in tmp_input_idx
            missing_vars = comp_input_names[tmp_input_idx .=== nothing]
            @warn "Input variables $missing_vars not found in available_vars"
        end

        push!(input_idx, filter(!isnothing, tmp_input_idx))

        comp_output_state = vcat(get_state_names(component), get_output_names(component))
        available_vars = vcat(available_vars, comp_output_state)
    end

    for name in var_names
        idx = findfirst(==(name), available_vars)
        if !isnothing(idx)
            push!(output_idx, idx)
        end
    end

    return input_idx, output_idx
end

"""
    _extract_component_states(initstates, comp_state_names, model_state_names, num_nodes)

Extract component-level initial states from model-level initstates for multi-node (3D) computation.

For ComponentVector (named), indexes by Symbol names directly.
For plain Vector (multi-node flat layout), extracts by positional index.
The flat vector layout is: [n1_s1, n2_s1, ..., nN_s1, n1_s2, ..., nN_s2, ...]
where each state variable occupies a contiguous block of `num_nodes` elements.
"""
function _extract_component_states(initstates::ComponentVector, comp_state_names, model_state_names, num_nodes)
    return initstates[comp_state_names]
end

function _extract_component_states(initstates::AbstractVector, comp_state_names, model_state_names, num_nodes)
    indices = Int[]
    for sname in comp_state_names
        global_idx = findfirst(==(sname), model_state_names)
        isnothing(global_idx) && error("State $sname not found in model states")
        append!(indices, ((global_idx - 1) * num_nodes + 1):(global_idx * num_nodes))
    end
    return initstates[indices]
end

"""
    _normalize_component_configs(config, n_components)

Normalize model configuration to per-component HydroConfig tuple.
"""
function _normalize_component_configs(config, n_components::Int)
    if config isa Tuple
        length(config) == n_components || error("Component configs length must equal components length")
        return map(normalize_config, config)
    elseif config isa NamedTuple && haskey(config, :components)
        comp_cfgs = config.components
        comp_cfgs isa Tuple || error("components field in config must be a Tuple")
        length(comp_cfgs) == n_components || error("Component configs length must equal components length")
        return map(normalize_config, comp_cfgs)
    else
        return ntuple(_ -> normalize_config(config), n_components)
    end
end

@inline function _slice_row(blocks::Vector, idx::Int, ::Val{2})
    row_idx = idx
    for block in blocks
        n_rows = size(block, 1)
        if row_idx <= n_rows
            return @view block[row_idx:row_idx, :]
        end
        row_idx -= n_rows
    end
    error("Row index $idx out of bounds for concatenated blocks")
end

@inline function _slice_row(blocks::Vector, idx::Int, ::Val{3})
    row_idx = idx
    for block in blocks
        n_rows = size(block, 1)
        if row_idx <= n_rows
            return @view block[row_idx:row_idx, :, :]
        end
        row_idx -= n_rows
    end
    error("Row index $idx out of bounds for concatenated blocks")
end

function _collect_rows(blocks::Vector, indices::AbstractVector{Int}, ::Val{2})
    first_block = blocks[1]
    T = eltype(first_block)
    isempty(indices) && return zeros(T, 0, size(first_block, 2))
    rows = map(i -> _slice_row(blocks, i, Val(2)), indices)
    return vcat(rows...)
end

function _collect_rows(blocks::Vector, indices::AbstractVector{Int}, ::Val{3})
    first_block = blocks[1]
    T = eltype(first_block)
    isempty(indices) && return zeros(T, 0, size(first_block, 2), size(first_block, 3))
    rows = map(i -> _slice_row(blocks, i, Val(3)), indices)
    return cat(rows...; dims=1)
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

It sequentially runs each component in the model, automatically handling the routing of data between components.
It supports both 2D (single-node) and 3D (multi-node) inputs.

Common kwargs include `initstates` and `config` (for component-specific settings).

# Arguments
- `input`: Input array with dimensions (variables, timesteps) or (variables, nodes, timesteps)
- `params`: Parameter vector (ComponentVector or AbstractVector)
- `config`: Configuration object or configuration tuple (for multi-component)
- `kwargs`: Additional keyword arguments

# Returns
- Output array containing all state and output variables
"""

function (model::HydroModel)(
    input::AbstractArray{T,2},
    params::AbstractVector,
    config::Union{ConfigType,Tuple}=default_config();
    kwargs...
) where {T}
    params = _as_componentvector(params)
    comp_configs = _normalize_component_configs(config, length(model.components))

    @assert size(input, 1) == length(get_input_names(model)) "Input variables length mismatch"

    initstates = get(kwargs, :initstates, model._defaultstates)
    outputs_blocks = Any[input]

    for (idx_, comp_, config_) in zip(model._varindices, model.components, comp_configs)
        tmp_input = _collect_rows(outputs_blocks, idx_, Val(2))

        comp_state_names = get_state_names(comp_)
        comp_initstates = !isempty(comp_state_names) ? initstates[comp_state_names] : nothing

        tmp_output = if !isnothing(comp_initstates)
            comp_(tmp_input, params, config_; initstates=comp_initstates)
        else
            comp_(tmp_input, params, config_)
        end

        push!(outputs_blocks, tmp_output)
    end

    return _collect_rows(outputs_blocks, model._outputindices, Val(2))
end

function (model::HydroModel)(
    input::AbstractArray{T,3},
    params::AbstractVector,
    config::Union{ConfigType,Tuple}=default_config();
    kwargs...
) where {T}
    params = _as_componentvector(params)
    comp_configs = _normalize_component_configs(config, length(model.components))

    @assert size(input, 1) == length(get_input_names(model)) "Input variables length mismatch"

    num_nodes = size(input, 2)
    initstates_raw = get(kwargs, :initstates, model._defaultstates)
    model_state_names = get_state_names(model)
    outputs_blocks = Any[input]

    for (idx_, comp_, config_) in zip(model._varindices, model.components, comp_configs)
        tmp_input = _collect_rows(outputs_blocks, idx_, Val(3))

        comp_state_names = get_state_names(comp_)
        comp_initstates = if !isempty(comp_state_names)
            _extract_component_states(initstates_raw, comp_state_names, model_state_names, num_nodes)
        else
            nothing
        end

        tmp_output = if !isnothing(comp_initstates)
            comp_(tmp_input, params, config_; initstates=comp_initstates)
        else
            comp_(tmp_input, params, config_)
        end

        push!(outputs_blocks, tmp_output)
    end

    return _collect_rows(outputs_blocks, model._outputindices, Val(3))
end

# Export interfaces
export HydroModel, @hydromodel
