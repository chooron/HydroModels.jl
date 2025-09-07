"""
    HydroModel{CS,NT,VI,OI,DS} <: AbstractModel

Represents a complete hydrological model by integrating a sequence of components. It manages the data flow between them automatically.

The model orchestrates how inputs and the outputs of preceding components are routed as inputs to subsequent components.

# Arguments
- `components::Tuple`: A tuple of computational components (e.g., `HydroBucket`, `UnitHydrograph`) that are executed in sequence.
- `name::Optional{Symbol}=nothing`: An optional identifier for the model. A unique name is generated if not provided.

# Fields
- `name::Symbol`: The identifier for the model.
- `components::Tuple`: The tuple of hydrological components that constitute the model.
- `infos::NamedTuple`: Aggregated metadata from all components, including the names of all `inputs`, `outputs`, `states`, `params`, and `nns`.
- `_varindices::Vector{Vector{Int}}`: Internal field. Stores the input variable indices for each component to manage data routing.
- `_outputindices::Vector{Int}`: Internal field. Stores the indices used to select and order the final model outputs from all intermediate variables.
- `_defaultstates::ComponentVector`: A `ComponentVector` holding the default initial values (zeros) for all state variables in the model.
"""
struct HydroModel{CS,NT,VI,OI,DS} <: AbstractModel
    "hydrological mode name"
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
        params = reduce(union, get_param_names.(components))
        nns = reduce(union, get_nn_names.(components))
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
"""
function _prepare_indices(components::CT, input_names::Vector{Symbol}, var_names::Vector{Symbol}) where {CT}
    input_idx, output_idx = Vector{Int}[], Vector{Int}()
    for component in components
        tmp_input_idx = map((nm) -> findfirst(varnm -> varnm == nm, input_names), get_input_names(component))
        if nothing in tmp_input_idx
            @warn "input variable $(get_input_names(component)) not found in input_names"
        else
            push!(input_idx, tmp_input_idx)
        end
        tmp_cpt_vcat_names = vcat(get_state_names(component), get_output_names(component))
        input_names = vcat(input_names, tmp_cpt_vcat_names)
    end
    for name in var_names
        push!(output_idx, findfirst(varnm -> varnm == name, input_names))
    end
    return input_idx, output_idx
end

"""
    @hydromodel name begin ... end

A macro to conveniently create a `HydroModel` from a sequence of components.

The components listed within the `begin...end` block are automatically collected into a tuple and passed to the `HydroModel` constructor.

# Arguments
- `name`: A `Symbol` that provides a name for the model.
- The block should contain a list of pre-defined component instances (e.g., `HydroBucket` objects).
"""
macro hydromodel(args...)
    name = length(args) == 1 ? nothing : args[1]
    expr = length(args) == 1 ? args[1] : args[2]

    @assert Meta.isexpr(expr, :block) "Expected a begin...end block after model name"
    components = filter(x -> !(x isa LineNumberNode), expr.args)
    return esc(:(HydroModel(; name=$(name), components=tuple($(components...)))))
end

"""
    (model::HydroModel)(input::AbstractArray, params::ComponentVector; kwargs...)

Executes the hydrological model simulation.

This function sequentially runs each component in the model, handling the flow of data from inputs and previous components to the next.

# Arguments
- `input::AbstractArray`: The input data for the model. Must be a 2D array `(variables, timesteps)` for a single-node simulation or a 3D array `(variables, nodes, timesteps)` for a multi-node simulation.
- `params::ComponentVector`: A `ComponentVector` containing all parameters required by the model's components.

# Keyword Arguments
- `initstates::AbstractVector=model._defaultstates`: A vector or `ComponentVector` specifying the initial states for any stateful components. Defaults to a zero vector for all states.
- `config::NamedTuple=NamedTuple()`: Configuration options (e.g., ODE solver settings) to be passed down to the components.

# Returns
- `AbstractArray`: The final computed model output, with dimensions corresponding to the input `(output_variables, timesteps)` or `(output_variables, nodes, timesteps)`.
"""
function (model::HydroModel)(
    input::AbstractArray{T,D},
    params::ComponentVector;
    initstates::AbstractVector=model._defaultstates,
    config::NamedTuple=NamedTuple(),
) where {T,D}
    comp_configs = config isa NamedTuple ? fill(config, length(model.components)) : config
    @assert D in (2, 3) "input array dimension must be 2 or 3"
    @assert length(comp_configs) == length(model.components) "component configs length must be equal to components length"
    @assert size(input, 1) == length(get_input_names(model)) "input variables length must be equal to input variables length"
    outputs = input
    for (idx_, comp_, config_) in zip(model._varindices, model.components, comp_configs)
        tmp_input = copy(view(outputs, idx_, ntuple(_ -> Colon(), D - 1)...))
        tmp_output = comp_(
            tmp_input, params;
            initstates=initstates[get_state_names(comp_)],
            config=config_
        )
        outputs = vcat(outputs, tmp_output)
    end
    return outputs[model._outputindices, ntuple(_ -> Colon(), D - 1)...]
end