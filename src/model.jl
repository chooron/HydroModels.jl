"""
    HydroModel{CS, NT, VI, OI, DS} <: AbstractModel

Represents a complete hydrological model composed of a sequence of components (e.g., `HydroBucket`, `HydroRoute`).

It automatically manages the data flow between components, routing the outputs of one component as inputs to the next.

$(FIELDS)
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
    @hydromodel [name] begin ... end

A macro to conveniently create a `HydroModel` from a sequence of components.

# Usage
The macro takes an optional name and a `begin...end` block containing an ordered list of component instances.

```julia
# Assuming bucket1 and route1 are pre-defined components
@hydromodel :my_full_model begin
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
    (model::HydroModel)(input, params; kwargs...)

Executes the full hydrological model simulation. This is the functor implementation for `HydroModel`.

It sequentially runs each component in the model, automatically handling the routing of data between them. It supports both 2D (single-node) and 3D (multi-node) inputs.

Common `kwargs` include `initstates` and `config` (for component-specific settings).
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