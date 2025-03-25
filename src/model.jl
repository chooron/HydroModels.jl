"""
	HydroModel <: AbstractModel

Represents a hydrological model composed of multiple components.

# Fields
- `infos::NamedTuple`: Contains metadata about the model, including name, input variables, all variables, output variables, state variables, and neural network variables.
- `components::Vector{<:AbstractComponent}`: A vector of hydrological computation elements (components) that make up the model.
- `varindices::Vector`: A vector of indices for each component's input, used to map overall model inputs to component-specific inputs.

# Constructor
	HydroModel(name; components::Vector{<:AbstractComponent})

Constructs a HydroModel with the given name and components.

# Description
HydroModel is a structure that encapsulates a complete hydrological model. It manages multiple hydrological components, 
handles the flow of data between these components, and provides methods for running simulations.

The model automatically determines the connections between components based on their input and output variables. 
It also keeps track of all variables in the system, including inputs, outputs, states, and any neural network parameters.

When called as a function, the HydroModel applies its components in sequence, passing the outputs of earlier components 
as inputs to later ones, effectively simulating the hydrological system over time.

Each component's kwargs may be different, include solver, interp
"""
struct HydroModel <: AbstractModel
    "name of hydrological model"
    name::Symbol
    "hydrological computation elements"
    components::Vector{<:AbstractComponent}
    "input variables index for each components"
    varindices::AbstractVector{<:AbstractVector{<:Integer}}
    "output variables index for sort output variables"
    outputindices::AbstractVector{<:Integer}
    "meta data of hydrological model"
    meta::ComponentVector

    function HydroModel(;
        components::Vector{C},
        name::Union{Symbol,Nothing}=nothing,
        sort_components::Bool=false,
    ) where {C<:AbstractComponent}
        components = sort_components ? sort_components(components) : components
        model_inputs, model_outputs, model_states = get_all_vars(components)
        model_params = reduce(union, get_param_vars.(components))
        model_nns_ntp = reduce(merge, map(flux -> NamedTuple(get_nn_vars(flux)), components))
        input_idx, output_idx = _prepare_indices(components, model_inputs, vcat(model_inputs, model_states, model_outputs))
        model_meta = ComponentVector(inputs=model_inputs, outputs=model_outputs, states=model_states, params=model_params, nns=model_nns_ntp)
        model_name = isnothing(name) ? Symbol("##model#", hash(model_meta)) : name
        new(model_name, components, input_idx, output_idx, model_meta)
    end
end

function _prepare_indices(components::Vector{<:AbstractComponent}, inputs::Vector{<:Num}, outputs::Vector{<:Num})
    input_idx, output_idx = Vector{Integer}[], Vector{Integer}()
    var_names = Symbolics.tosymbol.(inputs)
    vcat_names = Symbolics.tosymbol.(outputs)
    for component in components
        #* extract input index
        tmp_input_idx = map((nm) -> findfirst(varnm -> varnm == nm, var_names), get_input_names(component))
        push!(input_idx, tmp_input_idx)
        #* extract output index
        tmp_cpt_vcat_names = vcat(get_state_names(component), get_output_names(component))
        var_names = vcat(var_names, tmp_cpt_vcat_names)
        tmp_output_idx = map((nm) -> findfirst(varnm -> varnm == nm, vcat_names), tmp_cpt_vcat_names)
        output_idx = vcat(output_idx, tmp_output_idx)
    end
    return input_idx, output_idx
end

function (model::HydroModel)(
    input::AbstractArray{T,2},
    pas::ComponentVector;
    config::Union{NamedTuple,Vector{<:NamedTuple}}=NamedTuple(),
    kwargs...,
) where {T<:Number}
    comp_configs = config isa NamedTuple ? fill(config, length(model.components)) : config
    @assert length(comp_configs) == length(model.components) "component configs length must be equal to components length"
    outputs = input
    params = haskey(pas, :params) ? view(pas, :params) : ComponentVector()
    initstates = haskey(pas, :initstates) ? view(pas, :initstates) : ComponentVector()
    nns = haskey(pas, :nns) ? view(pas, :nns) : ComponentVector()
    for (idx_, (comp_, config_)) in enumerate(zip(model.components, comp_configs))
        extract_params = view(params, get_param_names(comp_))
        extract_initstates = view(initstates, get_state_names(comp_))
        extract_nns = view(nns, get_nn_names(comp_))
        tmp_pas = ComponentVector(params=extract_params, initstates=extract_initstates, nns=extract_nns)
        tmp_outputs = comp_(view(outputs, model.varindices[idx_], :), tmp_pas; config=config_)
        outputs = cat(outputs, tmp_outputs, dims=1)
    end
    return view(outputs, model.outputindices, :)
end

function (model::HydroModel)(
    input::AbstractArray{T,3},
    pas::ComponentVector;
    config::Union{<:NamedTuple,Vector{<:NamedTuple}}=NamedTuple(),
    kwargs...,
) where {T<:Number}
    comp_configs = config isa NamedTuple ? fill(config, length(model.components)) : config
    @assert length(comp_configs) == length(model.components) "component configs length must be equal to components length"
    outputs = input
    params = haskey(pas, :params) ? view(pas, :params) : ComponentVector()
    initstates = haskey(pas, :initstates) ? view(pas, :initstates) : ComponentVector()
    nns = haskey(pas, :nns) ? view(pas, :nns) : ComponentVector()
    for (idx_, (comp_, config_)) in enumerate(zip(model.components, comp_configs))
        extract_params = view(params, get_param_names(comp_))
        extract_initstates = view(initstates, get_state_names(comp_))
        extract_nns = view(nns, get_nn_names(comp_))
        tmp_pas = ComponentVector(params=extract_params, initstates=extract_initstates, nns=extract_nns)
        tmp_outputs = comp_(view(outputs, model.varindices[idx_], :, :), tmp_pas; config=config_)
        outputs = cat(outputs, tmp_outputs, dims=1)
    end
    return view(outputs, model.outputindices, :, :)
end