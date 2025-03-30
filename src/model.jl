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
    varindices
    "output variables index for sort output variables"
    outputindices::AbstractVector{<:Integer}
    "meta data of hydrological model"
    infos::NamedTuple

    function HydroModel(;
        components::Vector{C},
        name::Union{Symbol,Nothing}=nothing,
        sort_components::Bool=false,
    ) where {C<:AbstractComponent}
        components = sort_components ? sort_components(components) : components
        input_names, output_names, state_names = get_var_names(components)
        param_names = reduce(union, get_param_names.(components))
        nn_names = reduce(union, get_nn_names.(components))
        input_idx, output_idx = _prepare_indices(components, input_names, vcat(input_names, state_names, output_names))
        infos = (;inputs=input_names, outputs=output_names, states=state_names, params=param_names, nns=nn_names)
        model_name = isnothing(name) ? Symbol("##model#", hash(infos)) : name
        new(model_name, components, input_idx, output_idx, infos)
    end
end

function _prepare_indices(components::Vector{<:AbstractComponent}, input_names::Vector{Symbol}, vcat_names::Vector{Symbol})
    input_idx, output_idx = Vector{Integer}[], Vector{Integer}()
    for component in components
        #* extract input index
        tmp_input_idx = map((nm) -> findfirst(varnm -> varnm == nm, input_names), get_input_names(component))
        push!(input_idx, tmp_input_idx)
        #* extract output index
        tmp_cpt_vcat_names = vcat(get_state_names(component), get_output_names(component))
        input_names = vcat(input_names, tmp_cpt_vcat_names)
        tmp_output_idx = map((nm) -> findfirst(varnm -> varnm == nm, vcat_names), tmp_cpt_vcat_names)
        output_idx = vcat(output_idx, tmp_output_idx)
    end
    return input_idx, output_idx
end

function (model::HydroModel)(
    input::AbstractArray{T,2},
    params::ComponentVector;
    initstates::ComponentVector=ComponentVector(),
    config::Union{NamedTuple,Vector{<:NamedTuple}}=NamedTuple(),
) where {T<:Number}
    initstates = length(initstates) == 0 ? ComponentVector(NamedTuple{Tuple(get_state_names(model))}(zeroes(T, length(get_state_names(model))))) : initstates
    comp_configs = config isa NamedTuple ? fill(config, length(model.components)) : config
    @assert length(comp_configs) == length(model.components) "component configs length must be equal to components length"
    outputs = input
    for (idx_, comp_, config_) in zip(model.varindices, model.components, comp_configs)
        tmp_outputs = comp_(outputs[idx_, :], params; initstates=initstates[get_state_names(comp_)], config=config_)
        outputs = cat(outputs, tmp_outputs, dims=1)
    end
    return view(outputs, model.outputindices, :)
end

function (model::HydroModel)(
    input::AbstractArray{T,3},
    params::PasDataType;
    initstates::ComponentVector=ComponentVector(), 
    config::Union{<:NamedTuple,Vector{<:NamedTuple}}=NamedTuple(),
    kwargs...,
) where {T<:Number}
    initstates = length(initstates) == 0 ? ComponentVector(NamedTuple{Tuple(get_state_names(model))}(fill(zeroes(T, size(input, 2)), length(get_state_names(model))))) : initstates
    comp_configs = config isa NamedTuple ? fill(config, length(model.components)) : config
    @assert length(comp_configs) == length(model.components) "component configs length must be equal to components length"
    outputs = input
    for (idx_, comp_, config_) in zip(model.varindices, model.components, comp_configs)
        tmp_outputs = comp_(view(outputs, idx_, :, :), params; initstates=initstates[get_state_names(comp_)], config=config_)
        outputs = cat(outputs, tmp_outputs, dims=1)
    end
    return view(outputs, model.outputindices, :, :)
end