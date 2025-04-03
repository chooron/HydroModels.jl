"""
    get(cpt::AbstractComponent)::Num

Get the name of a component from its type parameters.
"""
get_name(cpt::AbstractComponent)::Num = typeof(cpt).parameters[1]

"""
    get_inputs(cpt::AbstractComponent)::AbstractVector{Num}

Get input variable names from component metadata.
"""
get_inputs(cpt::AbstractComponent)::AbstractVector{Num} = haskey(cpt.infos, :inputs) ? cpt.infos.inputs : Num[]
get_input_names(cpt::AbstractComponent) = length(get_inputs(cpt)) == 0 ? Symbol[] : tosymbol.(get_inputs(cpt))

"""
    get_outputs(cpt::AbstractComponent)::AbstractVector{Num}

Get output variable names from component metadata.
"""
get_outputs(cpt::AbstractComponent)::AbstractVector{Num} = haskey(cpt.infos, :outputs) ? cpt.infos.outputs : Num[]
get_output_names(cpt::AbstractComponent) = length(get_outputs(cpt)) == 0 ? Symbol[] : tosymbol.(get_outputs(cpt))

"""
    get_states(cpt::AbstractComponent)::AbstractVector{Num}

Get state variable names from component metadata.
"""
get_states(cpt::AbstractComponent)::AbstractVector{Num} = haskey(cpt.infos, :states) ? cpt.infos.states : Num[]
get_state_names(cpt::AbstractComponent) = length(get_states(cpt)) == 0 ? Symbol[] : tosymbol.(get_states(cpt))

"""
    get_params(cpt::AbstractComponent)::AbstractVector{Num}

Get parameter names from component metadata.
"""
get_params(cpt::AbstractComponent)::AbstractVector{Num} = haskey(cpt.infos, :params) ? cpt.infos.params : Num[]
get_param_names(cpt::AbstractComponent) = length(get_params(cpt)) == 0 ? Symbol[] : tosymbol.(get_params(cpt))

"""
    get_nns(cpt::AbstractComponent)::AbstractVector{Num}

Get neural network names from component metadata.
"""
get_nns(cpt::AbstractComponent) = haskey(cpt.infos, :nns) ? cpt.infos.nns : Num[]
get_nn_names(cpt::AbstractComponent) = length(get_nns(cpt)) == 0 ? Symbol[] : tosymbol.(unwrap.(get_nns(cpt)))

"""
    get_exprs(cpt::AbstractComponent)

Get expressions defined in component.
"""
get_exprs(cpt::AbstractComponent) = cpt.exprs

"""
    get_vars(comps::AbstractComponent)

Get all variable names (inputs, outputs, states) from a component.
Returns (inputs, outputs, states).
"""
get_vars(comps::AbstractComponent) = get_inputs(comps), get_outputs(comps), get_states(comps)
get_var_names(comps::AbstractComponent) = get_input_names(comps), get_output_names(comps), get_state_names(comps)

"""
    get_var_names(funcs::Vector{<:AbstractFlux}, dfuncs::Vector{<:AbstractStateFlux})

Get all variable names from flux and state flux functions.
Returns (input_names, output_names, state_names).
"""
function get_vars(funcs::Vector{<:AbstractFlux}, dfuncs::Vector{<:AbstractStateFlux})
    inputs = Vector{Num}()
    outputs = Vector{Num}()
    states = reduce(union, get_states.(dfuncs))
    for func in vcat(funcs, dfuncs)
        tmp_inputs = setdiff(setdiff(get_inputs(func), outputs), states)
        tmp_outputs = setdiff(get_outputs(func), inputs)
        union!(inputs, tmp_inputs)
        union!(outputs, tmp_outputs)
    end
    inputs, outputs, states
end

"""
    get_vars(components::Vector{<:AbstractComponent})

Get all variable names from components, handling dependencies.
Returns (input_names, output_names, state_names).
"""
function get_vars(components::Vector{<:AbstractComponent})
    inputs = Vector{Num}()
    outputs = Vector{Num}()
    states = Vector{Num}()
    for comp in components
        tmp_inputs, tmp_outputs, tmp_states = get_vars(comp)
        tmp_inputs = setdiff(tmp_inputs, outputs)
        tmp_inputs = setdiff(tmp_inputs, states)
        union!(inputs, tmp_inputs)
        union!(outputs, tmp_outputs)
        union!(states, tmp_states)
    end
    inputs, outputs, states
end
