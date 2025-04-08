"""
    get_name(cpt)::Symbol

Get the name of a component from its type parameters.
"""
get_name(cpt)::Symbol = typeof(cpt).parameters[1]

"""
    get_input_names(cpt)::AbstractVector{Symbol}

Get input variable names from component metadata.
"""
get_input_names(cpt)::AbstractVector{Symbol} = haskey(cpt.infos, :inputs) ? cpt.infos.inputs : Symbol[]

"""
    get_output_names(cpt)::AbstractVector{Symbol}

Get output variable names from component metadata.
"""
get_output_names(cpt)::AbstractVector{Symbol} = haskey(cpt.infos, :outputs) ? cpt.infos.outputs : Symbol[]

"""
    get_state_names(cpt)::AbstractVector{Symbol}

Get state variable names from component metadata.
"""
get_state_names(cpt)::AbstractVector{Symbol} = get(cpt.infos, :states, Symbol[])

"""
    get_param_names(cpt)::AbstractVector{Symbol}

Get parameter names from component metadata.
"""
get_param_names(cpt)::AbstractVector{Symbol} = get(cpt.infos, :params, Symbol[])

"""
    get_nn_names(cpt)::AbstractVector{Symbol}

Get neural network names from component metadata.
"""
get_nn_names(cpt)::AbstractVector{Symbol} = get(cpt.infos, :nns, Symbol[])