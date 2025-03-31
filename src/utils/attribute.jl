"""
    get_name(cpt::AbstractComponent)::Symbol

Get the name of a component from its type parameters.
"""
get_name(cpt::AbstractComponent)::Symbol = typeof(cpt).parameters[1]

"""
    get_input_names(cpt::AbstractComponent)::AbstractVector{Symbol}

Get input variable names from component metadata.
"""
get_input_names(cpt::AbstractComponent)::AbstractVector{Symbol} = haskey(cpt.infos, :inputs) ? cpt.infos.inputs : Symbol[]

"""
    get_output_names(cpt::AbstractComponent)::AbstractVector{Symbol}

Get output variable names from component metadata.
"""
get_output_names(cpt::AbstractComponent)::AbstractVector{Symbol} = haskey(cpt.infos, :outputs) ? cpt.infos.outputs : Symbol[]

"""
    get_state_names(cpt::AbstractComponent)::AbstractVector{Symbol}

Get state variable names from component metadata.
"""
get_state_names(cpt::AbstractComponent)::AbstractVector{Symbol} = get(cpt.infos, :states, Symbol[])

"""
    get_param_names(cpt::AbstractComponent)::AbstractVector{Symbol}

Get parameter names from component metadata.
"""
get_param_names(cpt::AbstractComponent) = get(cpt.infos, :params, Symbol[])

"""
    get_nn_names(cpt::AbstractComponent)::AbstractVector{Symbol}

Get neural network names from component metadata.
"""
get_nn_names(cpt::AbstractComponent)::AbstractVector{Symbol} = get(cpt.infos, :nns, Symbol[])

"""
    get_exprs(cpt::AbstractComponent)

Get expressions defined in component.
"""
get_exprs(cpt::AbstractComponent) = cpt.exprs

"""
    get_var_names(comps::AbstractComponent)

Get all variable names (inputs, outputs, states) from a component.
Returns (input_names, output_names, state_names).
"""
get_var_names(comps::AbstractComponent) = get_input_names(comps), get_output_names(comps), get_state_names(comps)

"""
    get_var_names(funcs::Vector{<:AbstractHydroFlux})

Get all variable names from flux functions, handling input/output dependencies.
Returns (input_names, output_names).
"""
function get_var_names(funcs::Vector{<:AbstractHydroFlux})
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    for func in funcs
        tmp_input_names = setdiff(get_input_names(func), output_names)
        tmp_output_names = setdiff(get_output_names(func), input_names)
        union!(input_names, tmp_input_names)
        union!(output_names, tmp_output_names)
    end
    input_names, output_names
end

"""
    get_var_names(funcs::Vector{<:AbstractFlux}, dfuncs::Vector{<:AbstractStateFlux})

Get all variable names from flux and state flux functions.
Returns (input_names, output_names, state_names).
"""
function get_var_names(funcs::Vector{<:AbstractFlux}, dfuncs::Vector{<:AbstractStateFlux})
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    state_names = reduce(union, get_state_names.(dfuncs))
    for func in vcat(funcs, dfuncs)
        tmp_input_names = setdiff(setdiff(get_input_names(func), output_names), state_names)
        tmp_output_names = setdiff(get_output_names(func), input_names)
        union!(input_names, tmp_input_names)
        union!(output_names, tmp_output_names)
    end
    input_names, output_names, state_names
end

"""
    get_var_names(components::Vector{<:AbstractComponent})

Get all variable names from components, handling dependencies.
Returns (input_names, output_names, state_names).
"""
function get_var_names(components::Vector{<:AbstractComponent})
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    state_names = Vector{Symbol}()
    for comp in components
        tmp_input_names = get_input_names(comp)
        tmp_output_names = get_output_names(comp)
        tmp_state_names = get_state_names(comp)
        tmp_input_names = setdiff(tmp_input_names, output_names)
        tmp_input_names = setdiff(tmp_input_names, state_names)
        union!(input_names, tmp_input_names)
        union!(output_names, tmp_output_names)
        union!(state_names, tmp_state_names)
    end
    input_names, output_names, state_names
end