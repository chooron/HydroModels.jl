"""
    get_name(cpt::AbstractComponent)::Symbol

Get the name of a component from its type parameters.
"""
get_name(cpt::AbstractComponent)::Symbol = cpt.name

"""
    get_input_vars(cpt::AbstractComponent)::AbstractVector{Num}
    get_input_names(cpt::AbstractComponent)::AbstractVector{Symbol}

Get the input variables or their names from a component's metadata.
Returns empty vector if no inputs are defined.
"""
get_input_names(cpt::AbstractComponent)::AbstractVector{Symbol} = haskey(cpt.infos, :inputs) ? cpt.infos.inputs : Symbol[]

"""
    get_output_vars(cpt::AbstractComponent)::AbstractVector{Num}
    get_output_names(cpt::AbstractComponent)::AbstractVector{Symbol}

Get the output variables or their names from a component's metadata.
Returns empty vector if no outputs are defined.
"""
get_output_names(cpt::AbstractComponent)::AbstractVector{Symbol} = haskey(cpt.infos, :outputs) ? cpt.infos.outputs : Symbol[]

"""
    get_state_vars(cpt::AbstractComponent)::AbstractVector{Num}
    get_state_names(cpt::AbstractComponent)::AbstractVector{Symbol}

Get the state variables or their names from a component's metadata.
Returns empty vector if no states are defined.
"""
get_state_names(cpt::AbstractComponent)::AbstractVector{Symbol} = get(cpt.infos, :states, Symbol[])

"""
    get_param_names(cpt::AbstractComponent)::AbstractVector{Symbol}

Get the parameter variables or their names from a component's metadata.
Returns empty tuple if no parameters are defined.
"""
get_param_names(cpt::AbstractComponent) = get(cpt.infos, :params, Symbol[])

"""
    get_nn_vars(cpt::AbstractComponent)::AbstractVector
    get_nn_names(cpt::AbstractComponent)::AbstractVector{Symbol}

Get the neural network variables or their names from a component's metadata.
Returns empty tuple/ComponentVector if no neural networks are defined.
"""
get_nn_names(cpt::AbstractComponent)::AbstractVector{Symbol} = get(cpt.infos, :nns, Symbol[])


"""
    get_exprs(cpt::AbstractComponent)

Get the expressions defined in a component.
"""
get_exprs(cpt::AbstractComponent) = cpt.exprs

"""
    get_var_names(comps::AbstractComponent)

Get all variable names (inputs, outputs, and states) from a component.
Returns a tuple of three vectors: (input_names, output_names, state_names).
"""
get_var_names(comps::AbstractComponent) = get_input_names(comps), get_output_names(comps), get_state_names(comps)

"""
    get_var_names(funcs::Vector{<:AbstractHydroFlux})

Get all variable names from a vector of flux functions, handling dependencies between inputs and outputs.
Returns a tuple of three vectors: (input_names, output_names, state_names).
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

Get all variable names from vectors of flux and state flux functions, handling dependencies between inputs, outputs, and states.
Returns a tuple of three vectors: (input_names, output_names, state_names).
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
    get_var_names(funcs::Vector{<:AbstractComponent})

Get all variable names from a vector of flux functions, handling dependencies between inputs and outputs.
Returns a tuple of three vectors: (input_names, output_names, state_names).
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


# """
# 获取demo输入数据
# """
# get_demo_pas(flux::AbstractFlux) = ComponentVector(
#     params=NamedTuple{Tuple(flux.infos[:param])}(ones(length(flux.infos[:param]))),
#     initstates=NamedTuple{Tuple(flux.infos[:state])}(ones(length(flux.infos[:state]))),
# )

# get_demo_pas(flux::AbstractNeuralFlux) = ComponentVector(
#     nn=NamedTuple{Tuple(flux.infos[:nn])}(zeros(flux.nnios[:paramlen]))
# )

# get_demo_pas(bucket::AbstractBucket) = ComponentVector(
#     params=NamedTuple{Tuple(bucket.infos[:param])}(ones(length(flux.infos[:param]))),
#     initstates=NamedTuple{Tuple(bucket.infos[:state])}(ones(length(bucket.infos[:state]))),
#     nn=NamedTuple{Tuple(bucket.infos[:nn])}([zeros(flux.nnios[:paramlen]) for flux in bucket.fluxes if flux isa AbstractNeuralFlux]),
# )

# get_demo_pas(route::AbstractRoute) = ComponentVector(
#     params=NamedTuple{Tuple(route.infos[:param])}(ones(length(route.infos[:param]))),
#     initstates=NamedTuple{Tuple(route.infos[:state])}(ones(length(route.infos[:state]))),
#     nn=NamedTuple{Tuple(route.infos[:nn])}([zeros(flux.nnios[:paramlen]) for flux in route.fluxes if flux isa AbstractNeuralFlux]),
# )