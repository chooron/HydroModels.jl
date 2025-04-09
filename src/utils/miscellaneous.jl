"""
    sort_fluxes(fluxes::AbstractVector{<:AbstractComponent})

Construct a directed calculation graph based on hydrological fluxes and return them in topological order.

# Arguments
- `fluxes::AbstractVector{<:AbstractComponent}`: A vector of flux components to be sorted.

# Returns
- `AbstractVector{<:AbstractComponent}`: A vector of flux components sorted in topological order for calculation.

# Description
This function creates a directed graph representing the dependencies between flux components,
where edges connect input variables to output variables. The function then performs a topological
sort to determine the correct calculation order that respects these dependencies.

The process involves:
1. Identifying all input and output variables across all flux components
2. Building a directed graph where nodes are variables and edges represent dependencies
3. Performing a topological sort on the graph
4. Extracting the flux components in the order determined by the sort

This ensures that when calculations are performed, all required inputs are available before
a flux component is evaluated.

# Examples
```julia
fluxes = [flux1, flux2, flux3]
sorted_fluxes = sort_fluxes(fluxes)
# Now sorted_fluxes contains the components in the correct calculation order
```
"""
function sort_fluxes(fluxes::AbstractVector{<:AbstractComponent})
    input_names = reduce(union, get_input_names.(fluxes))
    output_names = reduce(union, get_output_names.(fluxes))
    input_names = setdiff(input_names, output_names)
    output_names = setdiff(output_names, input_names)

    # Build a named tuple mapping output names to their corresponding flux instances
    fluxes_ntp = reduce(merge, map(fluxes) do flux
        tmp_output_names = get_output_names(flux)
        NamedTuple{Tuple(tmp_output_names)}(repeat([flux], length(tmp_output_names)))
    end)

    # Construct a directed calculation graph
    var_names = vcat(input_names, output_names)
    var_names_ntp = NamedTuple{Tuple(var_names)}(1:length(var_names))
    digraph = SimpleDiGraph(length(var_names))
    for flux in fluxes
        tmp_input_names, tmp_output_names = get_input_names(flux), get_output_names(flux)
        for ipnm in tmp_input_names
            for opnm in tmp_output_names
                println((ipnm => opnm))
                add_edge!(digraph, var_names_ntp[ipnm], var_names_ntp[opnm])
            end
        end
    end
    
    # Sort fluxes based on the topological order of the directed graph
    sorted_fluxes = AbstractComponent[]
    for idx in topological_sort(digraph)
        tmp_var_nm = var_names[idx]
        if (tmp_var_nm in output_names)
            tmp_flux = fluxes_ntp[tmp_var_nm]
            if !(tmp_flux in sorted_fluxes)
                push!(sorted_fluxes, tmp_flux)
            end
        end
    end
    sorted_fluxes
end

"""
    sort_components(components::AbstractVector{<:AbstractComponent})

Construct a directed calculation graph based on hydrological components and return them in topological order.

# Arguments
- `components::AbstractVector{<:AbstractComponent}`: A vector of hydrological components to be sorted.

# Returns
- `AbstractVector{<:AbstractComponent}`: A vector of components sorted in topological order for calculation.

# Description
This function creates a directed graph representing the dependencies between hydrological components,
where edges connect input variables to output and state variables. The function then performs a 
topological sort to determine the correct calculation order that respects these dependencies.

The process involves:
1. Identifying all input, output, and state variables across all components
2. Building a directed graph where nodes are variables and edges represent dependencies
3. Performing a topological sort on the graph
4. Extracting the components in the order determined by the sort

This ensures that when calculations are performed, all required inputs are available before
a component is evaluated.

# Examples
```julia
components = [bucket1, flux1, route1]
sorted_components = sort_components(components)
# Now sorted_components contains the components in the correct calculation order
```
"""
function sort_components(components::AbstractVector{<:AbstractComponent})
    input_names, output_names, state_names = get_var_names(components)
    components_ntp = reduce(merge, map(components) do component
        tmp_input_names, tmp_output_names, tmp_state_names = get_var_names(component)
        tmp_output_state_names = vcat(tmp_output_names, tmp_state_names)
        NamedTuple{Tuple(tmp_output_state_names)}(repeat([component], length(tmp_output_state_names)))
    end)
    var_names = reduce(union, [input_names, output_names, state_names])
    var_names_ntp = namedtuple(var_names, collect(1:length(var_names)))
    digraph = SimpleDiGraph(length(var_names))
    for component in components
        tmp_input_names, tmp_output_names, tmp_state_names = get_var_names(component)
        tmp_output_names = vcat(tmp_output_names, tmp_state_names)
        for ipnm in tmp_input_names
            for opnm in tmp_output_names
                add_edge!(digraph, var_names_ntp[ipnm], var_names_ntp[opnm])
            end
        end
    end
    sorted_components = AbstractComponent[]
    for idx in topological_sort(digraph)
        tmp_var_nm = var_names[idx]
        if (tmp_var_nm in output_names)
            tmp_component = components_ntp[tmp_var_nm]
            if !(tmp_component in sorted_components)
                push!(sorted_components, tmp_component)
            end
        end
    end
    sorted_components
end

"""
Expand the parameters of a component vector based on the provided index.

# Arguments
- `pas::ComponentVector`: The component vector to be expanded.
- `ptyidx::AbstractVector`: The index of the parameters to be expanded.

# Returns
- `new_pas::ComponentVector`: The expanded component vector.
"""
function expand_component_params(pas::ComponentVector, ptyidx::AbstractVector)
    params = view(pas, :params)
    expand_params = NamedTuple{Tuple(keys(params))}([params[p][ptyidx] for p in keys(params)])
    return if haskey(pas, :nns)
        ComponentVector(params=expand_params, nns=pas[:nns])
    else
        ComponentVector(params=expand_params)
    end
end

"""
Expand the initial states of a component vector based on the provided index.

# Arguments
- `initstates::ComponentVector`: The component vector to be expanded.
- `styidx::AbstractVector`: The index of the initial states to be expanded.
- `num_nodes::Int`: The number of nodes.

# Returns
- `new_pas::ComponentVector`: The expanded component vector.
"""
function expand_component_initstates(initstates::ComponentVector, styidx::AbstractVector)
    num_states = length(keys(initstates))
    initstates_arr = reshape(Vector(initstates), :, num_states)'
    expand_component_initstates(initstates_arr, styidx)
end

function expand_component_initstates(initstates::AbstractMatrix, styidx::AbstractVector)
    initstates[:, styidx]
end

function get_default_states(component::AbstractComponent, input::AbstractArray{T,2}) where {T}
    state_names = get_state_names(component)
    return ComponentVector(NamedTuple{Tuple(state_names)}(fill(zero(T), length(state_names))))
end

function get_default_states(component::AbstractComponent, input::AbstractArray{T,3}) where {T}
    state_names = get_state_names(component)
    return ComponentVector(NamedTuple{Tuple(state_names)}(fill(zeros(T, size(input, 2)), length(state_names))))
end