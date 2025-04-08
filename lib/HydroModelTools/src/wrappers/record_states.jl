"""
	RecordComponentState <: AbstractHydroWrapper

A wrapper component that records the state values of a hydrological component during simulation.

# Fields
- `component::AbstractComponent`: The wrapped hydrological component
- `states::ComponentVector`: Current state values of the component
- `meta::HydroMeta`: Metadata about the component

# Constructors
	RecordComponentState(component::AbstractComponent, initstates::Union{Nothing,ComponentVector}=nothing)

# Description
`RecordComponentState` wraps a hydrological component and tracks its state values during simulation.
It updates and stores the latest state values after each simulation step.

The wrapper requires that the wrapped component has state variables. If no initial states are provided,
it will use default states from the component.

When called, it:
1. Updates the input parameters with current state values
2. Runs the wrapped component simulation
3. Extracts and stores the final state values
4. Returns the original component output

The state values can be accessed through the `states` field at any time.
"""
struct RecordComponentState{N, M <: ComponentVector} <: AbstractHydroWrapper
	component::AbstractComponent
	states::ComponentVector
	meta::M

	function RecordComponentState(component::C, initstates::Union{Nothing, ComponentVector} = nothing, name::Union{Symbol,Nothing}=nothing) where {C <: AbstractComponent}
		@assert length(get_state_names(component)) != 0 "Component $component must have state flux"
		initstates = isnothing(initstates) ? get_default_states(component) : initstates
		name = isnothing(name) ? Symbol("##wrapper#", hash(initstates)) : name
		new{name,typeof(component.meta)}(component, initstates, component.meta)
	end
end

function (wrapper::RecordComponentState)(input::Any, pas::ComponentVector; kwargs...)
	reset_state = get(kwargs, :reset_state, false)
	if reset_state
		wrapper.states = pas[:initstates]
	end
	pas_replace_states = update_ca(pas, ComponentVector(initstates = wrapper.states))
	output = wrapper.component(input, pas_replace_states; kwargs...)
	state_names = get_state_names(wrapper.component)
	if get(kwargs, :convert_to_ntp, false)
		component_nodes = get_node_names(wrapper.component)
		if isnothing(component_nodes)
			states_ca = ComponentVector(NamedTuple{Tuple(state_names)}([output[nm][end] for nm in state_names]))
		else
			states_ca = ComponentVector(NamedTuple{Tuple(component_nodes)}(map(component_nodes) do node
				NamedTuple{Tuple(state_names)}([output[node][nm][end] for nm in state_names])
			end))
		end
	else
		if output isa AbstractArray{<:Number, 2}
			var_names = wrapper.component.varnames
			state_indices = map(state_names) do sname
				findfirst(vname -> vname == sname, var_names)
			end
			latest_states = output[state_indices, end]
			states_ca = ComponentVector(NamedTuple{Tuple(state_names)}(latest_states))
		elseif output isa AbstractArray{<:Number, 3}
			latest_states = output[state_indices, :, end]
			component_nodes = get_node_names(wrapper.component)
			states_ca = ComponentVector(NamedTuple{Tuple(component_nodes)}(map(component_nodes) do node
				NamedTuple{Tuple(state_names)}(latest_states[:, node])
			end
			))
		end
	end
	@reset wrapper.states = states_ca
	return output
end
