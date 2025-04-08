"""
    EstimateComponentParams <: AbstractHydroWrapper

A wrapper component that estimates parameters for a hydrological component during simulation.

# Fields
- `component::AbstractComponent`: The wrapped hydrological component
- `estfuncs::AbstractVector{<:AbstractHydroFlux}`: Vector of estimation functions
- `ptypes::AbstractVector{Symbol}`: Parameter types to be estimated
- `meta::HydroMeta`: Metadata about the component

# Constructors
    EstimateComponentParams(component::AbstractComponent, estfuncs::AbstractVector{<:AbstractHydroFlux}, ptypes::AbstractVector{Symbol}=Symbol[])

# Description
`EstimateComponentParams` wraps a hydrological component and dynamically estimates its parameters
using provided estimation functions during simulation.

The wrapper requires at least one estimation function. Each estimation function should implement:
- `get_input_names`: Returns required input variable names
- `get_param_names`: Returns parameter names it can estimate
- `get_var_names`: Returns both input and output variable names

When called, it:
1. Applies each estimation function to calculate parameters
2. Updates the parameter set with estimated values
3. Runs the wrapped component with updated parameters
4. Returns the component output

The estimation process uses the parameter types specified in `ptypes` to organize
and update the correct parameter groups in the component.
"""
struct EstimateComponentParams{N,M<:ComponentVector} <: AbstractHydroWrapper
    component::N
    estfuncs::AbstractVector{<:AbstractHydroFlux}
    pkeys::Union{AbstractVector{Symbol},Nothing}
    meta::M

    function EstimateComponentParams(
        component::C, estfuncs::AbstractVector{E}; pkeys::P=nothing, name::Union{Symbol,Nothing}=nothing
    ) where {C,E<:AbstractHydroFlux,P<:Union{AbstractVector{Symbol},Nothing}}
        @assert length(estfuncs) != 0 "At least one estimation function is required"
        comp_meta = component.meta
        est_raw_inputs = reduce(union, get_input_vars.(estfuncs))
        est_outputs = reduce(union, get_output_vars.(estfuncs))
        est_inputs = setdiff(est_raw_inputs, est_outputs)
        est_params = reduce(union, get_param_vars.(estfuncs))
        new_params = setdiff(reduce(union, [comp_meta.params, est_inputs, est_params]), est_outputs)
        new_meta = ComponentVector(inputs=comp_meta.inputs, params=new_params, outputs=comp_meta.outputs, states=comp_meta.states, nns=comp_meta.nns)
        name = isnothing(name) ? Symbol("##wrapper#", hash(new_meta)) : name
        return new{name,typeof(new_meta)}(component, estfuncs, pkeys, new_meta)
    end
end

function (wrapper::EstimateComponentParams)(input::AbstractArray{T,2}, pas::ComponentVector; kwargs...) where {T}
    est_params_ntp = reduce(merge, map(wrapper.estfuncs) do estfunc
        tmp_input = collect(pas[:params][get_input_names(estfunc)])
        tmp_pas = collect(pas[:params][get_param_names(estfunc)])
        NamedTuple{Tuple(get_output_names(estfunc))}(estfunc(tmp_input, pas))
    end)
    new_pas = update_ca(pas, ComponentVector(params=est_params_ntp))
    output = wrapper.component(input, new_pas; kwargs...)
    return output
end

function (wrapper::EstimateComponentParams)(input::AbstractArray{T,3}, pas::ComponentVector; kwargs...) where {T}
    est_config = (ptypes=wrapper.pkeys,)
    est_params_mat = reduce(vcat, map(wrapper.estfuncs) do estfunc
        tmp_input_mat = reduce(vcat, map(wrapper.pkeys) do pkey
            tmp_input_vec = collect([pas[:params][pkey][nm] for nm in get_input_names(estfunc)])
            reshape(tmp_input_vec, 1, length(tmp_input_vec))
        end)
        tmp_input_arr = permutedims(reshape(tmp_input_mat, 1, size(tmp_input_mat)...), (3, 2, 1))
        tmp_output_mat = estfunc(tmp_input_arr, pas, config=est_config)[:, :, 1]
        tmp_output_mat
    end)
    est_param_names = reduce(union, get_output_names.(wrapper.estfuncs))
    est_params_ntp = NamedTuple{Tuple(wrapper.pkeys)}([NamedTuple{Tuple(est_param_names)}(est_params_mat[:, i]) for i in eachindex(wrapper.pkeys)])
    new_pas = update_ca(pas, ComponentVector(params=est_params_ntp))
    output = wrapper.component(input, new_pas; kwargs...)
    return output
end