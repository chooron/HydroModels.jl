"""
    ComputeComponentOutlet <: AbstractHydroWrapper

A wrapper component that computes the outlet values of a hydrological component during simulation.

# Fields
- `component::AbstractComponent`: The wrapped hydrological component
- `outlet::Symbol`: The name of the outlet node to compute values for
- `meta::HydroMeta`: Metadata about the component

# Constructors
    ComputeComponentOutlet(component::AbstractComponent, outlet::Symbol, meta::HydroMeta)

# Description
`ComputeComponentOutlet` wraps a hydrological component and extracts the values at a specific outlet node
during simulation.

The wrapper requires that the outlet name exists in the component's node names. When called, it:
1. Runs the wrapped component simulation
2. Extracts the values at the specified outlet node
3. Returns just the outlet values, either as a named tuple or array depending on configuration

This is useful for focusing analysis on a particular outlet point in a hydrological network.
"""
struct SelectComponentOutlet{N} <: AbstractDataPostprocessor
    outlet::Int
    infos::NamedTuple

    function SelectComponentOutlet(component::C, outlet::I; name::Union{Symbol,Nothing}=nothing) where {C,I<:Int}
        wrapper_name = isnothing(name) ? Symbol("##select-outlet#", get_name(component)) : name
        new{wrapper_name}(outlet, component.infos)
    end
end

function (wrapper::SelectComponentOutlet)(output::AbstractArray)
    return output[:, wrapper.outlet, :]
end

# """
#     ComputeComponentOutlet <: AbstractHydroWrapper

# A wrapper component that computes the outlet values of a hydrological component during simulation.

# # Fields
# - `component::AbstractComponent`: The wrapped hydrological component
# - `outlet::Symbol`: The name of the outlet node to compute values for
# - `meta::HydroMeta`: Metadata about the component

# # Constructor
#     ComputeComponentOutlet(component::AbstractComponent, outlet::Symbol)

# # Description
# `ComputeComponentOutlet` wraps a hydrological component and extracts the values at a specific outlet node
# during simulation.

# The wrapper requires that the outlet name exists in the component's node names. When called, it:
# 1. Runs the wrapped component simulation
# 2. Extracts the values at the specified outlet node
# 3. Returns just the outlet values, either as a named tuple or array depending on configuration

# This is useful for focusing analysis on a particular outlet point in a hydrological network.
# """
# struct WeightSumComponentOutlet{N, M<:ComponentVector} <: AbstractHydroWrapper
#     component::AbstractComponent
#     vars::AbstractVector{<:Number}
#     meta::M
#     weight_names::Symbol

#     function WeightSumComponentOutlet(component::C, vars::AbstractVector{T}; weight_names::Symbol=:weight) where {C<:AbstractComponent,T<:Num}
#         weight_ps = first(@parameters $(weight_names))
#         new_params = vcat(component.meta.params, [weight_ps])
#         var_names = Symbolics.tosymbol.(vars)
#         new_meta = ComponentVector(inputs=component.meta.inputs, params=new_params, outputs=var_names, states=component.meta.states, nns=component.meta.nns)
#         name = isnothing(name) ? Symbol("##wrapper#", hash(new_meta)) : name
#         new{name,typeof(new_meta)}(component, vars, new_meta, weight_names)
#     end
# end

# function (wrapper::WeightSumComponentOutlet)(input::AbstractArray{<:Number,3}, pas::ComponentVector; kwargs...)
#     output = wrapper.component(input, pas; kwargs...)
#     var_names = Symbolics.tosymbol.(wrapper.vars)
#     weight_vec = [pas[:params][ptype][wrapper.weight_names] for ptype in wrapper.ptypes]

#     var_idx = [findfirst(nm -> nm == var_nm, wrapper.component.varnames) for var_nm in var_names]
#     var_vec = map(var_idx) do idx
#         sum([output[idx, i, :] .* weight_vec[i] for i in eachindex(weight_vec)])
#     end
#     return reduce(hcat, var_vec)
# end