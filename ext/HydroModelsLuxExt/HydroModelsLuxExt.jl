module HydroModelsLuxExt

using ComponentArrays
using ComponentArrays: ComponentVector, getaxes
using HydroModelCore: HydroInfos
using HydroModels
using Lux
using LuxCore
using Random

function HydroModels.NeuralFlux(
    inputs::Vector{T},
    outputs::Vector{T},
    chain::LuxCore.AbstractLuxLayer;
    norm::Function=identity,
    name::HydroModels.Optional{Symbol}=nothing,
    st=LuxCore.initialstates(Random.default_rng(), chain),
    chain_name::HydroModels.Optional{Symbol}=nothing,
) where {T<:HydroModels.Num}
    default_chain_name = hasproperty(chain, :name) ? getproperty(chain, :name) : nothing
    resolved_chain_name = isnothing(chain_name) ? default_chain_name : chain_name
    @assert !isnothing(resolved_chain_name) "`chain_name` must be provided for NeuralFlux, or set `name` in chain"

    ps = LuxCore.initialparameters(Random.default_rng(), chain)
    ps_axes = getaxes(ComponentVector(ps))
    nn_func = (x, p) -> LuxCore.apply(chain, x, ComponentVector(p, ps_axes), st)[1]

    infos = HydroInfos(
        inputs=!isempty(inputs) ? HydroModels.tosymbol.(inputs) : Symbol[],
        outputs=!isempty(outputs) ? HydroModels.tosymbol.(outputs) : Symbol[],
        nns=[resolved_chain_name],
    )
    flux_name = isnothing(name) ? Symbol("##neural_flux#", hash(infos)) : name

    return HydroModels.NeuralFlux(
        flux_name,
        chain,
        nn_func,
        norm,
        infos,
    )
end

function HydroModels._initial_neural_bucket_states(
    bucket::HydroModels.NeuralBucket{FN,SN,ON},
    rng,
) where {
    FN<:LuxCore.AbstractLuxLayer,
    SN<:LuxCore.AbstractLuxLayer,
    ON<:LuxCore.AbstractLuxLayer,
}
    return (
        flux=LuxCore.initialstates(rng, bucket.flux_network),
        state=LuxCore.initialstates(rng, bucket.state_network),
        output=LuxCore.initialstates(rng, bucket.output_network),
    )
end

function HydroModels._get_nn_params(
    component::HydroModels.NeuralFlux{C},
    nn_names,
    rng,
) where {C<:LuxCore.AbstractLuxLayer}
    ps = LuxCore.initialparameters(rng, component.chain)
    return NamedTuple{Tuple(nn_names)}((ComponentVector(ps),))
end

function HydroModels._get_nn_params(
    component::HydroModels.NeuralBucket{FN,SN,ON},
    nn_names,
    rng,
) where {
    FN<:LuxCore.AbstractLuxLayer,
    SN<:LuxCore.AbstractLuxLayer,
    ON<:LuxCore.AbstractLuxLayer,
}
    flux_ps = LuxCore.initialparameters(rng, component.flux_network)
    state_ps = LuxCore.initialparameters(rng, component.state_network)
    output_ps = LuxCore.initialparameters(rng, component.output_network)

    return NamedTuple{Tuple(nn_names)}((
        ComponentVector(flux_ps),
        ComponentVector(state_ps),
        ComponentVector(output_ps),
    ))
end

function HydroModels.create_neural_bucket(
    ::Type{Val{:lux}};
    name::Symbol,
    n_inputs::Int,
    n_states::Int,
    n_outputs::Int,
    n_fluxes::Int=n_states,
    hidden_size::Int=16,
    inputs::Vector{Symbol},
    states::Vector{Symbol},
    outputs::Vector{Symbol},
    htypes::HydroModels.Optional{Vector{Int}}=nothing,
    flux_activation=tanh,
    state_activation=identity,
    output_activation=identity,
)
    flux_network = Lux.Chain(
        Lux.Dense(n_states + n_inputs => hidden_size, flux_activation),
        Lux.Dense(hidden_size => n_fluxes),
        name=Symbol(name, :_flux),
    )
    state_network = Lux.Chain(
        Lux.Dense(n_states + n_fluxes => hidden_size, state_activation),
        Lux.Dense(hidden_size => n_states),
        name=Symbol(name, :_state),
    )
    output_network = Lux.Chain(
        Lux.Dense(n_fluxes => n_outputs, output_activation),
        name=Symbol(name, :_output),
    )

    return HydroModels.NeuralBucket(
        name=name,
        flux_network=flux_network,
        state_network=state_network,
        output_network=output_network,
        n_inputs=n_inputs,
        n_states=n_states,
        n_outputs=n_outputs,
        inputs=inputs,
        states=states,
        outputs=outputs,
        htypes=htypes,
    )
end

function HydroModels.create_simple_neural_bucket(
    ::Type{Val{:lux}};
    name::Symbol,
    n_inputs::Int,
    n_states::Int,
    n_outputs::Int,
    inputs::Vector{Symbol},
    states::Vector{Symbol},
    outputs::Vector{Symbol},
    htypes::HydroModels.Optional{Vector{Int}}=nothing,
)
    flux_network = Lux.Dense(
        n_states + n_inputs => n_states,
        name=Symbol(name, :_flux),
    )
    state_network = Lux.Dense(
        n_states + n_states => n_states,
        name=Symbol(name, :_state),
    )
    output_network = Lux.Dense(
        n_states => n_outputs,
        name=Symbol(name, :_output),
    )

    return HydroModels.NeuralBucket(
        name=name,
        flux_network=flux_network,
        state_network=state_network,
        output_network=output_network,
        n_inputs=n_inputs,
        n_states=n_states,
        n_outputs=n_outputs,
        inputs=inputs,
        states=states,
        outputs=outputs,
        htypes=htypes,
    )
end

end # module
