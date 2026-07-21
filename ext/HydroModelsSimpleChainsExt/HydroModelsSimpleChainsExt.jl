module HydroModelsSimpleChainsExt

using ComponentArrays: ComponentVector
using HydroModels
using HydroModels: HydroInfos
using SimpleChains

"""
    SimpleChainsNetworkAdapter

Wraps a stateless `SimpleChain` and exposes HydroModels' three-argument
network protocol. SimpleChains parameters remain separate from the chain
structure and are held in the `params` child of the ComponentVector.
"""
struct SimpleChainsNetworkAdapter{C}
    name::Symbol
    chain::C
end

function (network::SimpleChainsNetworkAdapter)(x, p, ::Nothing)
    return network.chain(x, p.params), nothing
end

function _simplechains_network_params(network::SimpleChainsNetworkAdapter; rng)
    # SimpleChains 0.4.8 requires an element type for deterministic parameter
    # storage. Float32 is its documented default; mixed Float64 inputs remain
    # valid because SimpleChains promotes input and parameter element types.
    ps = SimpleChains.init_params(network.chain, Float32; rng)
    return ComponentVector(; params=collect(ps))
end

function HydroModels.NeuralFlux(
    inputs::Vector{T},
    outputs::Vector{T},
    chain::SimpleChains.SimpleChain;
    norm=identity,
    name::HydroModels.Optional{Symbol}=nothing,
    chain_name::HydroModels.Optional{Symbol}=nothing,
) where {T<:HydroModels.Num}
    isnothing(chain_name) && throw(
        ArgumentError("`chain_name` must be provided for SimpleChains-backed NeuralFlux"),
    )

    network = SimpleChainsNetworkAdapter(chain_name, chain)
    infos = HydroInfos(;
        inputs=HydroModels.tosymbol.(inputs),
        outputs=HydroModels.tosymbol.(outputs),
        nns=[chain_name],
    )
    flux_name = isnothing(name) ? Symbol("##neural_flux#", hash(infos)) : name
    nn_func = (x, p) -> network(x, p, nothing)[1]

    return HydroModels.NeuralFlux(flux_name, network, nn_func, norm, infos)
end

function HydroModels._get_nn_params(
    component::HydroModels.NeuralFlux{A}, nn_names, rng
) where {A<:SimpleChainsNetworkAdapter}
    length(nn_names) == 1 ||
        throw(ArgumentError("SimpleChains NeuralFlux expects exactly one network name"))
    params = _simplechains_network_params(component.chain; rng)
    return NamedTuple{Tuple(nn_names)}((params,))
end

function HydroModels._initial_neural_bucket_states(
    bucket::HydroModels.NeuralBucket{FN,SN,ON}, rng
) where {
    FN<:SimpleChainsNetworkAdapter,
    SN<:SimpleChainsNetworkAdapter,
    ON<:SimpleChainsNetworkAdapter,
}
    return (flux=nothing, state=nothing, output=nothing)
end

function HydroModels._get_nn_params(
    component::HydroModels.NeuralBucket{FN,SN,ON}, nn_names, rng
) where {
    FN<:SimpleChainsNetworkAdapter,
    SN<:SimpleChainsNetworkAdapter,
    ON<:SimpleChainsNetworkAdapter,
}
    length(nn_names) == 3 || throw(
        ArgumentError(
            "SimpleChains NeuralBucket expects flux, state, and output network names"
        ),
    )
    return NamedTuple{Tuple(nn_names)}((
        _simplechains_network_params(component.flux_network; rng),
        _simplechains_network_params(component.state_network; rng),
        _simplechains_network_params(component.output_network; rng),
    ))
end

function _simplechains_bucket(
    name,
    n_inputs,
    n_states,
    n_outputs,
    n_fluxes,
    hidden_size,
    inputs,
    states,
    outputs,
    htypes,
    flux_activation,
    state_activation,
    output_activation,
)
    flux_network = SimpleChainsNetworkAdapter(
        Symbol(name, :_flux),
        SimpleChain(
            SimpleChains.static(n_states + n_inputs),
            TurboDense(flux_activation, hidden_size),
            TurboDense(identity, n_fluxes),
        ),
    )
    state_network = SimpleChainsNetworkAdapter(
        Symbol(name, :_state),
        SimpleChain(
            SimpleChains.static(n_states + n_fluxes),
            TurboDense(state_activation, hidden_size),
            TurboDense(identity, n_states),
        ),
    )
    output_network = SimpleChainsNetworkAdapter(
        Symbol(name, :_output),
        SimpleChain(
            SimpleChains.static(n_fluxes), TurboDense(output_activation, n_outputs)
        ),
    )

    return HydroModels.NeuralBucket(;
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

function HydroModels.create_neural_bucket(
    ::Val{:simplechains};
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
    return _simplechains_bucket(
        name,
        n_inputs,
        n_states,
        n_outputs,
        n_fluxes,
        hidden_size,
        inputs,
        states,
        outputs,
        htypes,
        flux_activation,
        state_activation,
        output_activation,
    )
end

function HydroModels.create_simple_neural_bucket(
    ::Val{:simplechains};
    name::Symbol,
    n_inputs::Int,
    n_states::Int,
    n_outputs::Int,
    inputs::Vector{Symbol},
    states::Vector{Symbol},
    outputs::Vector{Symbol},
    htypes::HydroModels.Optional{Vector{Int}}=nothing,
)
    return _simplechains_bucket(
        name,
        n_inputs,
        n_states,
        n_outputs,
        n_states,
        n_states,
        inputs,
        states,
        outputs,
        htypes,
        identity,
        identity,
        identity,
    )
end

end
