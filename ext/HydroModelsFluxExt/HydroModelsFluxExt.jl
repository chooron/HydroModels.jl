module HydroModelsFluxExt

using ComponentArrays: ComponentVector
using Flux
using HydroModels
using HydroModels: HydroInfos

"""
    FluxNetworkAdapter

Stores a Flux model, its `destructure` rebuild closure, and the model's
current flat parameter vector. The adapter is the backend boundary used by
`NeuralBucket`; its three-argument call matches HydroModels' network protocol.
"""
struct FluxNetworkAdapter{M,R,P}
    name::Symbol
    model::M
    rebuild::R
    parameters::P
end

function FluxNetworkAdapter(model, name::Symbol)
    _validate_flux_model(model)
    parameters, rebuild = Flux.destructure(model)
    return FluxNetworkAdapter(name, model, rebuild, copy(parameters))
end

function _validate_flux_model(model)
    model isa Flux.Chain ||
        model isa Flux.Dense ||
        throw(
            ArgumentError(
                "HydroModels Flux backend supports Flux.Chain and Flux.Dense models " *
                "in the forward-only stateless subset",
            ),
        )
    _validate_flux_layer(model)
    return nothing
end

function _validate_flux_layer(layer)
    layer_name = nameof(typeof(layer))
    if layer_name in (:BatchNorm, :Dropout, :RNNCell, :LSTMCell, :GRUCell, :Recur)
        throw(
            ArgumentError(
                "Flux layer $(layer_name) is stateful or mode-dependent and is " *
                "not supported by HydroModelsFluxExt; use a stateless Dense/Chain model",
            ),
        )
    end
    if layer isa Flux.Chain
        foreach(_validate_flux_layer, layer.layers)
    end
    return nothing
end

function (network::FluxNetworkAdapter)(x, p, ::Nothing)
    runtime_model = network.rebuild(p.params)
    return runtime_model(x), nothing
end

function _flux_network_params(network::FluxNetworkAdapter)
    return ComponentVector(; params=copy(network.parameters))
end

function HydroModels.NeuralFlux(
    inputs::Vector{T},
    outputs::Vector{T},
    chain::Union{Flux.Chain,Flux.Dense};
    norm=identity,
    name::HydroModels.Optional{Symbol}=nothing,
    chain_name::HydroModels.Optional{Symbol}=nothing,
) where {T<:HydroModels.Num}
    isnothing(chain_name) && throw(
        ArgumentError(
            "`chain_name` must be provided for Flux-backed NeuralFlux because " *
            "Flux models do not carry HydroModels network names",
        ),
    )

    network = FluxNetworkAdapter(chain, chain_name)
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
) where {A<:FluxNetworkAdapter}
    length(nn_names) == 1 ||
        throw(ArgumentError("Flux NeuralFlux expects exactly one network name"))
    return NamedTuple{Tuple(nn_names)}((_flux_network_params(component.chain),))
end

function HydroModels._initial_neural_bucket_states(
    bucket::HydroModels.NeuralBucket{FN,SN,ON}, rng
) where {FN<:FluxNetworkAdapter,SN<:FluxNetworkAdapter,ON<:FluxNetworkAdapter}
    return (flux=nothing, state=nothing, output=nothing)
end

function HydroModels._get_nn_params(
    component::HydroModels.NeuralBucket{FN,SN,ON}, nn_names, rng
) where {FN<:FluxNetworkAdapter,SN<:FluxNetworkAdapter,ON<:FluxNetworkAdapter}
    length(nn_names) == 3 || throw(
        ArgumentError("Flux NeuralBucket expects flux, state, and output network names")
    )
    return NamedTuple{Tuple(nn_names)}((
        _flux_network_params(component.flux_network),
        _flux_network_params(component.state_network),
        _flux_network_params(component.output_network),
    ))
end

function HydroModels.create_neural_bucket(
    ::Val{:flux};
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
    flux_network = FluxNetworkAdapter(
        Flux.Chain(
            Flux.Dense(n_states + n_inputs => hidden_size, flux_activation),
            Flux.Dense(hidden_size => n_fluxes),
        ),
        Symbol(name, :_flux),
    )
    state_network = FluxNetworkAdapter(
        Flux.Chain(
            Flux.Dense(n_states + n_fluxes => hidden_size, state_activation),
            Flux.Dense(hidden_size => n_states),
        ),
        Symbol(name, :_state),
    )
    output_network = FluxNetworkAdapter(
        Flux.Chain(Flux.Dense(n_fluxes => n_outputs, output_activation)),
        Symbol(name, :_output),
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

function HydroModels.create_simple_neural_bucket(
    ::Val{:flux};
    name::Symbol,
    n_inputs::Int,
    n_states::Int,
    n_outputs::Int,
    inputs::Vector{Symbol},
    states::Vector{Symbol},
    outputs::Vector{Symbol},
    htypes::HydroModels.Optional{Vector{Int}}=nothing,
)
    flux_network = FluxNetworkAdapter(
        Flux.Chain(Flux.Dense(n_states + n_inputs => n_states)), Symbol(name, :_flux)
    )
    state_network = FluxNetworkAdapter(
        Flux.Chain(Flux.Dense(n_states + n_states => n_states)), Symbol(name, :_state)
    )
    output_network = FluxNetworkAdapter(
        Flux.Chain(Flux.Dense(n_states => n_outputs)), Symbol(name, :_output)
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

end
