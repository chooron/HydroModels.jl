"""
Neural component interfaces for HydroModels.

Lux-backed construction and state initialization are provided by the
`HydroModelsLuxExt` package extension. The core package keeps the public
component types and macros so the user-facing API remains stable while Lux is
loaded on demand.
"""

"""
    _require_lux_extension(feature)

Raise a consistent error for Lux-backed features when the Lux extension has not
been loaded.
"""
function _require_lux_extension(feature::AbstractString)
    throw(ArgumentError(
        "$feature requires the HydroModelsLuxExt package extension. " *
        "Load Lux before using it:\n\nusing HydroModels\nusing Lux"
    ))
end

# ============================================================================
# NeuralFlux - Neural Network Flux Component
# ============================================================================

"""
    NeuralFlux{C, CF, NF, NT} <: AbstractNeuralFlux

Represents a flux component driven by a neural network or neural-style callable.

The Lux-backed constructor is added by `HydroModelsLuxExt`. A pure functional
constructor remains available in core for lightweight testing or custom
integration code.

$(FIELDS)
"""
struct NeuralFlux{C,CF,NF,NT} <: AbstractNeuralFlux
    "neural flux name"
    name::Symbol
    "wrapped neural object"
    chain::C
    "compiled function that calculates the flux"
    chain_func::CF
    "input normalization function"
    norm_func::NF
    "metadata about inputs, outputs, and neural network names"
    infos::NT
end

"""
    NeuralFlux(inputs, outputs, chain; kwargs...)

Construct a Lux-backed `NeuralFlux`.

This method is provided by the `HydroModelsLuxExt` extension. Load `Lux` before
calling it.
"""
function NeuralFlux(
    inputs::Vector{T},
    outputs::Vector{T},
    chain;
    norm::Function=identity,
    name::Optional{Symbol}=nothing,
    kwargs...,
) where {T<:Num}
    _require_lux_extension("NeuralFlux construction from a Lux layer")
end

"""
    NeuralFlux(func; inputs, outputs, name=nothing, norm=identity)

Construct a `NeuralFlux` directly from a Julia function. The wrapped function is
called as `func(x)` and does not require neural-network parameters.
"""
function NeuralFlux(
    func::Function;
    inputs::Vector{Symbol},
    outputs::Vector{Symbol},
    name::Optional{Symbol}=nothing,
    norm::Function=identity,
)
    infos = HydroInfos(
        inputs=inputs,
        outputs=outputs,
        nns=Symbol[],
    )
    flux_name = isnothing(name) ? Symbol("##neural_flux_func#", hash(infos)) : name
    wrapped_func = (x, _) -> func(x)

    return NeuralFlux(
        flux_name,
        nothing,
        wrapped_func,
        norm,
        infos,
    )
end

"""
    @neuralflux [name] eq

Macro to conveniently create a `NeuralFlux` from an equation.

When the right-hand side is a Lux layer call, the actual constructor comes from
`HydroModelsLuxExt`, which is loaded automatically after `using Lux`.
"""
macro neuralflux(args...)
    name = length(args) == 1 ? nothing : args[1]
    eq_expr = length(args) == 1 ? args[1] : args[2]

    for var_name in extract_variables(eq_expr)
        if !@isdefined(var_name)
            expr_str = string(eq_expr)
            return :(error("Undefined variable '", $(string(var_name)), "' detected in expression: `", $expr_str, "`"))
        end
    end

    @assert eq_expr.head == :call && eq_expr.args[1] == :~ "Expected equation in the form: outputs ~ chain(inputs)"
    lhs, rhs = eq_expr.args[2], eq_expr.args[3]
    @assert rhs.head == :call "The right-hand side of `~` must be a function call, e.g., my_chain([x, y])"
    @assert length(rhs.args) >= 2 "The chain call must have at least one argument for the inputs"

    chain_expr, inputs_expr = rhs.args[1], rhs.args[2]

    return esc(quote
        local outputs = $lhs isa AbstractVector ? $lhs : [$lhs]
        NeuralFlux($inputs_expr, outputs, $chain_expr; name=$(name))
    end)
end

@inline function _get_neural_flux_params(flux::NeuralFlux, params::ComponentVector)
    nn_names = get_nn_names(flux)
    isempty(nn_names) && return nothing

    haskey(params, :nns) ||
        throw(ArgumentError("Missing `nns` parameters for neural flux $(flux.name)"))

    return params[:nns][only(nn_names)]
end

# 2D computation (single node)
function (flux::NeuralFlux)(
    input::AbstractArray{T,2},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...,
)::AbstractArray{T,2} where {T}
    params_cv = _as_componentvector(params)
    nn_params = _get_neural_flux_params(flux, params_cv)
    flux.chain_func(flux.norm_func(input), nn_params)
end

# 3D computation (multi-node)
function (flux::NeuralFlux)(
    input::AbstractArray{T,3},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...,
)::AbstractArray{T,3} where {T}
    params_cv = _as_componentvector(params)
    nn_params = _get_neural_flux_params(flux, params_cv)
    norm_input = flux.norm_func(input)

    node_outputs = map(1:size(input, 2)) do i
        flux.chain_func(norm_input[:, i, :], nn_params)
    end

    stack(node_outputs, dims=2)
end

# ============================================================================
# NeuralBucket - Neural Network Hydrological Bucket
# ============================================================================

"""
    NeuralBucket{FN, SN, ON, HT, I} <: AbstractHydroBucket

A neural-style hydrological bucket with three callables:
- a flux network
- a state update network
- an output network

The default Lux-based builders and network-state initialization are provided by
`HydroModelsLuxExt`.

$(FIELDS)
"""
struct NeuralBucket{FN,SN,ON,HT,I} <: AbstractHydroBucket
    "bucket name"
    name::Symbol
    "network for computing fluxes from states and inputs"
    flux_network::FN
    "network for updating states from previous states and fluxes"
    state_network::SN
    "network for computing outputs from fluxes"
    output_network::ON
    "number of input variables"
    n_inputs::Int
    "number of state variables"
    n_states::Int
    "number of output variables"
    n_outputs::Int
    "HRU types (Nothing = 2D, Vector{Int} = 3D)"
    htypes::HT
    "metadata about inputs, outputs, states, and neural network names"
    infos::I
end

"""
    NeuralBucket(; name, flux_network, state_network, output_network, ...)

Construct a `NeuralBucket` from pre-built network-like callables.
"""
function NeuralBucket(;
    name::Symbol,
    flux_network,
    state_network,
    output_network,
    n_inputs::Int,
    n_states::Int,
    n_outputs::Int,
    inputs::Vector{Symbol},
    states::Vector{Symbol},
    outputs::Vector{Symbol},
    htypes::Optional{Vector{Int}}=nothing,
)
    @assert length(inputs) == n_inputs "Number of input names must match n_inputs"
    @assert length(states) == n_states "Number of state names must match n_states"
    @assert length(outputs) == n_outputs "Number of output names must match n_outputs"

    flux_name = hasproperty(flux_network, :name) ? getproperty(flux_network, :name) : :flux_net
    state_name = hasproperty(state_network, :name) ? getproperty(state_network, :name) : :state_net
    output_name = hasproperty(output_network, :name) ? getproperty(output_network, :name) : :output_net

    infos = HydroInfos(
        inputs=inputs,
        states=states,
        outputs=outputs,
        nns=[flux_name, state_name, output_name],
    )

    return NeuralBucket(
        name,
        flux_network,
        state_network,
        output_network,
        n_inputs,
        n_states,
        n_outputs,
        htypes,
        infos,
    )
end

"""
    _initial_neural_bucket_states(bucket, rng)

Return backend-specific network state containers for a `NeuralBucket`.
"""
function _initial_neural_bucket_states(bucket::NeuralBucket, rng)
    _require_lux_extension("NeuralBucket execution")
end

"""
    _neural_bucket_step(bucket, x_t, hydro_state, nn_params, network_states)

Single timestep computation for `NeuralBucket`. Returns
`(output, new_hydro_state, new_network_states)`.
"""
function _neural_bucket_step(bucket::NeuralBucket, x_t, hydro_state, nn_params, network_states)
    flux_input = vcat(hydro_state, x_t)
    fluxes, flux_state_new = bucket.flux_network(flux_input, nn_params.flux, network_states.flux)

    state_input = vcat(hydro_state, fluxes)
    state_delta, state_state_new = bucket.state_network(state_input, nn_params.state, network_states.state)
    new_hydro_state = hydro_state .+ state_delta

    y_t, output_state_new = bucket.output_network(fluxes, nn_params.output, network_states.output)

    new_network_states = (
        flux=flux_state_new,
        state=state_state_new,
        output=output_state_new,
    )

    return y_t, new_hydro_state, new_network_states
end

@inline function _get_neural_bucket_params(bucket::NeuralBucket, params::ComponentVector)
    haskey(params, :nns) ||
        throw(ArgumentError("Missing `nns` parameters for neural bucket $(bucket.name)"))

    flux_name, state_name, output_name = bucket.infos.nns
    return (
        flux=params[:nns][flux_name],
        state=params[:nns][state_name],
        output=params[:nns][output_name],
    )
end

# 2D computation (single node, htypes = Nothing)
function (bucket::NeuralBucket{FN,SN,ON,Nothing,I})(
    input::AbstractArray{T,2},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...,
)::AbstractArray{T,2} where {FN,SN,ON,I,T}
    params_cv = _as_componentvector(params)
    n_inputs, n_steps = size(input)
    @assert n_inputs == bucket.n_inputs "Input size mismatch: expected $(bucket.n_inputs), got $n_inputs"

    nn_params = _get_neural_bucket_params(bucket, params_cv)
    network_states = _initial_neural_bucket_states(bucket, Random.default_rng())

    initstates = get(kwargs, :initstates, zeros(T, bucket.n_states))
    hydro_state = T.(Vector(initstates))

    all_states = Vector{AbstractVector{T}}(undef, n_steps)
    all_outputs = Vector{AbstractVector{T}}(undef, n_steps)

    for t in 1:n_steps
        x_t = input[:, t]
        y_t, hydro_state, network_states = _neural_bucket_step(
            bucket,
            x_t,
            hydro_state,
            nn_params,
            network_states,
        )
        all_states[t] = hydro_state
        all_outputs[t] = y_t
    end

    states_matrix = reduce(hcat, all_states)
    outputs_matrix = reduce(hcat, all_outputs)
    vcat(states_matrix, outputs_matrix)
end

# 3D computation (multi-node, htypes = Vector{Int})
function (bucket::NeuralBucket{FN,SN,ON,Vector{Int},I})(
    input::AbstractArray{T,3},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...,
)::AbstractArray{T,3} where {FN,SN,ON,I,T}
    params_cv = _as_componentvector(params)
    n_inputs, n_nodes, _ = size(input)
    @assert n_inputs == bucket.n_inputs "Input size mismatch: expected $(bucket.n_inputs), got $n_inputs"

    initstates_kw = get(kwargs, :initstates, nothing)
    node_outputs = map(1:n_nodes) do node_idx
        node_initstates = if !isnothing(initstates_kw)
            state_vals = [initstates_kw[s][node_idx] for s in get_state_names(bucket)]
            reduce(vcat, state_vals)
        else
            zeros(T, bucket.n_states)
        end

        bucket(input[:, node_idx, :], params_cv, config; initstates=node_initstates)
    end

    stack(node_outputs, dims=2)
end

# Error: single-node NeuralBucket receiving 3D input
function (bucket::NeuralBucket{FN,SN,ON,Nothing,I})(
    input::AbstractArray{T,3},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...,
) where {FN,SN,ON,I,T}
    error(
        "NeuralBucket without htypes only accepts 2D input (variables x time).\n" *
        "For multi-node computation, provide htypes.\n" *
        "Got input shape: $(size(input))",
    )
end

# Error: multi-node NeuralBucket receiving 2D input
function (bucket::NeuralBucket{FN,SN,ON,Vector{Int},I})(
    input::AbstractArray{T,2},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...,
) where {FN,SN,ON,I,T}
    error(
        "NeuralBucket with htypes only accepts 3D input (variables x nodes x time).\n" *
        "For single-node computation, omit htypes.\n" *
        "Got input shape: $(size(input))",
    )
end

"""
    create_neural_bucket(; kwargs...)

Convenience constructor for Lux-backed `NeuralBucket` components. Implemented by
`HydroModelsLuxExt`.
"""
function create_neural_bucket(; kwargs...)
    ext = Base.get_extension(@__MODULE__, :HydroModelsLuxExt)
    if !isnothing(ext)
        return create_neural_bucket(Val(:lux); kwargs...)
    end
    _require_lux_extension("create_neural_bucket")
end

"""
    create_simple_neural_bucket(; kwargs...)

Minimal Lux-backed `NeuralBucket` constructor. Implemented by
`HydroModelsLuxExt`.
"""
function create_simple_neural_bucket(; kwargs...)
    ext = Base.get_extension(@__MODULE__, :HydroModelsLuxExt)
    if !isnothing(ext)
        return create_simple_neural_bucket(Val(:lux); kwargs...)
    end
    _require_lux_extension("create_simple_neural_bucket")
end

export NeuralFlux, NeuralBucket
export @neuralflux
export create_neural_bucket, create_simple_neural_bucket
