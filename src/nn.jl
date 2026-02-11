"""
Neural network components module - integrates Lux.jl neural networks with HydroModels framework.

This module provides:
- `NeuralFlux`: Neural network-driven flux component
- `NeuralBucket`: Neural network-based hydrological bucket with RNN-like recurrence structure

# Design Philosophy

The NeuralBucket leverages the structural similarity between hydrological ODEs and RNNs:
- **Hydrological ODE**: `dS/dt = f(S, input, params)` → iterate over time
- **RNN**: `h_t = f(h_{t-1}, x_t, params)` → iterate over time

NeuralBucket implements the same interface as HydroBucket (`AbstractHydroBucket`),
so it can be directly embedded into HydroModel.
"""

# ============================================================================
# NeuralFlux - Neural Network Flux Component
# ============================================================================

"""
    NeuralFlux{C, CF, NF, NT} <: AbstractNeuralFlux

Represents a flux component driven by a neural network.

It wraps a Lux.AbstractLuxLayer and connects it to symbolic variables for integration
into a hydrological model. Supports both 2D (single-node) and 3D (multi-node) input.

$(FIELDS)
"""
struct NeuralFlux{C,CF,NF,NT} <: AbstractNeuralFlux
    "neural flux name"
    name::Symbol
    "chain of the neural network"
    chain::C
    "Compiled function that calculates the flux using the neural network"
    chain_func::CF
    "input normalization functions"
    norm_func::NF
    "Information about the neural network's input and output structure"
    infos::NT

    function NeuralFlux(
        inputs::Vector{T},
        outputs::Vector{T},
        chain::LuxCore.AbstractLuxLayer;
        norm::Function=(x) -> x,
        name::Optional{Symbol}=nothing,
        st=LuxCore.initialstates(Random.default_rng(), chain),
        chain_name::Optional{Symbol}=nothing,
    ) where {T<:Num}
        chain_name = chain_name === nothing ? chain.name : chain_name
        @assert !isnothing(chain_name) "`chain_name` must be provided for NeuralFlux, or set `name` in chain"

        ps = LuxCore.initialparameters(Random.default_rng(), chain)
        ps_axes = getaxes(ComponentVector(ps))
        nn_func = (x, p) -> LuxCore.apply(chain, x, ComponentVector(p, ps_axes), st)[1]

        infos = HydroInfos(
            inputs=!isempty(inputs) ? tosymbol.(inputs) : Symbol[],
            outputs=!isempty(outputs) ? tosymbol.(outputs) : Symbol[],
            nns=[chain_name]
        )
        flux_name = isnothing(name) ? Symbol("##neural_flux#", hash(infos)) : name

        new{typeof(chain),typeof(nn_func),typeof(norm),typeof(infos)}(
            flux_name, chain, nn_func, norm, infos
        )
    end

    # Functional constructor - directly use Julia function
    function NeuralFlux(
        func::Function;
        inputs::Vector{Symbol},
        outputs::Vector{Symbol},
        name::Optional{Symbol}=nothing,
        norm::Function=(x) -> x,
    )
        infos = HydroInfos(
            inputs=inputs,
            outputs=outputs,
            nns=Symbol[]
        )
        flux_name = isnothing(name) ? Symbol("##neural_flux_func#", hash(infos)) : name

        wrapped_func = (x, p) -> func(x)

        new{Nothing,typeof(wrapped_func),typeof(norm),typeof(infos)}(
            flux_name, nothing, wrapped_func, norm, infos
        )
    end
end

# ============================================================================
# @neuralflux macro
# ============================================================================

"""
    @neuralflux [name] eq

Macro to conveniently create a NeuralFlux from an equation.

# Examples
```julia
@variables x, y, z
chain = Chain(Dense(2 => 10, relu), Dense(10 => 1), name=:my_net)
flux = @neuralflux z ~ chain([x, y])
```
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

# ============================================================================
# NeuralFlux functor methods
# ============================================================================

# NeuralFlux computation for 2D input (single node)
function (flux::NeuralFlux)(
    input::AbstractArray{T,2},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,2} where {T}
    params = _as_componentvector(params)
    nn_params = params[:nns][get_nn_names(flux)[1]]
    flux.chain_func(flux.norm_func(input), nn_params)
end

# NeuralFlux computation for 3D input (multi-node)
function (flux::NeuralFlux)(
    input::AbstractArray{T,3},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,3} where {T}
    params = _as_componentvector(params)
    nn_params = params[:nns][get_nn_names(flux)[1]]
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

A neural network-based hydrological bucket that mimics ODE structure using RNN-like recurrence.

Implements the same interface as `HydroBucket` (`AbstractHydroBucket`), so it can be
directly embedded into `HydroModel`.

# Architecture
Three neural networks:
1. **Flux Network**: `(n_states + n_inputs) → hidden → n_fluxes`
2. **State Network**: `(n_states + n_fluxes) → hidden → n_states` (outputs delta_S)
3. **Output Network**: `n_fluxes → n_outputs`

# Recurrence
```
flux_t = flux_network(vcat(S_{t-1}, input_t))
delta_S_t = state_network(vcat(S_{t-1}, flux_t))
S_t = S_{t-1} + delta_S_t  # Residual connection
output_t = output_network(flux_t)
```

$(FIELDS)
"""
struct NeuralBucket{FN,SN,ON,HT,I} <: AbstractHydroBucket
    "Bucket name"
    name::Symbol
    "Neural network for computing fluxes from states and inputs"
    flux_network::FN
    "Neural network for updating states from previous states and fluxes"
    state_network::SN
    "Neural network for computing outputs from fluxes"
    output_network::ON
    "Number of input variables"
    n_inputs::Int
    "Number of state variables"
    n_states::Int
    "Number of output variables"
    n_outputs::Int
    "HRU types (Nothing = 2D, Vector{Int} = 3D)"
    htypes::HT
    "Metadata about inputs, outputs, states"
    infos::I
end

"""
    NeuralBucket(; name, flux_network, state_network, output_network, ...)

Construct a NeuralBucket component.
"""
function NeuralBucket(;
    name::Symbol,
    flux_network::LuxCore.AbstractLuxLayer,
    state_network::LuxCore.AbstractLuxLayer,
    output_network::LuxCore.AbstractLuxLayer,
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

    flux_name = hasproperty(flux_network, :name) ? flux_network.name : :flux_net
    state_name = hasproperty(state_network, :name) ? state_network.name : :state_net
    output_name = hasproperty(output_network, :name) ? output_network.name : :output_net

    infos = HydroInfos(
        inputs=inputs,
        states=states,
        outputs=outputs,
        nns=[flux_name, state_name, output_name]
    )

    return NeuralBucket(
        name, flux_network, state_network, output_network,
        n_inputs, n_states, n_outputs,
        htypes, infos
    )
end

# ============================================================================
# Internal step function (not exposed as Lux interface)
# ============================================================================

"""
    _neural_bucket_step(bucket, x_t, hydro_state, nn_params, lux_states)

Single timestep computation for NeuralBucket. Internal function.

Returns `(output, new_hydro_state, new_lux_states)`.
"""
function _neural_bucket_step(bucket::NeuralBucket, x_t, hydro_state, nn_params, lux_states)
    # Step 1: Compute fluxes from states and inputs
    flux_input = vcat(hydro_state, x_t)
    fluxes, flux_st_new = bucket.flux_network(flux_input, nn_params.flux, lux_states.flux)

    # Step 2: Update states with residual connection (S_new = S_old + delta)
    state_input = vcat(hydro_state, fluxes)
    state_delta, state_st_new = bucket.state_network(state_input, nn_params.state, lux_states.state)
    new_hydro_state = hydro_state .+ state_delta

    # Step 3: Compute outputs from fluxes
    y_t, output_st_new = bucket.output_network(fluxes, nn_params.output, lux_states.output)

    new_lux_states = (
        flux = flux_st_new,
        state = state_st_new,
        output = output_st_new,
    )

    return y_t, new_hydro_state, new_lux_states
end

# ============================================================================
# NeuralBucket functor - HydroModel compatible interface
# ============================================================================

# 2D computation (single-node, htypes = Nothing)
function (bucket::NeuralBucket{FN,SN,ON,Nothing,I})(
    input::AbstractArray{T,2},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,2} where {FN,SN,ON,I,T}
    params = _as_componentvector(params)
    n_inputs, n_steps = size(input)
    @assert n_inputs == bucket.n_inputs "Input size mismatch: expected $(bucket.n_inputs), got $n_inputs"

    # Extract neural network parameters
    flux_name = bucket.infos.nns[1]
    state_name = bucket.infos.nns[2]
    output_name = bucket.infos.nns[3]

    nn_params = (
        flux = params[:nns][flux_name],
        state = params[:nns][state_name],
        output = params[:nns][output_name]
    )

    # Initialize Lux states
    rng = Random.default_rng()
    lux_states = (
        flux = LuxCore.initialstates(rng, bucket.flux_network),
        state = LuxCore.initialstates(rng, bucket.state_network),
        output = LuxCore.initialstates(rng, bucket.output_network),
    )

    # Initialize hydrological states
    initstates = get(kwargs, :initstates, zeros(T, bucket.n_states))
    hydro_state = T.(Vector(initstates))

    # Process sequence (RNN-style recurrence)
    all_states = Vector{AbstractVector{T}}(undef, n_steps)
    all_outputs = Vector{AbstractVector{T}}(undef, n_steps)

    for t in 1:n_steps
        x_t = input[:, t]
        y_t, hydro_state, lux_states = _neural_bucket_step(bucket, x_t, hydro_state, nn_params, lux_states)
        all_states[t] = hydro_state
        all_outputs[t] = y_t
    end

    # Stack and return [states; outputs] to match HydroBucket interface
    states_matrix = reduce(hcat, all_states)
    outputs_matrix = reduce(hcat, all_outputs)
    vcat(states_matrix, outputs_matrix)
end

# 3D computation (multi-node, htypes = Vector{Int})
function (bucket::NeuralBucket{FN,SN,ON,Vector{Int},I})(
    input::AbstractArray{T,3},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,3} where {FN,SN,ON,I,T}
    params = _as_componentvector(params)
    n_inputs, n_nodes, n_steps = size(input)
    @assert n_inputs == bucket.n_inputs "Input size mismatch: expected $(bucket.n_inputs), got $n_inputs"

    initstates_kw = get(kwargs, :initstates, nothing)

    # Process each node independently
    node_outputs = map(1:n_nodes) do node_idx
        node_initstates = if !isnothing(initstates_kw)
            # Extract per-node initial states
            state_vals = [initstates_kw[s][node_idx] for s in get_state_names(bucket)]
            reduce(vcat, state_vals)
        else
            zeros(T, bucket.n_states)
        end
        bucket(input[:, node_idx, :], params, config; initstates=node_initstates)
    end

    # Stack node outputs (n_states+n_outputs × n_nodes × n_timesteps)
    stack(node_outputs, dims=2)
end

# Error: single-node NeuralBucket receiving 3D input
function (bucket::NeuralBucket{FN,SN,ON,Nothing,I})(
    input::AbstractArray{T,3},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
) where {FN,SN,ON,I,T}
    error("NeuralBucket without htypes only accepts 2D input (variables × time).\n" *
          "For multi-node computation, provide htypes.\n" *
          "Got input shape: $(size(input))")
end

# Error: multi-node NeuralBucket receiving 2D input
function (bucket::NeuralBucket{FN,SN,ON,Vector{Int},I})(
    input::AbstractArray{T,2},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
) where {FN,SN,ON,I,T}
    error("NeuralBucket with htypes only accepts 3D input (variables × nodes × time).\n" *
          "For single-node computation, omit htypes.\n" *
          "Got input shape: $(size(input))")
end

# ============================================================================
# Utility Functions
# ============================================================================

"""
    create_neural_bucket(; name, n_inputs, n_states, n_outputs, ...)

Convenience function to create a NeuralBucket with default architecture.

# Architecture
- **Flux Network**: `(n_states + n_inputs) → hidden → n_fluxes`
- **State Network**: `(n_states + n_fluxes) → hidden → n_states`
- **Output Network**: `n_fluxes → n_outputs`
"""
function create_neural_bucket(;
    name::Symbol,
    n_inputs::Int,
    n_states::Int,
    n_outputs::Int,
    n_fluxes::Int=n_states,
    hidden_size::Int=16,
    inputs::Vector{Symbol},
    states::Vector{Symbol},
    outputs::Vector{Symbol},
    htypes::Optional{Vector{Int}}=nothing,
    flux_activation=tanh,
    state_activation=identity,
    output_activation=identity
)
    flux_network = Chain(
        Dense(n_states + n_inputs => hidden_size, flux_activation),
        Dense(hidden_size => n_fluxes),
        name = Symbol(name, :_flux)
    )
    state_network = Chain(
        Dense(n_states + n_fluxes => hidden_size, state_activation),
        Dense(hidden_size => n_states),
        name = Symbol(name, :_state)
    )
    output_network = Chain(
        Dense(n_fluxes => n_outputs, output_activation),
        name = Symbol(name, :_output)
    )

    return NeuralBucket(;
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
        htypes=htypes
    )
end

"""
    create_simple_neural_bucket(; name, n_inputs, n_states, n_outputs, ...)

Create a minimal NeuralBucket with direct linear transformations (no hidden layers).
"""
function create_simple_neural_bucket(;
    name::Symbol,
    n_inputs::Int,
    n_states::Int,
    n_outputs::Int,
    inputs::Vector{Symbol},
    states::Vector{Symbol},
    outputs::Vector{Symbol},
    htypes::Optional{Vector{Int}}=nothing,
)
    flux_network = Dense(n_states + n_inputs => n_states, name=Symbol(name, :_flux))
    state_network = Dense(n_states + n_states => n_states, name=Symbol(name, :_state))
    output_network = Dense(n_states => n_outputs, name=Symbol(name, :_output))

    return NeuralBucket(;
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
        htypes=htypes
    )
end

# Export interfaces
export NeuralFlux, NeuralBucket
export @neuralflux
export create_neural_bucket, create_simple_neural_bucket
