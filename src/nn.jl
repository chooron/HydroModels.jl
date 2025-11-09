"""
Neural network components module - integrates Lux.jl neural networks with HydroModels framework.

This module provides:
- `NeuralFlux`: Neural network-driven flux component
- `NeuralBucket`: Lux-based hydrological bucket with RNN-like recurrence structure

# Design Philosophy

The NeuralBucket leverages the structural similarity between hydrological ODEs and RNNs:
- **Hydrological ODE**: `dS/dt = f(S, input, params)` → iterate over time
- **RNN**: `h_t = f(h_{t-1}, x_t, params)` → iterate over time

Both involve:
1. Maintain states across time steps
2. Compute outputs from states and inputs
3. Update states based on computations

By implementing HydroBucket logic using Lux's recurrence mechanism, we gain:
- Native Lux parameter management
- Seamless automatic differentiation
- Composability with other Lux layers
- Efficient GPU acceleration
"""

# ============================================================================
# NeuralFlux - Neural Network Flux Component
# ============================================================================

"""
    NeuralFlux{C, CF, NF, NT} <: AbstractNeuralFlux

Represents a flux component driven by a neural network.

It wraps a Lux.AbstractLuxLayer and connects it to symbolic variables for integration 
into a hydrological model.

$(FIELDS)

# Examples
```jldoctest
julia> @variables x, y, z
julia> using Lux
julia> chain = Chain(Dense(2 => 10, relu), Dense(10 => 1), name=:my_net)
julia> flux = NeuralFlux([x, y], [z], chain)
```
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
end

"""
    @neuralflux [name] eq

Macro to conveniently create a NeuralFlux from an equation.

# Usage
The macro takes an optional name and an equation of the form `output ~ chain(inputs)`.

# Examples
```jldoctest
julia> @variables x, y, z, z₁, z₂
julia> chain = Chain(Dense(2 => 10, relu), Dense(10 => 1), name=:my_net)
julia> # Single output
julia> flux1 = @neuralflux z ~ chain([x, y])
julia> # With an optional name
julia> flux2 = @neuralflux :my_flux z ~ chain([x, y])
julia> # Multiple outputs
julia> chain2 = Chain(Dense(2 => 16, relu), Dense(16 => 2), name=:multi_net)
julia> flux3 = @neuralflux [z₁, z₂] ~ chain2([x, y])
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

# NeuralFlux computation for 2D input (single node)
function (flux::NeuralFlux)(
    input::AbstractArray{T,2},
    params::ComponentVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,2} where {T}
    nn_params = params[:nns][get_nn_names(flux)[1]]
    flux.chain_func(flux.norm_func(input), nn_params)
end

# NeuralFlux computation for 3D input (multi-node)
function (flux::NeuralFlux)(
    input::AbstractArray{T,3},
    params::ComponentVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,3} where {T}
    nn_params = params[:nns][get_nn_names(flux)[1]]
    norm_input = flux.norm_func(input)
    # Use map instead of ntuple for better performance
    node_outputs = map(1:size(input, 2)) do i
        flux.chain_func(norm_input[:, i, :], nn_params)
    end
    stack(node_outputs, dims=2)
end

# ============================================================================
# NeuralBucket - Lux-based Hydrological Bucket with RNN-like Structure
# ============================================================================

"""
    NeuralBucket{FN, SN, ON, HT, I} <: LuxCore.AbstractLuxLayer

A Lux-based hydrological bucket that mimics the ODE structure using RNN-like recurrence.

# Design Concept

Traditional HydroBucket uses ODEs:
```
dS/dt = f(S, input, params)  →  S_t = S_{t-1} + dt * f(S_{t-1}, input_t, params)
output_t = g(S_t, input_t, params)
```

NeuralBucket reformulates this as RNN-like recurrence:
```
flux_t = flux_network(S_{t-1}, input_t, params)
S_t = state_update_network(S_{t-1}, flux_t, params)
output_t = output_network(flux_t, S_t, params)
```

$(FIELDS)

# Architecture

The bucket consists of three neural networks:
1. **Flux Network**: Computes fluxes from states and inputs
2. **State Update Network**: Updates states based on fluxes (can be simple addition)
3. **Output Network**: Extracts outputs from fluxes and states

# Examples
```julia
using Lux, HydroModels

# Define networks for flux computation
flux_net = Chain(
    Dense(4 => 16, tanh),   # 2 inputs + 2 states
    Dense(16 => 3),          # 3 fluxes
    name = :flux_net
)

# State update can be simple or complex
state_net = Chain(
    Dense(5 => 8, relu),     # 2 states + 3 fluxes
    Dense(8 => 2),           # 2 new states
    name = :state_net
)

# Output network
output_net = Chain(
    Dense(3 => 1),           # 3 fluxes → 1 output
    name = :output_net
)

neural_bucket = NeuralBucket(
    name = :my_bucket,
    flux_network = flux_net,
    state_network = state_net,
    output_network = output_net,
    n_inputs = 2,
    n_states = 2,
    n_outputs = 1,
    inputs = [:prcp, :temp],
    states = [:snowpack, :soilwater],
    outputs = [:runoff]
)
```

# Notes
- Mimics HydroBucket structure but uses neural networks
- No explicit ODE solver needed (discrete time stepping)
- States are hydrological states (not RNN hidden states)
- Compatible with Lux ecosystem
- Supports multi-node computation
"""
struct NeuralBucket{FN,SN,ON,HT,I} <: LuxCore.AbstractLuxLayer
    "Bucket name"
    name::Symbol
    "Neural network for computing fluxes from states and inputs"
    flux_network::FN
    "Neural network for updating states from previous states and fluxes"
    state_network::SN
    "Neural network for computing outputs from fluxes/states"
    output_network::ON
    "Number of input variables"
    n_inputs::Int
    "Number of state variables"
    n_states::Int
    "Number of output variables"
    n_outputs::Int
    "HRU types for multi-node scenarios"
    hru_types::HT
    "Metadata about inputs, outputs, states"
    infos::I
end

"""
    NeuralBucket(;
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
        hru_types::Vector{Int}=Int[]
    )

Construct a NeuralBucket component.

# Arguments
- `name`: Name identifier for the bucket
- `flux_network`: Neural network for flux computation (input: states+inputs, output: fluxes)
- `state_network`: Neural network for state updates (input: states+fluxes, output: new states)
- `output_network`: Neural network for output computation (input: fluxes/states, output: outputs)
- `n_inputs`: Number of input variables
- `n_states`: Number of state variables
- `n_outputs`: Number of output variables
- `inputs`: Names of input variables
- `states`: Names of state variables
- `outputs`: Names of output variables
- `hru_types`: HRU type indices for parameter sharing in multi-node scenarios

# Returns
A NeuralBucket instance that can be used in HydroModel
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
    hru_types::Vector{Int}=Int[]
)
    @assert length(inputs) == n_inputs "Number of input names must match n_inputs"
    @assert length(states) == n_states "Number of state names must match n_states"
    @assert length(outputs) == n_outputs "Number of output names must match n_outputs"
    
    # Get neural network names from layers
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
        hru_types, infos
    )
end

# Implement Lux interface for NeuralBucket
function LuxCore.initialparameters(rng::AbstractRNG, bucket::NeuralBucket)
    flux_params = LuxCore.initialparameters(rng, bucket.flux_network)
    state_params = LuxCore.initialparameters(rng, bucket.state_network)
    output_params = LuxCore.initialparameters(rng, bucket.output_network)
    
    return (
        flux = flux_params,
        state = state_params,
        output = output_params
    )
end

function LuxCore.initialstates(rng::AbstractRNG, bucket::NeuralBucket)
    flux_states = LuxCore.initialstates(rng, bucket.flux_network)
    state_states = LuxCore.initialstates(rng, bucket.state_network)
    output_states = LuxCore.initialstates(rng, bucket.output_network)
    
    return (
        flux = flux_states,
        state = state_states,
        output = output_states,
        hydro_state = nothing  # Will be initialized on first call
    )
end

"""
    (bucket::NeuralBucket)(x_t, ps, st)

Standard Lux layer call interface for NeuralBucket (single time step).

Mimics the hydrological ODE structure:
1. Compute fluxes from current states and inputs
2. Update states based on fluxes
3. Compute outputs from fluxes/states

# Arguments
- `x_t`: Input at time t, shape (n_inputs,) or (n_inputs, batch_size)
- `ps`: Parameters (from LuxCore.initialparameters)
- `st`: States (from LuxCore.initialstates), includes `hydro_state`

# Returns
- `y_t`: Output at time t, shape (n_outputs,) or (n_outputs, batch_size)
- `st_new`: Updated states including new `hydro_state`
"""
function (bucket::NeuralBucket)(x_t::AbstractArray, ps::NamedTuple, st::NamedTuple)
    # Initialize hydrological state if needed
    hydro_state = st.hydro_state
    if isnothing(hydro_state)
        # Initialize states with zeros
        batch_size = ndims(x_t) == 1 ? 1 : size(x_t, 2)
        if batch_size == 1
            hydro_state = zeros(eltype(x_t), bucket.n_states)
        else
            hydro_state = zeros(eltype(x_t), bucket.n_states, batch_size)
        end
    end
    
    # Step 1: Compute fluxes from states and inputs
    # Concatenate states and inputs
    if ndims(x_t) == 1
        flux_input = vcat(hydro_state, x_t)
    else
        flux_input = vcat(hydro_state, x_t)
    end
    
    fluxes, flux_st_new = bucket.flux_network(flux_input, ps.flux, st.flux)
    
    # Step 2: Update states based on previous states and fluxes
    # Concatenate states and fluxes
    if ndims(x_t) == 1
        state_input = vcat(hydro_state, fluxes)
    else
        state_input = vcat(hydro_state, fluxes)
    end
    
    new_hydro_state, state_st_new = bucket.state_network(state_input, ps.state, st.state)
    
    # Step 3: Compute outputs from fluxes (or fluxes + states)
    # Here we use fluxes as input, but could also use [fluxes; new_hydro_state]
    y_t, output_st_new = bucket.output_network(fluxes, ps.output, st.output)
    
    # Package updated states
    st_new = (
        flux = flux_st_new,
        state = state_st_new,
        output = output_st_new,
        hydro_state = new_hydro_state
    )
    
    return y_t, st_new
end

"""
    (bucket::NeuralBucket)(
        input::AbstractArray{T,2},
        params::ComponentVector,
        config::ConfigType=default_config();
        initstates::ComponentVector=ComponentVector(),
        kwargs...
    ) where {T}

HydroModels interface for single-node NeuralBucket computation.

Processes a sequence of inputs (variables × time) and returns outputs for all time steps.
This mimics the HydroBucket interface but uses neural networks for computation.

# Arguments
- `input`: Input matrix (n_inputs × n_timesteps)
- `params`: ComponentVector containing neural network parameters under `:nns` key
- `config`: HydroConfig (not used, for interface compatibility)
- `initstates`: Initial hydrological states (optional)

# Returns
Output matrix (n_outputs × n_timesteps)

# Example
```julia
# Create bucket
bucket = NeuralBucket(...)

# Prepare parameters
ps_lux = LuxCore.initialparameters(rng, bucket)
params = ComponentVector(nns = (flux_net = ps_lux.flux, ...))

# Run
input = rand(Float32, 2, 100)  # 2 inputs, 100 timesteps
output = bucket(input, params, HydroConfig())
```
"""
function (bucket::NeuralBucket)(
    input::AbstractArray{T,2},
    params::ComponentVector,
    config::ConfigType=default_config();
    initstates::ComponentVector=ComponentVector(),
    kwargs...
) where {T}
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
    st = LuxCore.initialstates(rng, bucket)
    
    # Initialize hydrological states from initstates if provided
    if !isempty(initstates) && haskey(initstates, :states)
        # Convert initstates to array format
        state_vals = [initstates[:states][s] for s in get_state_names(bucket)]
        st = merge(st, (hydro_state = reduce(vcat, state_vals),))
    end
    
    # Process sequence (RNN-style recurrence)
    outputs = Vector{AbstractVector{T}}(undef, n_steps)
    
    for t in 1:n_steps
        x_t = input[:, t]
        y_t, st = bucket(x_t, nn_params, st)
        outputs[t] = y_t
    end
    
    # Stack outputs to matrix
    return reduce(hcat, outputs)
end

"""
    (bucket::NeuralBucket)(
        input::AbstractArray{T,3},
        params::ComponentVector,
        config::ConfigType=default_config();
        initstates::ComponentVector=ComponentVector(),
        kwargs...
    ) where {T}

HydroModels interface for multi-node NeuralBucket computation.

Processes sequences for multiple nodes in parallel (variables × nodes × time).
Each node maintains its own hydrological state, similar to traditional HydroBucket.

# Arguments
- `input`: Input array (n_inputs × n_nodes × n_timesteps)
- `params`: ComponentVector with parameters (supports hru_types for parameter sharing)
- `config`: HydroConfig (not used, for interface compatibility)
- `initstates`: Initial hydrological states (optional)

# Returns
Output array (n_outputs × n_nodes × n_timesteps)

# Example
```julia
# Multi-node input
input = rand(Float32, 2, 5, 100)  # 2 inputs, 5 nodes, 100 timesteps

# Run bucket (each node has independent states)
output = bucket(input, params, config)  # (n_outputs × 5 × 100)
```
"""
function (bucket::NeuralBucket)(
    input::AbstractArray{T,3},
    params::ComponentVector,
    config::ConfigType=default_config();
    initstates::ComponentVector=ComponentVector(),
    kwargs...
) where {T}
    n_inputs, n_nodes, n_steps = size(input)
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
    
    # Initialize states
    rng = Random.default_rng()
    
    # Process each node separately (each maintains its own hydrological state)
    node_outputs = map(1:n_nodes) do node_idx
        st = LuxCore.initialstates(rng, bucket)
        
        # Initialize hydrological states from initstates if provided
        if !isempty(initstates) && haskey(initstates, :states)
            state_vals = [initstates[:states][s][node_idx] for s in get_state_names(bucket)]
            st = merge(st, (hydro_state = reduce(vcat, state_vals),))
        end
        
        outputs_node = Vector{AbstractVector{T}}(undef, n_steps)
        
        for t in 1:n_steps
            x_t = input[:, node_idx, t]
            y_t, st = bucket(x_t, nn_params, st)
            outputs_node[t] = y_t
        end
        
        reduce(hcat, outputs_node)
    end
    
    # Stack node outputs (n_outputs × n_nodes × n_timesteps)
    return stack(node_outputs, dims=2)
end

# ============================================================================
# Utility Functions
# ============================================================================

"""
    create_neural_bucket(;
        name::Symbol,
        n_inputs::Int,
        n_states::Int,
        n_outputs::Int,
        n_fluxes::Int=n_states,
        hidden_size::Int=16,
        inputs::Vector{Symbol},
        states::Vector{Symbol},
        outputs::Vector{Symbol},
        hru_types::Vector{Int}=Int[],
        flux_activation=tanh,
        state_activation=identity,
        output_activation=identity
    )

Convenience function to create a NeuralBucket with default architecture.

# Architecture
- **Flux Network**: `(n_states + n_inputs) → hidden → n_fluxes`
- **State Network**: `(n_states + n_fluxes) → hidden → n_states`
- **Output Network**: `n_fluxes → n_outputs`

# Example
```julia
bucket = create_neural_bucket(
    name = :rainfall_runoff,
    n_inputs = 2,           # prcp, temp
    n_states = 2,           # snowpack, soilwater
    n_outputs = 1,          # runoff
    n_fluxes = 3,           # snowfall, melt, runoff
    hidden_size = 16,
    inputs = [:prcp, :temp],
    states = [:snowpack, :soilwater],
    outputs = [:runoff]
)
```
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
    hru_types::Vector{Int}=Int[],
    flux_activation=tanh,
    state_activation=identity,
    output_activation=identity
)
    @assert length(inputs) == n_inputs "Number of input names must match n_inputs"
    @assert length(states) == n_states "Number of state names must match n_states"
    @assert length(outputs) == n_outputs "Number of output names must match n_outputs"
    
    # Flux network: (states + inputs) → fluxes
    flux_network = Chain(
        Dense(n_states + n_inputs => hidden_size, flux_activation),
        Dense(hidden_size => n_fluxes),
        name = Symbol(name, :_flux)
    )
    
    # State network: (states + fluxes) → new_states
    state_network = Chain(
        Dense(n_states + n_fluxes => hidden_size, state_activation),
        Dense(hidden_size => n_states),
        name = Symbol(name, :_state)
    )
    
    # Output network: fluxes → outputs
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
        hru_types=hru_types
    )
end

"""
    create_simple_neural_bucket(;
        name::Symbol,
        n_inputs::Int,
        n_states::Int,
        n_outputs::Int,
        inputs::Vector{Symbol},
        states::Vector{Symbol},
        outputs::Vector{Symbol},
        hru_types::Vector{Int}=Int[]
    )

Create a minimal NeuralBucket with direct linear transformations (no hidden layers).

This is useful for debugging or as a baseline model.

# Example
```julia
bucket = create_simple_neural_bucket(
    name = :simple,
    n_inputs = 2,
    n_states = 1,
    n_outputs = 1,
    inputs = [:prcp, :temp],
    states = [:storage],
    outputs = [:flow]
)
```
"""
function create_simple_neural_bucket(;
    name::Symbol,
    n_inputs::Int,
    n_states::Int,
    n_outputs::Int,
    inputs::Vector{Symbol},
    states::Vector{Symbol},
    outputs::Vector{Symbol},
    hru_types::Vector{Int}=Int[]
)
    # Simple linear transformations
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
        hru_types=hru_types
    )
end

# Export interfaces
export NeuralFlux, NeuralBucket
export @neuralflux
export create_neural_bucket, create_simple_neural_bucket

