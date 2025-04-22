"""
    NeuralFlux{N} <: AbstractNeuralFlux

Represents a neural network-based flux component in a hydrological model. The type parameter 
`N` is used to encode the component name at the type level for better type stability.

# Arguments
- `inputs::Vector{Num}`: Vector of input variables
- `outputs::Vector{Num}`: Vector of output variables
- `chain::LuxCore.AbstractLuxLayer`: Neural network model
- `name::Union{Symbol,Nothing}=nothing`: Optional flux identifier
- `chain_name::Union{Symbol,Nothing}=nothing`: Optional neural network identifier

# Fields
- `chain::LuxCore.AbstractLuxLayer`: Neural network model
- `func::Function`: Generated function for flux calculations
- `infos::NamedTuple`: Component metadata

# Model Structure
The neural flux model consists of:
- Input variables: External variables fed to neural network
- Output variables: Neural network predictions
- Neural network: Machine learning model
- Generated function: Optimized implementation for calculations

## Metadata Components
The `infos` field tracks:
- `inputs`: Input variable names
- `outputs`: Output variable names
- `nns`: Neural network names
- `nn_inputs`: Neural network input names
- `nn_outputs`: Neural network output names
"""
struct NeuralFlux{N} <: AbstractNeuralFlux
    "chain of the neural network"
    chain::LuxCore.AbstractLuxLayer
    "Compiled function that calculates the flux using the neural network"
    func::Function
    "Information about the neural network's input and output structure"
    infos::NamedTuple

    function NeuralFlux(
        inputs::Vector{T},
        outputs::Vector{T},
        chain::LuxCore.AbstractLuxLayer;
        name::Union{Symbol,Nothing}=nothing,
        st=LuxCore.initialstates(StableRNG(42), chain),
        chain_name::Union{Symbol,Nothing}=nothing,
    ) where {T<:Num}
        #* Check chain name
        chain_name = chain_name === nothing ? chain.name : chain_name
        ps = LuxCore.initialparameters(StableRNG(42), chain)
        ps_axes = getaxes(ComponentVector(ps))
        nn_func = (x, p) -> LuxCore.apply(chain, x, ComponentVector(p, ps_axes), st)[1]
        nn_ps = @parameters $chain_name[1:length(ComponentVector(ps))]

        nn_input_name, nn_output_name = Symbol(chain_name, :_input), Symbol(chain_name, :_output)
        infos = (; inputs=inputs, outputs=outputs, nns=nn_ps, nn_inputs=nn_input_name, nn_outputs=nn_output_name)
        flux_name = isnothing(name) ? Symbol("##neural_flux#", hash(infos)) : name
        new{flux_name}(chain, nn_func, infos)
    end

    #* construct neural flux with input fluxes and output fluxes
    function NeuralFlux(fluxes::Pair{Vector{Num},Vector{Num}}, chain, name::Union{Symbol,Nothing}=nothing)
        return NeuralFlux(fluxes[1], fluxes[2], chain, name=name)
    end
end

(chain::LuxCore.AbstractLuxLayer)(inputs::AbstractVector{Num}) = (chain=chain, inputs=inputs, name=hasproperty(chain, :name) ? chain.name : nothing)

"""
    @neuralflux(eq::Expr)

Create a `NeuralFlux` using the syntax: `output ~ chain(inputs)`, where:
- `output` is the output variable or a vector of output variables
- `~` is the separator
- `chain` is a Lux neural network chain
- `inputs` is a vector of input variables

# Examples
```julia
@variables x, y, z
chain = Chain(Dense(2 => 10, relu), Dense(10 => 1), name=:my_net)

# Create a neural flux with a single output
flux1 = @neuralflux z ~ chain([x, y])

# Create a neural flux with multiple outputs
chain2 = Chain(Dense(2 => 16, relu), Dense(16 => 2), name=:multi_net)
flux2 = @neuralflux [z₁, z₂] ~ chain2([x, y])
```
"""
macro neuralflux(args...)
    name = length(args) == 1 ? nothing : args[1]
    eqs_expr = length(args) == 1 ? args[1] : args[2]
    @assert eqs_expr.head == :call && eqs_expr.args[1] == :~ "Expected equation in the form: outputs ~ chain(inputs)"
    lhs, rhs = eqs_expr.args[2], eqs_expr.args[3]  # Output variable(s) and Chain info expression
    return esc(quote
        outputs = $lhs isa Vector ? $lhs : [$lhs]
        chain_info = $rhs
        NeuralFlux(
            chain_info.inputs, outputs, chain_info.chain;
            name=$(name), chain_name=chain_info.name
        )
    end)
end


"""
    (flux::NeuralFlux)(input::AbstractArray, params::ComponentVector; kwargs...)

Apply a neural network-based flux model to input data for calculating water fluxes.

# Arguments
- `input::AbstractArray`: Input data array with the following possible dimensions:
  - `Matrix`: Time series data with shape (variables, timesteps)
  - `Array{3}`: Distributed data with shape (variables, nodes, timesteps)
- `params::ComponentVector`: Model parameters containing neural network weights and biases
- `kwargs`: Additional keyword arguments:
  - `ptyidx::AbstractVector{Int}`: Parameter type indices for distributed runs (default: all nodes)

# Returns
- `Matrix`: For 2D input, returns matrix of shape (output_variables, timesteps)
- `Array{3}`: For 3D input, returns array of shape (output_variables, nodes, timesteps)

# Description
This function applies neural network-based flux calculations to input time series data, 
supporting both single-node and distributed (multi-node) simulations. The neural network 
parameters are automatically retrieved from the appropriate component of the parameter vector.

## Input Data Structure
- Variables must be arranged along the first dimension
- For distributed runs, nodes are along the second dimension
- Time steps are always in the last dimension

## Parameter Handling
- Neural network parameters are accessed via `params[:nns][chain_name]`
- Parameters maintain their structure as defined in the neural network
- For distributed runs, the same network is applied to each node

## Example
```julia
# Define neural network flux
flux = NeuralFlux(
    [P, E] => [Q],           # input/output structure
    Dense(2, 1, tanh)        # neural network architecture
)

# Single-node simulation
output = flux(input, params)  # input: (2, timesteps)

# Multi-node simulation
output = flux(input, params)  # input: (2, n_nodes, timesteps)
```

# Notes
- Neural network parameters must be properly initialized before use
- The neural network architecture must match the input/output dimensions
- For distributed runs, the same network is shared across all nodes
"""
function (flux::NeuralFlux{N})(input::AbstractArray{T,2}, params::ComponentVector; kwargs...) where {T,N}
    nn_params = params[:nns][get_nn_names(flux)[1]]
    flux.func(input, nn_params)
end

function (flux::NeuralFlux{N})(input::AbstractArray{T,3}, params::ComponentVector; kwargs...) where {T,N}
    nn_params = params[:nns][get_nn_names(flux)[1]]
    #* array dims: (ts_len * node_names * var_names)
    stack(ntuple(i -> flux.func(input[:, i, :], nn_params), size(input)[2]), dims=2)
end

"""
    HydroLuxLayer(layer_func::Function)

Create a `HydroLuxLayer` from a function, where:
- `layer_func` is the function that calculates the output and state
"""
struct HydroLuxLayer <: Lux.AbstractRecurrentCell
    layer_func::Function
end

function (l::HydroLuxLayer)(input::AbstractArray, ps, st::NamedTuple)
    return LuxCore.apply(l, (input, eachslice(ps.initstates, dims=1)), ps, st)
end

function (l::HydroLuxLayer)(input::Tuple, ps, st::NamedTuple)
    forcing, states = input
    output, dstates = l.layer_func(forcing, ps, states)
    new_states = map(.+, states, dstates)
    return (stack(vcat(new_states, output), dims=1), new_states), st
end

"""    NeuralBucket{N} <: AbstractBucket

A recurrent neural network-based bucket for hydrological modeling that leverages Lux.jl's recurrence layers.

# Fields
- `fluxes::Vector{<:AbstractHydroFlux}`: List of process functions to be applied.
- `dfluxes::Vector{<:AbstractStateFlux}`: List of state derivative functions for state updates.
- `func::Function`: Precompiled recurrent neural network function for efficient evaluation.
- `infos::NamedTuple`: Metadata including input, output, state, parameter, and neural network names.

# Constructor
```julia
NeuralBucket(; name=nothing, fluxes, dfluxes)
```
- `name::Union{Symbol,Nothing}`: Optional identifier for the bucket. Defaults to an autogenerated name.
- `fluxes::Vector{<:AbstractHydroFlux}`: Required. Main computational elements.
- `dfluxes::Vector{<:AbstractStateFlux}`: Required. State derivatives for recurrent updates.

# Usage Example
```julia
bucket = NeuralBucket(
    name = :neural_bucket,
    fluxes = [HydroFlux(...), NeuralFlux(...)],
    dfluxes = [StateFlux(...)]
)
```

The NeuralBucket differs from HydroBucket by using Lux.jl's recurrence layers to handle state
updates over time, making it particularly suitable for sequence modeling tasks.
"""
struct NeuralBucket{N} <: AbstractBucket
    fluxes::Vector{<:AbstractHydroFlux}
    dfluxes::Vector{<:AbstractStateFlux}
    func::Function
    infos::NamedTuple

    function NeuralBucket(;
        fluxes::Vector{<:AbstractHydroFlux},
        dfluxes::Vector{<:AbstractStateFlux},
        name::Union{Symbol,Nothing}=nothing,
    )
        inputs, outputs, states = get_vars(fluxes, dfluxes)
        params, nns = reduce(union, get_params.(vcat(fluxes, dfluxes))), reduce(union, get_nns.(fluxes))
        infos = (; inputs=inputs, outputs=outputs, states=states, params=params, nns=nns)
        name = isnothing(name) ? Symbol("##neural_flux#", hash(infos)) : name
        rec_layer = HydroLuxLayer(build_nnlayer_func(fluxes, dfluxes, infos))
        rec_model = Recurrence(rec_layer; return_sequence=true, ordering=TimeLastIndex())
        ps, st = Lux.setup(Random.default_rng(), rec_model)
        rec_func(x, ps) = LuxCore.apply(rec_model, x, ps, st)[1]
        return new{name}(fluxes, dfluxes, rec_func, infos)
    end
end

"""
    @neuralbucket name begin
        fluxes = begin
            ...
        end
        dfluxes = begin
            ...
        end
    end

Creates a HydroBucket with the specified name, fluxes, and dfluxes.

# Arguments
- `name`: Symbol for the bucket name
- `fluxes`: Array of HydroFlux or NeuralFlux objects
- `dfluxes`: Array of StateFlux objects

# Example
```julia
@hydrobucket :bucket1 begin
    fluxes = [flux1, flux2]
    dfluxes = [dflux1, dflux2]
end
```
"""
macro neuralbucket(args...)
    name = length(args) == 1 ? nothing : args[1]
    expr = length(args) == 1 ? args[1] : args[2]

    @assert Meta.isexpr(expr, :block) "Expected a begin...end block after bucket name"
    fluxes_expr, dfluxes_expr = nothing, nothing
    for assign in filter(x -> !(x isa LineNumberNode), expr.args)
        @assert Meta.isexpr(assign, :(=)) "Expected assignments in the form 'fluxes = begin...end'"
        lhs, rhs = assign.args
        if lhs == :fluxes
            @assert Meta.isexpr(rhs, :block) "Expected 'fluxes' to be defined in a begin...end block"
            fluxes_expr = Expr(:vect, filter(x -> !(x isa LineNumberNode), rhs.args)...)
        elseif lhs == :dfluxes
            @assert Meta.isexpr(rhs, :block) "Expected 'dfluxes' to be defined in a begin...end block"
            dfluxes_expr = Expr(:vect, filter(x -> !(x isa LineNumberNode), rhs.args)...)
        else
            error("Unknown assignment: $lhs. Expected 'fluxes' or 'dfluxes'")
        end
    end
    @assert !isnothing(fluxes_expr) "'fluxes' must be specified"
    @assert !isnothing(dfluxes_expr) "'dfluxes' must be specified"

    return esc(:(NeuralBucket(name=$(name), fluxes=$fluxes_expr, dfluxes=$dfluxes_expr)))
end

"""    (bucket::NeuralBucket)(input::AbstractArray{T,3}, pas::ComponentVector; kwargs...) where {T,N}

Apply a neural bucket model to input data for calculating water fluxes and state updates over time.

# Arguments
- `input::AbstractArray{T,3}`: Input data array with shape (variables, nodes, timesteps).
- `pas::ComponentVector`: Model parameters containing both regular parameters and initial states.
- `kwargs`: Additional keyword arguments:
  - `ptyidx::AbstractVector{Int}`: Parameter type indices for distributed runs (default: all nodes).
  - `styidx::AbstractVector{Int}`: State type indices for distributed runs (default: all nodes).
  - `initstates::ComponentVector`: Initial state values (default: zeros for all states).

# Returns
- `Array{T,3}`: Output array with shape (output_variables, nodes, timesteps).

# Description
This function applies the neural bucket's recurrent model to input time series data,
supporting distributed (multi-node) simulations. The function handles state updates
internally through the recurrent neural network architecture.

The input array must have variables arranged along the first dimension, nodes along
the second dimension, and time steps along the third dimension.
"""
function (bucket::NeuralBucket{N})(input::AbstractArray{T,3}, pas::ComponentVector; kwargs...) where {T,N}
    ptyidx = get(kwargs, :ptyidx, 1:size(input, 2))
    styidx = get(kwargs, :styidx, 1:size(input, 2))
    state_names = get_state_names(bucket)
    initstates = get(kwargs, :initstates, ComponentVector(NamedTuple{Tuple(state_names)}(fill(zeros(T, size(input, 2)), length(state_names)))))
    new_pas = expand_component_params_and_initstates(pas, initstates, ptyidx, styidx)
    output = bucket.func(input, new_pas)
    return stack(output, dims=3)
end
