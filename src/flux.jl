"""
    HydroFlux{N} <: AbstractHydroFlux

Represents a simple flux component with mathematical formulas in a hydrological model. 
The type parameter `N` is used to encode the component name at the type level for better 
type stability.

# Arguments
- `inputs::Vector{Num}`: Vector of input variables
- `outputs::Vector{Num}`: Vector of output variables
- `params::Vector{Num}=Num[]`: Vector of parameter variables
- `exprs::Vector{Num}`: Vector of expressions for output calculations
- `name::Union{Symbol,Nothing}=nothing`: Optional flux identifier

# Fields
- `exprs::Vector{Num}`: Vector of mathematical expressions
- `func::Function`: Generated function for flux calculations
- `infos::NamedTuple`: Component metadata

# Description
HydroFlux is a type-stable implementation of flux calculations in a hydrological model. 
It automatically compiles mathematical expressions into efficient functions for computing 
flux values from inputs and parameters.

## Model Structure
The flux model consists of:
- Input variables: External variables used in calculations
- Output variables: Computed flux results
- Parameters: Model parameters used in calculations
- Expressions: Mathematical formulas defining the relationships

## Metadata Components
The `infos` field tracks:
- `inputs`: Input variable names
- `outputs`: Output variable names
- `params`: Parameter names

## Usage Notes
1. Construct with explicit variables:
   ```julia
   flux = HydroFlux(
       inputs=[x, y],     # input variables
       outputs=[z],       # output variables
       params=[k],        # parameters
       exprs=[k*x + y]    # expressions
   )
   ```

2. Construct with flux pairs:
   ```julia
   flux = HydroFlux(
       [x, y] => [z],    # input => output fluxes
       params=[k],        # parameters
       exprs=[k*x + y]    # expressions
   )
   ```

3. The number of expressions must match the number of output variables
4. Names are auto-generated from expression hashes if not provided
"""
struct HydroFlux{N} <: AbstractHydroFlux
    "Vector of expressions describing the formulas for output variables"
    exprs::Vector
    "Compiled function that calculates the flux"
    func::Function
    "Metadata about the flux, including input, output, and parameter names"
    infos::NamedTuple

    function HydroFlux(
        inputs::Vector{T},
        outputs::Vector{T},
        params::Vector{T};
        exprs::Vector{T},
        name::Union{Symbol,Nothing}=nothing,
    ) where {T}
        #* construct meta
        input_names = length(inputs) == 0 ? Symbol[] : tosymbol.(inputs)
        output_names = length(outputs) == 0 ? Symbol[] : tosymbol.(outputs)
        param_names = length(params) == 0 ? Symbol[] : tosymbol.(params)
        @assert length(exprs) == length(outputs) "The number of expressions and outputs must match, but got expressions: $(length(exprs)) and outputs: $(length(outputs))"
        #* build flux function
        flux_func = build_flux_func(inputs, outputs, params, exprs)
        #* use hash of exprs to name the flux
        infos = (; inputs=input_names, outputs=output_names, params=param_names)
        flux_name = isnothing(name) ? Symbol("##hydro_flux#", hash(infos)) : name
        return new{flux_name}(exprs, flux_func, infos)
    end

    #* construct hydro flux with input fluxes and output fluxes
    function HydroFlux(
        fluxes::Pair{Vector{Num},Vector{Num}},
        params::Vector{Num}=Num[];
        exprs::Vector,
        name::Union{Symbol,Nothing}=nothing,
    )
        return HydroFlux(fluxes[1], fluxes[2], params, exprs=exprs, name=name)
    end
end

"""
    (flux::AbstractHydroFlux)(input::AbstractArray, params::ComponentVector; kwargs...)

Apply a HydroFlux model to input data for calculating water fluxes.

# Arguments
- `input::AbstractArray`: Input data array with the following possible dimensions:
  - `Matrix`: Time series data with shape (variables, timesteps)
  - `Array{3}`: Distributed data with shape (variables, nodes, timesteps)
- `params::ComponentVector`: Model parameters organized in a component vector
- `kwargs`: Additional keyword arguments:
  - `ptyidx::AbstractVector{Int}`: Parameter type indices for distributed runs (default: all nodes)

# Returns
- `Matrix`: For 2D input, returns matrix of shape (output_variables, timesteps)
- `Array{3}`: For 3D input, returns array of shape (output_variables, nodes, timesteps)

# Description
This function applies the flux calculations to input time series data, supporting both 
single-node and distributed (multi-node) simulations. For distributed runs, parameters 
are automatically expanded to match the number of nodes.

## Input Data Structure
- Variables must be arranged along the first dimension
- For distributed runs, nodes are along the second dimension
- Time steps are always in the last dimension

## Parameter Handling
- Single-node: Parameters are used as is
- Multi-node: Parameters are expanded using `ptyidx` to match node count
- Component parameters maintain type stability

## Example
```julia
# Single-node simulation
flux = HydroFlux([P, E] => [Q], params=[k], exprs=[k*P - E])
output = flux(input, params)  # input: (2, timesteps)

# Multi-node simulation
output = flux(input, params, ptyidx=1:n_nodes)  # input: (2, n_nodes, timesteps)
```
"""
function (flux::HydroFlux{N})(input::AbstractArray{T,2}, params::ComponentVector; kwargs...) where {T,N}
    reduce(hcat, flux.func(eachslice(input, dims=1), params)) |> permutedims
end

function (flux::HydroFlux{N})(input::AbstractArray{T,3}, params::ComponentVector; kwargs...) where {T,N}
    ptyidx = get(kwargs, :ptyidx, 1:size(input, 2))
    expand_params = expand_component_params(params, ptyidx)
    output = flux.func(eachslice(input, dims=1), expand_params)
    return if length(output) == 1
        reshape(output[1], 1, size(input, 2), size(input, 3))
    else
        permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), output), (3, 1, 2))
    end
end

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

# Description
NeuralFlux is a type-stable implementation of neural network-based flux calculations in a 
hydrological model. It combines symbolic mathematics with neural networks to create flexible 
representations of complex hydrological processes.

## Model Structure
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

## Usage Notes
1. Construct with explicit variables:
   ```julia
   flux = NeuralFlux(
       inputs=[x, y],     # input variables
       outputs=[z],       # output variables
       chain=Dense(2, 1)  # neural network
   )
   ```

2. Construct with flux pairs:
   ```julia
   flux = NeuralFlux(
       [x, y] => [z],    # input => output fluxes
       chain=Dense(2, 1)  # neural network
   )
   ```

3. Neural network parameters are accessed via params[:nns][chain_name]
4. Names are auto-generated from metadata hashes if not provided
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
        chain_name::Union{Symbol,Nothing}=nothing,
    ) where {T<:Num}
        #* Check chain name
        chain_name = chain_name === nothing ? chain.name : chain_name
        input_names = length(inputs) == 0 ? Symbol[] : tosymbol.(inputs)
        output_names = length(outputs) == 0 ? Symbol[] : tosymbol.(outputs)

        ps, st = Lux.setup(StableRNG(42), chain)
        ps_axes = getaxes(ComponentVector(ps))
        nn_func = (x, p) -> LuxCore.apply(chain, x, ComponentVector(p, ps_axes), st)[1]

        nn_input_name, nn_output_name = Symbol(chain_name, :_input), Symbol(chain_name, :_output)
        infos = (; inputs=input_names, outputs=output_names, nns=[chain_name], nn_inputs=nn_input_name, nn_outputs=nn_output_name)
        flux_name = isnothing(name) ? Symbol("##neural_flux#", hash(infos)) : name
        new{flux_name}(chain, nn_func, infos)
    end

    #* construct neural flux with input fluxes and output fluxes
    function NeuralFlux(fluxes::Pair{Vector{Num},Vector{Num}}, chain, name::Union{Symbol,Nothing}=nothing)
        return NeuralFlux(fluxes[1], fluxes[2], chain, name=name)
    end
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
    output_arr = flux.func(input, nn_params)
    return output_arr
end

function (flux::NeuralFlux{N})(input::AbstractArray{T,3}, params::ComponentVector; kwargs...) where {T,N}
    nn_params = params[:nns][get_nn_names(flux)[1]]
    #* array dims: (ts_len * node_names * var_names)
    flux_output_vec = [flux.func(input[:, i, :], nn_params) for i in 1:size(input)[2]]
    flux_output_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), flux_output_vec)
    permutedims(flux_output_arr, (1, 3, 2))
end

"""
    StateFlux{N} <: AbstractStateFlux

Represents a state flux component in a hydrological model that describes how state variables 
change over time. The type parameter `N` is used to encode the component name at the type 
level for better type stability.

# Arguments
- `inputs::Vector{Num}`: Vector of input variables affecting state change
- `state::Num`: State variable being updated
- `params::Vector{Num}=Num[]`: Optional vector of parameter variables
- `expr::Num`: Expression defining the state derivative
- `name::Union{Symbol,Nothing}=nothing`: Optional flux identifier

# Fields
- `exprs::Vector{Num}`: Vector of state derivative expressions
- `infos::NamedTuple`: Component metadata including input, state, and parameter information

# Description
StateFlux is a type-stable implementation for representing state changes in hydrological models.
It defines how state variables evolve over time based on inputs, current state, and parameters.
The component can be constructed in three ways to accommodate different modeling needs:

## Model Structure
The state flux model consists of:
- Input variables: External variables affecting state change
- State variable: The variable being updated
- Parameters: Model parameters used in calculations
- Expression: Mathematical formula defining state derivative

## Metadata Components
The `infos` field tracks:
- `inputs`: Input variable names
- `states`: State variable names
- `params`: Parameter names
- `inflows`: Input flux names (when constructed from fluxes)
- `outflows`: Output flux names (when constructed from fluxes)

## Construction Methods
1. Explicit Definition:
   ```julia
   # Define state change as a function of inputs and parameters
   flux = StateFlux(
       inputs=[precip, evap],  # input variables
       state=soil_moisture,    # state variable
       params=[k₁, k₂],       # parameters
       expr=k₁*precip - k₂*evap  # state derivative
   )
   ```

2. Flux-Based Definition:
   ```julia
   # Define state change as difference between inflows and outflows
   flux = StateFlux(
       [inflow₁, inflow₂] => [outflow₁, outflow₂],  # input/output fluxes
       soil_moisture                                  # state variable
   )
   ```

3. State Transition:
   ```julia
   # Define direct state transition
   flux = StateFlux(
       old_state => new_state  # state transition pair
   )
   ```

## Usage Notes
1. State derivatives are automatically calculated from expressions
2. Names are auto-generated from state variables if not provided
3. StateFlux components must be used within a HydroBucket
4. Flux pairs automatically create mass-balance equations
"""
struct StateFlux{N} <: AbstractStateFlux
    "flux expressions to descripe the formula of the state variable"
    exprs::Vector{Num}
    "bucket information: keys contains: input, output, param, state"
    infos::NamedTuple

    function StateFlux(
        inputs::Vector{T},
        state::T,
        params::Vector{T}=T[];
        expr::T,
        name::Union{Symbol,Nothing}=nothing,
    ) where {T<:Num}
        #* Convert to a symbol based on the variable
        input_names = length(inputs) == 0 ? Symbol[] : tosymbol.(inputs)
        state_names = length(state) == 0 ? Symbol[] : [tosymbol(state)]
        param_names = length(params) == 0 ? Symbol[] : tosymbol.(params)
        infos = (; inputs=input_names, states=state_names, params=param_names)
        flux_name = isnothing(name) ? Symbol("##state_flux#", hash(infos)) : name
        return new{flux_name}([expr], infos)
    end
    #* construct state flux with input fluxes and output fluxes
    function StateFlux(fluxes::Pair{Vector{Num},Vector{Num}}, state::Num, name::Union{Symbol,Nothing}=nothing)
        expr = sum(fluxes[1]) - sum(fluxes[2])
        inflow_names, outflow_names, state_name = tosymbol.(fluxes[1]), tosymbol.(fluxes[2]), tosymbol(state)
        infos = (; inputs=vcat(inflow_names, outflow_names), states=[state_name], params=Symbol[], inflows=inflow_names, outflows=outflow_names)
        flux_name = isnothing(name) ? Symbol("##state_flux#", hash(infos)) : name
        return new{flux_name}([expr], infos)
    end
    #* construct state flux with state variables
    StateFlux(states::Pair{Num,Num}, name::Union{Symbol,Nothing}=nothing) = StateFlux([states[2]], states[1], expr=states[2] - states[1], name=name)
end

(::AbstractStateFlux)(::AbstractArray{T,2}, ::ComponentVector; kwargs...) where {T} = @error "State Flux cannot run directly, please using HydroFlux to run"
(::AbstractStateFlux)(::AbstractArray{T,3}, ::ComponentVector; kwargs...) where {T} = @error "State Flux cannot run directly, please using HydroFlux to run"