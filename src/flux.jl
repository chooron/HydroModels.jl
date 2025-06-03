"""
    HydroFlux{N} <: AbstractHydroFlux

Represents a simple flux component with mathematical formulas in a hydrological model. 

# Arguments
- `inputs::Vector{Num}`: Vector of input variables
- `outputs::Vector{Num}`: Vector of output variables
- `params::Vector{Num}=Num[]`: Vector of parameter variables
- `exprs::Vector{Num}`: Vector of expressions for output calculations
- `name::Union{Symbol,Nothing}=nothing`: Optional flux identifier

# Fields
- `exprs::Vector{Num}`: Vector of mathematical expressions
- `func::Function`: Generated function for flux calculations
- `infos::NamedTuple`: Metadata (inputs::Vector{Num}, outputs::Vector{Num}, params::Vector{Num})

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
"""
struct HydroFlux{N} <: AbstractHydroFlux
    "Vector of expressions describing the formulas for output variables"
    exprs::Vector
    "Compiled function that calculates the flux"
    func::Function
    "Metadata about the flux, including input, output, and parameter names"
    infos::NamedTuple

    function HydroFlux(
        inputs::Vector{T}, outputs::Vector{T}, params::Vector{T};
        exprs::Vector{T}, name::Union{Symbol,Nothing}=nothing,
    ) where {T}
        @assert length(exprs) == length(outputs) "The number of expressions and outputs must match, but got expressions: $(length(exprs)) and outputs: $(length(outputs))"
        flux_func = build_flux_func(inputs, outputs, params, exprs)
        infos = (; inputs=inputs, outputs=outputs, params=params)
        flux_name = isnothing(name) ? Symbol("##hydro_flux#", hash(infos)) : name
        return new{flux_name}(exprs, flux_func, infos)
    end

    function HydroFlux(
        fluxes::Pair{Vector{Num},Vector{Num}}, params::Vector{Num}=Num[];
        exprs::Vector, name::Union{Symbol,Nothing}=nothing,
    )
        return HydroFlux(fluxes[1], fluxes[2], params, exprs=exprs, name=name)
    end
end

"""
    @hydroflux(name, eqs...)

Create a `HydroFlux` from a set of equations, where:
- Left sides of equations are output variables
- Right sides contain input variables and parameters
- Parameters are identified using ModelingToolkit's `isparameter`
- `name` is an optional name for the flux (provide as first argument)

# Examples
```julia
@variables x, y, z
@parameters a, b

# Create a named flux with one equation
flux1 = @hydroflux :my_flux z = a*x + b*y

# Create a flux with multiple equations
flux2 = @hydroflux begin
    z₁ = a*x + b*y
    z₂ = x^2 + y^2
end
```
"""
macro hydroflux(args...)
    name = length(args) == 1 ? nothing : args[1]
    eqs_expr = length(args) == 1 ? args[1] : args[2]
    vect_eqs_expr = if Meta.isexpr(eqs_expr, :block)
        Expr(:vect, filter(arg -> !(arg isa LineNumberNode) && !Meta.isexpr(arg, :line), eqs_expr.args)...)
    elseif Meta.isexpr(eqs_expr, :for)
        loop_var = eqs_expr.args[1].args[1]
        range_expr = eqs_expr.args[1].args[2]
        loop_body = eqs_expr.args[2]
        range_val = if Meta.isexpr(range_expr, :call) && range_expr.args[1] == :(:)
            if length(range_expr.args) == 3
                range_expr.args[2]:range_expr.args[3]
            elseif length(range_expr.args) == 4
                range_expr.args[2]:range_expr.args[3]:range_expr.args[4]
            end
        else
            eval(range_expr)
        end
        all_equations = []
        if !Meta.isexpr(loop_body, :block)
            loop_body = Expr(:block, loop_body)
        end
        equations = filter(x -> !(x isa LineNumberNode) && !Meta.isexpr(x, :line), loop_body.args)
        for i_val in range_val
            for eq in equations
                new_eq = deepcopy(eq)
                replace_loop_var!(new_eq, loop_var, i_val)
                push!(all_equations, new_eq)
            end
        end
        Expr(:vect, all_equations...)
    else
        Expr(:vect, eqs_expr)
    end

    return esc(quote
        Num = HydroModels.Num
        equations = $vect_eqs_expr
        lhs_terms = Num.([eq.lhs for eq in equations])
        rhs_terms = Num.([eq.rhs for eq in equations])

        all_vars = Num.(mapreduce(HydroModels.get_variables, union, rhs_terms, init=Set{Num}()))
        inputs = Num.(filter(x -> !HydroModels.isparameter(x), collect(all_vars)))
        params = Num.(filter(x -> HydroModels.isparameter(x), collect(all_vars)))
        inputs = setdiff(inputs, lhs_terms)

        HydroFlux(inputs, lhs_terms, params, exprs=Num.(rhs_terms), name=$(name))
    end)
end

function replace_loop_var!(expr, var_name, value)
    if expr isa Expr
        for i in 1:length(expr.args)
            if expr.args[i] == var_name
                expr.args[i] = value
            else
                replace_loop_var!(expr.args[i], var_name, value)
            end
        end
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
    ptyidx = get(kwargs, :ptyidx, collect(1:size(input, 2)))
    expand_params = expand_component_params(params, get_param_names(flux), ptyidx)
    output = flux.func(eachslice(input, dims=1), expand_params)
    stack(output, dims=1)
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
        infos = (; inputs=inputs, states=[state], params=params)
        flux_name = isnothing(name) ? Symbol("##state_flux#", hash(infos)) : name
        return new{flux_name}([expr], infos)
    end
    #* construct state flux with input fluxes and output fluxes
    function StateFlux(fluxes::Pair{Vector{Num},Vector{Num}}, state::Num, name::Union{Symbol,Nothing}=nothing)
        expr = sum(fluxes[1]) - sum(fluxes[2])
        infos = (; inputs=vcat(fluxes[1], fluxes[2]), states=[state], params=Num[])
        flux_name = isnothing(name) ? Symbol("##state_flux#", hash(infos)) : name
        return new{flux_name}([expr], infos)
    end
    #* construct state flux with state variables
    StateFlux(states::Pair{Num,Num}, name::Union{Symbol,Nothing}=nothing) = StateFlux([states[2]], states[1], expr=states[2] - states[1], name=name)
end

"""
    stateflux(name, eq)

Create a `StateFlux` from an equation, where:
- Left side of the equation is the state variable
- Right side contains the expression for the state variable's change
- `name` is an optional name for the flux (provide as first argument)

# Examples
```julia
@variables x, y, z
@parameters a, b

# Create a named state flux
flux1 = @stateflux_build :my_flux z = a*x + b*y

# Create a state flux without a name
flux2 = @stateflux_build z = a*x + b*y
```
"""
macro stateflux(args...)
    # Process arguments
    name = length(args) == 1 ? nothing : args[1]
    eq_expr = length(args) == 1 ? args[1] : args[2]

    # Handle both = and ~ operators
    if Meta.isexpr(eq_expr, :(=))
        lhs, rhs = eq_expr.args[1], eq_expr.args[2]
    elseif Meta.isexpr(eq_expr, :call) && eq_expr.args[1] == :~
        lhs, rhs = eq_expr.args[2], eq_expr.args[3]
    else
        error("Expected equation (using = or ~), got: $eq_expr")
    end

    return esc(quote
        Num = HydroModels.Num
        all_vars = Num.(HydroModels.get_variables($rhs))
        inputs = Num.(filter(x -> !HydroModels.isparameter(x), collect(all_vars)))
        params = Num.(filter(x -> HydroModels.isparameter(x), collect(all_vars)))
        inputs = setdiff(inputs, only($lhs))

        StateFlux(inputs, only($lhs), params, expr=Num($rhs), name=$(name))
    end)
end

(::StateFlux)(::AbstractArray, ::ComponentVector; kwargs...) = @error "State Flux cannot run directly, please using HydroFlux to run"

