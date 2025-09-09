"""
    HydroFlux{MS, E, F, HT, I} <: AbstractHydroFlux

Represents a simple flux component defined by mathematical expressions.

It automatically parses expressions to determine inputs, outputs, and parameters, and compiles them into a callable function for efficient computation.

$(FIELDS)
"""
struct HydroFlux{MS,E,F,HT,I} <: AbstractHydroFlux
    "flux name"
    name::Symbol
    "Vector of expressions describing the formulas for output variables"
    exprs::E
    "Compiled function that calculates the flux"
    func::F
    "nodes type"
    hru_types::HT
    "Metadata about the flux, including input, output, and parameter names"
    infos::I

    function HydroFlux(;
        exprs::E,
        hru_types::Vector{Int}=Int[],
        name::Optional{Symbol}=nothing,
    ) where {E}
        # parse expressions and extract variables
        outputs = Num.([eq.lhs for eq in exprs])
        eqs = Num.([eq.rhs for eq in exprs])
        all_vars = Num.(mapreduce(get_variables, union, eqs, init=Set{Num}()))
        inputs = setdiff(Num.(filter(x -> !isparameter(x), collect(all_vars))), outputs)
        params = Num.(filter(x -> isparameter(x), collect(all_vars)))
        assert_msg = "The number of expressions and outputs must match" *
                     " but got expressions: $(length(exprs)) and outputs: $(length(outputs))"
        @assert length(exprs) == length(outputs) assert_msg

        # build the function for flux calculations
        infos = HydroInfos(
            inputs=!isempty(inputs) ? tosymbol.(inputs) : Symbol[],
            params=!isempty(params) ? tosymbol.(params) : Symbol[],
            outputs=!isempty(outputs) ? tosymbol.(outputs) : Symbol[]
        )
        flux_func = build_flux_func(eqs, infos)
        flux_name = isnothing(name) ? Symbol("##hydro_flux#", hash(infos)) : name

        return new{length(hru_types) > 1,typeof(eqs),typeof(flux_func),typeof(hru_types),typeof(infos)}(
            flux_name, eqs, flux_func, hru_types, infos
        )
    end
end

"""
    @hydroflux(name, eqs...)

A macro to conveniently create a `HydroFlux` component from a set of equations.

It parses the given equations to identify output variables (left-hand sides), input variables, and parameters (right-hand sides). Parameters are distinguished from variables using `ModelingToolkit's` `isparameter` function.

# Arguments
- `name`: An optional `Symbol` to name the flux component. If it's the only argument, it's treated as the expression.
- `eqs...`: One or more equations defining the flux. These can be provided as a single equation, multiple equations, or within a `begin...end` block.
"""
macro hydroflux(args...)
    name = length(args) == 1 ? nothing : args[1]
    eqs_expr = length(args) == 1 ? args[1] : args[2]
    vect_eqs_expr = if Meta.isexpr(eqs_expr, :block)
        Expr(:tuple, filter(arg -> !(arg isa LineNumberNode) && !Meta.isexpr(arg, :line) && !Meta.isexpr(arg, :(=)), eqs_expr.args)...)
    else
        Expr(:tuple, eqs_expr)
    end

    hru_types_val = :(Int[])
    if Meta.isexpr(eqs_expr, :block)
        for assign in filter(x -> Meta.isexpr(x, :(=)), eqs_expr.args)
            # @assert Meta.isexpr(assign, :(=)) "Expected assignments in the form 'fluxes = begin...end'"
            lhs, rhs = assign.args
            if lhs == :hru_types
                hru_types_val = rhs
            end
        end
    end

    for var_name in extract_variables(vect_eqs_expr)
        if !@isdefined(var_name)
            expr_str = string(vect_eqs_expr)
            return :(error("Undefined variable '", $(string(var_name)), "' detected in expression: `", $expr_str, "`"))
        end
    end
    return esc(:(HydroFlux(exprs=$vect_eqs_expr, name=$name, hru_types=$hru_types_val)))
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
$(TYPEDSIGNATURES)

Generates a runtime function for flux calculations based on symbolic expressions.
The function is specialized for the given expressions and variable names in `infos`.
"""
function build_flux_func(exprs::Vector{Num}, infos::HydroModelCore.HydroInfos)
    flux_exprs = map(expr -> :(@. $(simplify_expr(toexpr(expr)))), exprs)
    input_assign_calls = generate_var_assignments(vars=infos.inputs, target=:inputs)
    params_assign_calls = generate_param_assignments(params=infos.params)
    compute_calls = [:($o = $expr) for (o, expr) in zip(infos.outputs, flux_exprs)]

    flux_func_expr = :(function (inputs, pas)
        Base.@_inline_meta
        $(input_assign_calls...)
        $(params_assign_calls...)
        $(compute_calls...)
        return [$((infos.outputs)...)]
    end)
    return @RuntimeGeneratedFunction(flux_func_expr)
end



function (flux::HydroFlux)(input::AbstractArray{T,2}, params::ComponentVector; kwargs...) where {T}
    stack(flux.func(eachslice(input, dims=1), params), dims=1)
end

function (flux::HydroFlux{true})(input::AbstractArray{T,3}, params::ComponentVector; kwargs...) where {T}
    expand_params = expand_component_params(params, get_param_names(flux), flux.hru_types)
    output = flux.func(eachslice(input, dims=1), expand_params)
    stack(output, dims=1)
end

"""
    StateFlux{N, E, I} <: AbstractStateFlux

Represents a state flux, defining the rate of change (derivative) for a state variable.

This is a declarative component used within a `HydroBucket` or `HydroRoute` to define the system's differential equations. It does not perform calculations itself.

$(FIELDS)
"""
struct StateFlux{N,E,I} <: AbstractStateFlux
    "Vector of expressions defining the state's rate of change."
    exprs::E
    "Metadata about the component, including input, state, and parameter names."
    infos::I

    function StateFlux(;
        exprs::E, name::Optional{Symbol}=nothing
    ) where {E}
        states = Num.([eq.lhs for eq in exprs])
        eqs = Num.([eq.rhs for eq in exprs])
        all_vars = Num.(mapreduce(HydroModels.get_variables, union, eqs, init=Set{Num}()))
        inputs = setdiff(Num.(filter(x -> !HydroModels.isparameter(x), collect(all_vars))), states)
        params = Num.(filter(x -> HydroModels.isparameter(x), collect(all_vars)))
        infos = HydroInfos(
            inputs=!isempty(inputs) ? tosymbol.(inputs) : Symbol[],
            states=!isempty(states) ? tosymbol.(states) : Symbol[],
            params=!isempty(params) ? tosymbol.(params) : Symbol[]
        )
        flux_name = isnothing(name) ? Symbol("##state_flux#", hash(infos)) : name
        return new{flux_name,typeof(eqs),typeof(infos)}(eqs, infos)
    end
end

"""
    @stateflux [name] eq

A macro to conveniently create a `StateFlux` component from a single equation.

# Usage
The left side of the equation is the state variable, and the right side is the expression for its rate of change. The `~` operator is typically used.

```julia
@variables S, P, Q
@stateflux :storage_change S ~ P - Q
```
"""
macro stateflux(args...)
    name = length(args) == 1 ? nothing : args[1]
    eq_expr = length(args) == 1 ? args[1] : args[2]
    for var_name in extract_variables(eq_expr)
        if !@isdefined(var_name)
            expr_str = string(eq_expr)
            return :(error("Undefined variable '", $(string(var_name)), "' detected in expression: `", $expr_str, "`"))
        end
    end
    return esc(:(StateFlux(exprs=[$eq_expr], name=$name)))
end

(::StateFlux)(::AbstractArray, ::ComponentVector; kwargs...) = @error "State Flux cannot run directly, please using HydroFlux to run"

"""
    NeuralFlux{C, F, NT} <: AbstractNeuralFlux

Represents a flux component driven by a neural network.

It wraps a `Lux.AbstractLuxLayer` and connects it to symbolic variables for integration into a hydrological model.

$(FIELDS)
"""
struct NeuralFlux{C,CF,NF,NT} <: AbstractNeuralFlux
    "neural flux name"
    name::Symbol
    "chain of the neural network"
    chain::C
    "Compiled function that calculates the flux using the neural network"
    chain_func::CF
    "input normalizatio functions"
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
        new{typeof(chain),typeof(nn_func),typeof(norm),typeof(infos)}(flux_name, chain, nn_func, norm, infos)
    end
end

"""
    @neuralflux [name] eq

A macro to conveniently create a `NeuralFlux` from an equation.

# Usage

The macro takes an optional name and an equation of the form `output ~ chain(inputs)`.

# Examples

```julia
@variables x, y, z, z₁, z₂
chain = Chain(Dense(2 => 10, relu), Dense(10 => 1), name=:my_net)

# Single output
flux1 = @neuralflux z ~ chain([x, y])

# With an optional name
flux2 = @neuralflux :my_flux z ~ chain([x, y])

# Multiple outputs
chain2 = Chain(Dense(2 => 16, relu), Dense(16 => 2), name=:multi_net)
flux3 = @neuralflux [z₁, z₂] ~ chain2([x, y])
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
    lhs, rhs = eq_expr.args[2], eq_expr.args[3]  # Output variable(s) and Chain info expression
    @assert rhs.head == :call "The right-hand side of `~` must be a function call, e.g., my_chain([x, y])"
    @assert length(rhs.args) >= 2 "The chain call must have at least one argument for the inputs."

    chain_expr, inputs_expr = rhs.args[1], rhs.args[2]
    return esc(quote
        local outputs = $lhs isa AbstractVector ? $lhs : [$lhs]
        NeuralFlux($inputs_expr, outputs, $chain_expr; name=$(name))
    end)
end

function (flux::NeuralFlux)(input::AbstractArray{T,2}, params::ComponentVector; kwargs...) where {T}
    nn_params = params[:nns][get_nn_names(flux)[1]]
    flux.chain_func(flux.norm_func(input), nn_params)
end

function (flux::NeuralFlux)(input::AbstractArray{T,3}, params::ComponentVector; kwargs...) where {T}
    nn_params = params[:nns][get_nn_names(flux)[1]]
    norm_input = flux.norm_func(input)
    stack(ntuple(i -> flux.chain_func(norm_input[:, i, :], nn_params), size(input)[2]), dims=2)
end