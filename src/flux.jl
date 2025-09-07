"""
    HydroFlux{E,F,I} <: AbstractHydroFlux

Represents a simple flux component with mathematical formulas in a hydrological model.

The component automatically determines inputs, outputs, and parameters from the provided expressions. It then builds a callable function to perform the calculations efficiently.

# Arguments
- `exprs::Vector{Equation}`: A vector of equations defining the flux calculations. The left-hand side of each equation is treated as an output variable, and the right-hand side is the expression used to compute it.
- `name::Union{Symbol,Nothing}=nothing`: An optional identifier for the flux component. If not provided, a unique name is generated automatically.

# Fields
- `name::Symbol`: The identifier for the flux.
- `exprs::Vector{Num}`: Vector of the right-hand-side mathematical expressions from the input equations.
- `func::Function`: A compiled function generated from the expressions for efficient computation. The function is designed to work with arrays of input data.
- `infos::NamedTuple`: Metadata about the component, containing the derived names of `inputs`, `outputs`, and `params`.
"""
struct HydroFlux{E,F,I} <: AbstractHydroFlux
    "flux name"
    name::Symbol
    "Vector of expressions describing the formulas for output variables"
    exprs::E
    "Compiled function that calculates the flux"
    func::F
    "Metadata about the flux, including input, output, and parameter names"
    infos::I

    function HydroFlux(;
        exprs::E, name::Optional{Symbol}=nothing,
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

        return new{typeof(eqs),typeof(flux_func),typeof(infos)}(
            flux_name, eqs, flux_func, infos
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
        Expr(:tuple, filter(arg -> !(arg isa LineNumberNode) && !Meta.isexpr(arg, :line), eqs_expr.args)...)
    elseif Meta.isexpr(eqs_expr, :for)
        loop_var = eqs_expr.args[1].args[1]
        range_expr = eqs_expr.args[1].args[2]
        range_val = if Meta.isexpr(range_expr, :call) && range_expr.args[1] == :(:)
            if length(range_expr.args) == 3
                range_expr.args[2]:range_expr.args[3]
            elseif length(range_expr.args) == 4
                range_expr.args[2]:range_expr.args[3]:range_expr.args[4]
            end
        else
            eval(range_expr)
        end
        loop_body = eqs_expr.args[2]
        if !Meta.isexpr(loop_body, :block)
            loop_body = Expr(:block, loop_body)
        end
        equations = filter(x -> !(x isa LineNumberNode) && !Meta.isexpr(x, :line), loop_body.args)
        all_equations = []
        for i_val in range_val
            for eq in equations
                new_eq = deepcopy(eq)
                replace_loop_var!(new_eq, loop_var, i_val)
                push!(all_equations, new_eq)
            end
        end
        Expr(:tuple, all_equations...)
    else
        Expr(:tuple, eqs_expr)
    end

    for var_name in extract_variables(vect_eqs_expr)
        if !@isdefined(var_name)
            expr_str = string(vect_eqs_expr)
            return :(error("Undefined variable '", $(string(var_name)), "' detected in expression: `", $expr_str, "`"))
        end
    end
    return esc(:(HydroFlux(exprs=$vect_eqs_expr, name=$name)))
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
    build_flux_func(inputs::Vector{Num}, outputs::Vector{Num}, params::Vector{Num}, exprs::Vector{Num})

Generates a runtime function for flux calculations based on symbolic expressions.
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


"""
    StateFlux{N,E,I} <: AbstractStateFlux

Represents a state flux component that symbolically defines the rate of change for a state variable in a hydrological model. The component's name is encoded in the type parameter `N` for enhanced type stability.

This component is primarily a declarative structure used within a larger model (like a `HydroBucket`) to define the system's differential equations. It does not perform calculations directly.

# Arguments
- `exprs::Vector{Equation}`: A vector of equations defining the state change. The left-hand side represents the state variable, and the right-hand side is the expression for its derivative (rate of change).
- `name::Union{Symbol,Nothing}=nothing`: An optional identifier for the state flux. If not provided, a name is generated based on the state variable.

# Fields
- `exprs::Vector{Num}`: A vector containing the right-hand-side expressions from the input equations, which define the state's rate of change.
- `infos::NamedTuple`: Metadata about the component, including the derived names of `inputs`, `states`, and `params` extracted from the expressions.
"""
function (flux::HydroFlux)(input::AbstractArray{T,2}, params::ComponentVector; kwargs...) where {T}
    stack(flux.func(eachslice(input, dims=1), params), dims=1)
end

function (flux::HydroFlux)(input::AbstractArray{T,3}, params::ComponentVector; kwargs...) where {T}
    ptyidx = get(kwargs, :ptyidx, collect(1:size(input, 2)))
    expand_params = expand_component_params(params, get_param_names(flux), ptyidx)
    output = flux.func(eachslice(input, dims=1), expand_params)
    stack(output, dims=1)
end

"""
    @stateflux(name, eq)

A macro to conveniently create a `StateFlux` component from a single equation.

The left side of the equation is interpreted as the state variable, and the right side is the expression defining its rate of change.

The standard convention is to use `~` to define the relationship (e.g., `S ~ P - E`), consistent with `ModelingToolkit` for differential equations.

# Arguments
- `name`: An optional `Symbol` to name the state flux component.
- `eq`: A single equation that defines the change in a state variable.
"""
struct StateFlux{N,E,I} <: AbstractStateFlux
    "flux expressions to descripe the formula of the state variable"
    exprs::E
    "bucket information: keys contains: input, output, param, state"
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
    @stateflux(name, eq)

A macro to conveniently create a `StateFlux` component from a single equation.

The left side of the equation is interpreted as the state variable, and the right side is the expression defining its rate of change.

The standard convention is to use `~` to define the relationship (e.g., `S ~ P - E`), consistent with `ModelingToolkit` for differential equations.

# Arguments
- `name`: An optional `Symbol` to name the state flux component.
- `eq`: A single equation that defines the change in a state variable.
"""
macro stateflux(args...)
    # Process arguments
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