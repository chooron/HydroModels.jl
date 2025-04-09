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

    equations_expr = if Meta.isexpr(eqs_expr, :block)
        eqs_list = filter(arg -> !(arg isa LineNumberNode) && !Meta.isexpr(arg, :line), eqs_expr.args)
        Expr(:vect, eqs_list...)
    else
        Expr(:vect, eqs_expr)
    end

    return esc(quote
        let
            equations = $equations_expr
            processed_eqs = map(equations) do eq
                if eq isa Equation
                    eq
                elseif Meta.isexpr(eq, :(=)) || Meta.isexpr(eq, :(~))
                    Symbolics.Equation(eq.args[1], eq.args[2])
                else
                    error("Expected equation (using = or ~), got: $eq")
                end
            end

            lhs_terms = [eq.lhs for eq in processed_eqs]
            rhs_terms = [eq.rhs for eq in processed_eqs]
            outputs = Num.(lhs_terms)

            all_vars = Num.(mapreduce(get_variables, union, rhs_terms, init=Set{Num}()))
            inputs = Num.(filter(x -> !ModelingToolkit.isparameter(x), collect(all_vars)))
            params = Num.(filter(x -> ModelingToolkit.isparameter(x), collect(all_vars)))
            inputs = setdiff(inputs, outputs)

            HydroFlux(inputs, outputs, params, exprs=Num.(rhs_terms), name=$(name))
        end
    end)
end


(chain::LuxCore.AbstractLuxLayer)(inputs::Vector{Num}) = (chain=chain, inputs=inputs, name=hasproperty(chain, :name) ? chain.name : nothing)

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

    # Get left and right sides of the tilde
    lhs = eqs_expr.args[2]  # Output variable(s)
    rhs = eqs_expr.args[3]  # Chain info expression

    # Create the expression for the macro expansion
    return esc(quote
        let
            # Extract outputs (already Num type)
            outputs = $lhs isa Vector ? $lhs : [$lhs]
            # Get the chain info
            chain_info = $rhs
            # Create the NeuralFlux
            NeuralFlux(
                chain_info.inputs,
                outputs,
                chain_info.chain;
                name=$(name),
                chain_name=chain_info.name
            )
        end
    end)
end

"""
    @stateflux_build(name, eq)

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
        let
            eq = Symbolics.Equation($lhs, $rhs)
            state = Num(eq.lhs)

            all_vars = Num.(get_variables(eq.rhs))
            inputs = Num.(filter(x -> !ModelingToolkit.isparameter(x), collect(all_vars)))
            params = Num.(filter(x -> ModelingToolkit.isparameter(x), collect(all_vars)))
            inputs = setdiff(inputs, state)

            StateFlux(inputs, only(state), params, expr=Num(eq.rhs), name=$(name))
        end
    end)
end

"""
    @hydrobucket name begin
        fluxes = [...]
        dfluxes = [...]
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
macro hydrobucket(name, expr)
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
    return esc(quote
        let
            fluxes = $fluxes_expr
            dfluxes = isnothing($dfluxes_expr) ? [] : $dfluxes_expr
            HydroBucket(name=$(name), fluxes=fluxes, dfluxes=dfluxes)
        end
    end)
end

"""
    @hydroroute name begin
        fluxes = [...]
        dfluxes = [...]
        proj_func = f(x)
    end

构建一个 HydroRoute 对象的宏。

# Arguments
- `name`: 路由的名称，一个 Symbol
- `expr`: 包含以下组件定义的代码块:
  - `fluxes`: HydroFlux 或 NeuralFlux 对象数组，定义路由的流量计算
  - `dfluxes`: StateFlux 对象数组，定义状态变量的变化
  - `proj_func`: 投影函数，类型为 Function，用于确保状态变量保持在有效范围内

# Example
```julia
route = @hydroroute :route1 begin
    fluxes = [
        @hydroflux a ~ k1 * b - k2 * c
        @hydroflux d ~ k1 * b - k2 * c
    ]
    dfluxes = [
        @stateflux c ~ b - a - d
    ]
    proj_func = x -> max(0, x)
end
```
"""
macro hydroroute(name, expr)
    @assert Meta.isexpr(expr, :block) "Expected a begin...end block after route name"
    fluxes_expr, dfluxes_expr, proj_func_expr = nothing, nothing, nothing
    for assign in filter(x -> !(x isa LineNumberNode), expr.args)
        @assert Meta.isexpr(assign, :(=)) "Expected assignments in the form 'fluxes = begin...end', 'dfluxes = begin...end', and 'proj_func = f(x)'"
        lhs, rhs = assign.args
        if lhs == :fluxes
            @assert Meta.isexpr(rhs, :block) "Expected 'fluxes' to be defined in a begin...end block"
            fluxes_expr = Expr(:vect, filter(x -> !(x isa LineNumberNode), rhs.args)...)
        elseif lhs == :dfluxes
            @assert Meta.isexpr(rhs, :block) "Expected 'dfluxes' to be defined in a begin...end block"
            dfluxes_expr = Expr(:vect, filter(x -> !(x isa LineNumberNode), rhs.args)...)
        elseif lhs == :proj_func
            proj_func_expr = rhs
        else
            error("Unknown assignment: $lhs. Expected 'fluxes', 'dfluxes', or 'proj_func'")
        end
    end
    @assert !isnothing(fluxes_expr) && !isnothing(dfluxes_expr) && !isnothing(proj_func_expr) "'fluxes', 'dfluxes', and 'proj_func' must all be specified"
    return esc(quote
        let
            fluxes = $fluxes_expr
            dfluxes = $dfluxes_expr
            proj_func = $proj_func_expr
            HydroRoute(
                rfluxes=fluxes,
                dfluxes=dfluxes,
                proj_func=proj_func,
                name=$(name)
            )
        end
    end)
end

"""
    @hydromodel name begin
        component1
        component2
        ...
    end

Creates a HydroModel with the specified name and components.

# Arguments
- `name`: Symbol for the model name
- Components can be:
  - HydroBucket instances
  - Flux definitions (using @hydroflux, @neuralflux, etc.)
  - Other model components

# Example
```julia
@hydromodel :model1 begin
    bucket1
    @neuralflux :flux1 [y] ~ chain([x1, x2])
    bucket2
end
```
"""
macro hydromodel(name, expr)
    @assert Meta.isexpr(expr, :block) "Expected a begin...end block after model name"
    
    # Filter out LineNumberNodes and get components
    components = filter(x -> !(x isa LineNumberNode), expr.args)

    # Create the model
    return esc(quote
        let
            # Create the model with all components
            HydroModel(; name=$(name), components=[$(components...)])
        end
    end)
end
