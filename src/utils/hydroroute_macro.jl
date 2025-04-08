"""
    @hydroroute name begin
        fluxes = [...]
        dfluxes = [...]
        proj_func = f(x)
    end

构建一个 HydroRoute 对象的宏。

# Arguments
- `name`: 路由的名称，一个 Symbol
- `expr`: 包含 fluxes、dfluxes 和 proj_func 定义的代码块

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
    if !Meta.isexpr(expr, :block)
        error("Expected a begin...end block after route name")
    end

    # Filter out LineNumberNodes and get assignments
    assignments = filter(x -> !(x isa LineNumberNode), expr.args)

    # Initialize containers
    fluxes_expr = nothing
    dfluxes_expr = nothing
    proj_func_expr = nothing

    # Process assignments
    for assign in assignments
        if !Meta.isexpr(assign, :(=))
            error("Expected assignments in the form 'fluxes = [...]', 'dfluxes = [...]', and 'proj_func = f(x)'")
        end

        lhs, rhs = assign.args
        if lhs == :fluxes
            fluxes_expr = rhs
        elseif lhs == :dfluxes
            dfluxes_expr = rhs
        elseif lhs == :proj_func
            proj_func_expr = rhs
        else
            error("Unknown assignment: $lhs. Expected 'fluxes', 'dfluxes', or 'proj_func'")
        end
    end

    # Ensure all required components are provided
    if isnothing(fluxes_expr) || isnothing(dfluxes_expr) || isnothing(proj_func_expr)
        error("'fluxes', 'dfluxes', and 'proj_func' must all be specified")
    end

    # Create the route
    return esc(quote
        let
            # Evaluate the flux arrays and projection function
            fluxes = $fluxes_expr
            dfluxes = $dfluxes_expr
            proj_func = $proj_func_expr

            # Validate flux types
            all(f -> f isa Union{HydroFlux,NeuralFlux}, fluxes) ||
                error("All elements in fluxes must be HydroFlux or NeuralFlux objects")
            all(f -> f isa StateFlux, dfluxes) ||
                error("All elements in dfluxes must be StateFlux objects")
            proj_func isa Function ||
                error("proj_func must be a function")

            # Create the route
            HydroRoute(
                rfluxes=fluxes,
                dfluxes=dfluxes,
                proj_func=proj_func,
                name=$(name)
            )
        end
    end)
end
