module HydroModelsOrdinaryDiffEqFunctionMapExt

using HydroModels
using SciMLBase
using OrdinaryDiffEqFunctionMap: FunctionMap

"""
    hydrosolve(::Val{HydroModels.DiscreteSolver}, du_func, params, initstates, timeidx, config)

Solve a discrete recurrence with `FunctionMap`.  This method is available
only after loading `OrdinaryDiffEqFunctionMap`; loading `OrdinaryDiffEq` or
`DifferentialEquations` alone does not provide it.
"""
function HydroModels.hydrosolve(
    ::Val{HydroModels.DiscreteSolver},
    du_func,
    params,
    initstates,
    timeidx,
    config,
)
    device = HydroModels.get_config_value(config, :device, identity)
    solve_alg = HydroModels.get_config_value(config, :solve_alg, FunctionMap{true}())
    sense_alg = HydroModels.get_config_value(config, :sense_alg, nothing)
    solve_cb = HydroModels.get_config_value(config, :solve_cb, nothing)

    function map_func!(du, u, p, t)
        du .= du_func(u, p, t)
        return nothing
    end

    prob = DiscreteProblem(
        map_func!,
        initstates,
        (first(timeidx), last(timeidx)),
        params,
    )

    solve_kwargs = (; saveat=timeidx)
    if !isnothing(sense_alg)
        solve_kwargs = merge(solve_kwargs, (; sensealg=sense_alg))
    end
    if !isnothing(solve_cb)
        solve_kwargs = merge(solve_kwargs, (; callback=solve_cb))
    end

    sol = solve(prob, solve_alg; solve_kwargs...)
    return device(Array(sol))
end

end
