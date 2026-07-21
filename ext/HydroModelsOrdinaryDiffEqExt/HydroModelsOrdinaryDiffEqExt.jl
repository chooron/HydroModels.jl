module HydroModelsOrdinaryDiffEqExt

using HydroModels
using SciMLBase
using OrdinaryDiffEq: Tsit5

@inline _constant_tstops(::Type{HydroModels.ConstantInterpolation}, timeidx) =
    length(timeidx) > 2 ? collect(@view timeidx[2:end-1]) : eltype(timeidx)[]

@inline _constant_tstops(::HydroModels.ConstantInterpolation, timeidx) =
    length(timeidx) > 2 ? collect(@view timeidx[2:end-1]) : eltype(timeidx)[]

@inline _constant_tstops(_, timeidx) = nothing

"""
    hydrosolve(::Val{HydroModels.ODESolver}, du_func, params, initstates, timeidx, config)

Solve a continuous ODE using the algorithm and optional callback/sensitivity
settings supplied by `config`.  This extension intentionally depends on the
topic package `OrdinaryDiffEq`; users selecting another solver family must
load that family's sublibrary themselves.
"""
function HydroModels.hydrosolve(
    ::Val{HydroModels.ODESolver},
    du_func,
    params,
    initstates,
    timeidx,
    config,
)
    device = HydroModels.get_config_value(config, :device, identity)
    solve_alg = HydroModels.get_config_value(config, :solve_alg, Tsit5())
    sense_alg = HydroModels.get_config_value(config, :sense_alg, nothing)
    solve_cb = HydroModels.get_config_value(config, :solve_cb, nothing)
    interp_type = HydroModels.get_config_value(
        config,
        :interpolator,
        HydroModels.ConstantInterpolation,
    )

    function ode_func!(du, u, p, t)
        du .= du_func(u, p, t)
        return nothing
    end

    prob = ODEProblem{true}(
        ode_func!,
        initstates,
        (first(timeidx), last(timeidx)),
        params,
    )

    tstops = _constant_tstops(interp_type, timeidx)

    solve_kwargs = (; saveat=timeidx)
    if !isnothing(sense_alg)
        solve_kwargs = merge(solve_kwargs, (; sensealg=sense_alg))
    end
    if !isnothing(tstops)
        solve_kwargs = merge(solve_kwargs, (; tstops))
    end
    if !isnothing(solve_cb)
        solve_kwargs = merge(solve_kwargs, (; callback=solve_cb))
    end

    sol = solve(prob, solve_alg; solve_kwargs...)
    return device(Array(sol))
end

end
