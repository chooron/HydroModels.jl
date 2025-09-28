module HydroModelsOptimizationExt

using SciMLBase
using Optimization
using HydroModels
using HydroModelCore
using ComponentArrays

function SciMLBase.OptimizationProblem(
    component::AbstractComponent,
    input::AbstractArray{T,2},
    target::AbstractVector;
    kwargs...
)::OptimizationProblem where {T}
    loss_func = get(kwargs, :loss_func, (y, y_hat) -> sum((y .- y_hat) .^ 2) ./ length(y))
    warm_up = get(kwargs, :warm_up, 1)
    adtype = get(kwargs, :adtype, nothing)
    interp = get(kwargs, :interpolator, Val(HydroModels.DirectInterpolation))
    timeidx = get(kwargs, :timeidx, collect(1:size(input, 2)))
    solver = get(kwargs, :solver, HydroModels.MutableSolver)
    solve_alg = get(kwargs, :solve_alg, nothing)
    sense_alg = get(kwargs, :sense_alg, nothing)
    solve_cb = get(kwargs, :solve_cb, nothing)

    param_names, state_names = HydroModels.get_param_names(component), HydroModels.get_state_names(component)
    default_initstates = get(kwargs, :default_initstates, ComponentVector(NamedTuple{Tuple(state_names)}(zeros(eltype(input), length(state_names)))))
    default_pas = get(kwargs, :default_pas, ComponentVector(params=NamedTuple{Tuple(param_names)}(ones(eltype(input), length(param_names)))))
    pas_axes = getaxes(default_pas)
    lb_pas, ub_pas = get(kwargs, :lb_pas, nothing), get(kwargs, :ub_pas, nothing)
    if isnothing(lb_pas) && isnothing(ub_pas)
        prob_kwargs = ()
    else
        prob_kwargs = (lb=lb_pas, ub=ub_pas)
    end
    config = (solver=solver, interpolator=interp, solve_alg=solve_alg, sense_alg=sense_alg, solve_cb=solve_cb)

    function objective(p, c)
        output = component(input, ComponentVector(p, pas_axes), c; timeidx=timeidx, initstates=default_initstates)
        loss_func(target[warm_up:end], output[end, warm_up:end])
    end
    if isnothing(adtype)
        optfunc = OptimizationFunction(objective)
    else
        optfunc = OptimizationFunction(objective, adtype)
    end
    prob = OptimizationProblem(optfunc, Vector(default_pas), config; prob_kwargs...)
    return prob
end
end
