module HydroModelOptimisers

using TOML
using Dates
using DataFrames
using ProgressMeter
using IterTools: ncycle
using NamedTupleTools
using ComponentArrays
using ComponentArrays: indexmap, getval
using Optimization

const version = VersionNumber(TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["version"])

abstract type AbstractHydroOptimizer end

@kwdef struct HydroOptimizer{C,S,A} <: AbstractHydroOptimizer
    component::C
    solve_alg::S
    adtype::A = nothing
    maxiters::Int = 1000
    warmup::Int = 100
    loss_func::Function = (obs, sim) -> sum((obs .- sim) .^ 2) / length(obs)
    callback_func::Function = (loss_recorder) -> get_callback_func(Progress(maxiters, desc="Training..."), loss_recorder)
    objective_func::Function = get_hydro_objective(component, loss_func, warmup)
end

function get_hydro_objective(component, loss_func, warmup)
    #* Constructing the objective function for optimization
    function objective(x::AbstractVector{T}, p) where {T}
        #* Optimization arguments: hydro component, input data, time index, ode solver,
        #*                         tunable parameters axes and default model params
        inputs, targets, kwargs, tunable_axes, default_model_pas = p
        #* Use merge_ca to replace the tunable parameters inner the model parameters
        tmp_tunable_pas = ComponentArray(x, tunable_axes)
        tmp_pas = update_ca(default_model_pas, tmp_tunable_pas)
        loss = mean(map(eachindex(inputs, targets, kwargs)) do i
            tmp_pred = component(inputs[i], tmp_pas, kwargs[i])
            tmp_loss = mean([loss_func(targets[i][key][warmup:end], tmp_pred[key][warmup:end]) for key in keys(targets[i])])
            tmp_loss
        end)
        loss
    end

    return objective
end


function (opt::HydroOptimizer{C,S,A})(
    input::Vector,
    target::Vector;
    tunable_pas::ComponentVector,
    const_pas::ComponentVector,
    run_kwargs::Vector=fill(Dict(), length(input)),
    kwargs...
) where {C,S,A}
    loss_recorder = NamedTuple[]
    callback_func = opt.callback_func(loss_recorder)
    tunable_axes = getaxes(tunable_pas)
    default_model_pas = ComponentArray(merge_recursive(NamedTuple(tunable_pas), NamedTuple(const_pas)))
    prob_args = (input, target, run_kwargs, tunable_axes, default_model_pas)
    #* Constructing and solving optimization problems
    @info "The size of tunable parameters is $(length(tunable_pas))"
    if A isa Nothing
        optf = Optimization.OptimizationFunction(opt.objective_func)
        lb = get(kwargs, :lb, zeros(length(tunable_pas)))
        ub = get(kwargs, :ub, ones(length(tunable_pas)) .* 100)
        optprob = Optimization.OptimizationProblem(optf, collect(tunable_pas), prob_args, lb=lb, ub=ub)
    else
        optf = Optimization.OptimizationFunction(opt.objective_func, opt.adtype)
        optprob = Optimization.OptimizationProblem(optf, collect(tunable_pas), prob_args)
    end
    sol = Optimization.solve(optprob, opt.solve_alg, callback=callback_func, maxiters=opt.maxiters)
    opt_pas = update_ca(default_model_pas, ComponentVector(sol.u, tunable_axes))
    if get(kwargs, :return_loss_df, false)
        loss_df = DataFrame(loss_recorder)
        return opt_pas, loss_df
    else
        return opt_pas
    end
end

@kwdef struct BatchOptimizer{C,S} <: AbstractHydroOptimizer
    component::C
    solve_alg::S
    maxiters::Int = 100
    warmup::Int = 100
    adtype::AbstractADType = AutoForwardDiff()
    loss_func::Function = (obs, sim) -> sum((obs .- sim) .^ 2) / length(obs)
    callback_func::Function = (batch_size, loss_recorder) -> get_batch_callback_func(batch_size, loss_recorder)
    objective_func::Function = get_batch_objective(wrapper(component), loss_func, warmup)
end

function get_batch_objective(component, loss_func, warmup)

    #* Constructing the objective function for optimization
    function objective(x::AbstractVector{T}, p, input, target, config) where {T}
        #* Optimization arguments: hydro component, input data, time index, ode solver,
        #*                         tunable parameters axes and default model params
        run_kwargs, tunable_axes, default_model_pas = p
        #* Use merge_ca to replace the tunable parameters inner the model parameters
        tmp_tunable_pas = ComponentArray(x, tunable_axes)
        tmp_pas = update_ca(default_model_pas, tmp_tunable_pas)
        tmp_pred = component(input, tmp_pas; config=config, run_kwargs...)
        loss = mean([loss_func(target[key][warmup:end], tmp_pred[key][warmup:end]) for key in keys(target)])
        loss
    end

    return objective
end

function (opt::BatchOptimizer{C,S})(
    input::Vector,
    target::Vector;
    tunable_pas::ComponentVector,
    const_pas::ComponentVector,
    config::Vector=fill(NamedTuple(), length(input)),
    kwargs...
) where {C,S}
    loss_recorder = NamedTuple[]
    callback_func = opt.callback_func(length(input), loss_recorder)
    default_model_pas = ComponentArray(merge_recursive(NamedTuple(tunable_pas), NamedTuple(const_pas)))

    run_kwargs = vcat([merge(get(kwargs, :run_kwargs, NamedTuple()), (reset_states=true,))],
        fill(get(kwargs, :run_kwargs, NamedTuple()), length(input) - 1))
    return_loss_df = get(kwargs, :return_loss_df, false)

    tunable_axes = getaxes(tunable_pas)

    #* prepare the batch data
    train_batch = [(input_i, target_i, cfg_i, run_kw_i) for (input_i, target_i, cfg_i, run_kw_i) in zip(input, target, config, run_kwargs)]
    #* Construct default model parameters based on tunbale parameters and constant parameters for subsequent merging
    default_model_pas = ComponentArray(merge_recursive(NamedTuple(tunable_pas), NamedTuple(const_pas)))

    prob_args = (tunable_axes, default_model_pas, ncycle(train_batch, opt.maxiters))

    #* Constructing and solving optimization problems
    optf = Optimization.OptimizationFunction(opt.objective_func, opt.adtype)
    @info "The size of tunable parameters is $(length(tunable_pas))"
    optprob = Optimization.OptimizationProblem(optf, collect(tunable_pas), prob_args)
    # TODO batch training is not supported for after Optimization.jl v0.4+
    sol = Optimization.solve(optprob, opt.solve_alg, ncycle(train_batch, opt.maxiters), callback=callback_func)
    #* Returns the optimized model parameters
    opt_pas = update_ca(default_model_pas, ComponentVector(sol.u, tunable_axes))
    if return_loss_df
        loss_df = DataFrame(loss_recorder)
        return opt_pas, loss_df
    else
        return opt_pas
    end
end


end # module HydroModelOptimisers
