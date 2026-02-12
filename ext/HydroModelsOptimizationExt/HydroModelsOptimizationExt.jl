module HydroModelsOptimizationExt

using SciMLBase
using Optimization
using HydroModels
using HydroModelCore
using ComponentArrays
using Random
using LuxCore
using Statistics

# Helper function to update ComponentArray with another ComponentArray
function update_ca(ca::ComponentArray{T1}, ca2::ComponentArray{T2}) where {T1,T2}
    ax = getaxes(ca)
    ax2 = getaxes(ca2)
    vks = valkeys(ax[1])
    vks2 = valkeys(ax2[1])
    _p = Vector{T2}()
    for vk in vks
        if length(getaxes(ca[vk])) > 0
            _p = vcat(_p, collect(update_ca(ca[vk], vk in vks2 ? getproperty(ca2, vk) : ComponentVector())))
        else
            if vk in vks2
                _p = vcat(_p, ca2[vk])
            else
                _p = vcat(_p, ca[vk])
            end
        end
    end
    ComponentArray(_p, ax)
end

# Metric functions for hydrological calibration
function _kge(obs::AbstractVector, sim::AbstractVector)
    r = cor(obs, sim)
    α = std(sim) / std(obs)
    β = mean(sim) / mean(obs)
    return 1 - sqrt((r - 1)^2 + (α - 1)^2 + (β - 1)^2)
end

function _nse(obs::AbstractVector, sim::AbstractVector)
    return 1 - sum((obs .- sim) .^ 2) / sum((obs .- mean(obs)) .^ 2)
end

function _log_kge(obs::AbstractVector, sim::AbstractVector)
    log_obs = log.(obs .+ 1e-6)
    log_sim = log.(sim .+ 1e-6)
    return _kge(log_obs, log_sim)
end

function _mse(obs::AbstractVector, sim::AbstractVector)
    return sum((obs .- sim) .^ 2) / length(obs)
end

function _get_metric_function(metric::String)
    metric_lower = lowercase(metric)
    if metric_lower == "kge"
        return (obs, sim) -> -_kge(obs, sim)  # Minimize negative KGE
    elseif metric_lower == "nse"
        return (obs, sim) -> -_nse(obs, sim)  # Minimize negative NSE
    elseif metric_lower == "logkge"
        return (obs, sim) -> -_log_kge(obs, sim)  # Minimize negative LogKGE
    elseif metric_lower == "mse"
        return _mse
    else
        error("Unknown metric: $metric. Supported metrics: KGE, NSE, LogKGE, MSE")
    end
end

function SciMLBase.OptimizationProblem(
    component::AbstractComponent,
    input::AbstractArray{T,2},
    target::AbstractVector;
    kwargs...
)::OptimizationProblem where {T}
    # Get metric or loss function
    metric = get(kwargs, :metric, nothing)
    if !isnothing(metric)
        loss_func = _get_metric_function(metric)
    else
        loss_func = get(kwargs, :loss_func, (y, y_hat) -> sum((y .- y_hat) .^ 2) ./ length(y))
    end

    warm_up = get(kwargs, :warm_up, 1)
    adtype = get(kwargs, :adtype, nothing)
    interp = get(kwargs, :interpolator, Val(HydroModels.ConstantInterpolation))
    timeidx = get(kwargs, :timeidx, collect(1:size(input, 2)))
    solver = get(kwargs, :solver, HydroModels.MutableSolver)
    solve_alg = get(kwargs, :solve_alg, nothing)
    sense_alg = get(kwargs, :sense_alg, nothing)
    solve_cb = get(kwargs, :solve_cb, nothing)

    param_names, state_names = HydroModels.get_param_names(component), HydroModels.get_state_names(component)
    nn_names = HydroModels.get_nn_names(component)
    default_initstates = get(kwargs, :default_initstates, ComponentVector(NamedTuple{Tuple(state_names)}(zeros(eltype(input), length(state_names)))))

    # Get initial parameters
    # User can provide initial_params, otherwise use get_initial_params
    initial_params = get(kwargs, :initial_params, nothing)
    if isnothing(initial_params)
        default_pas = HydroModels.get_initial_params(component; eltype=eltype(input))
    else
        default_pas = initial_params isa ComponentVector ? initial_params : ComponentVector(initial_params)
    end

    # Handle fixed parameters (can be ComponentVector with params and/or nns)
    fixed_params = get(kwargs, :fixed_params, nothing)

    # Prepare parameters and bounds based on whether fixed_params is provided
    local objective, lb_pas, ub_pas
    
    if !isnothing(fixed_params)
        # Convert fixed_params to ComponentVector if it's a NamedTuple
        fixed_pas = fixed_params isa ComponentVector ? fixed_params : ComponentVector(fixed_params)

        # Determine which params need calibration
        # Extract fixed param names from the params layer (if exists)
        fixed_param_names = Symbol[]
        fixed_axes = getaxes(fixed_pas)[1]
        if :params in keys(fixed_axes)
            fixed_param_names = collect(keys(fixed_pas[:params]))
        end

        # Filter out fixed parameters for calibration
        calibratable_params = filter(p -> !(p in fixed_param_names), param_names)

        # Create calibratable parameter vector (only non-fixed params)
        calibratable_pas = if !isempty(calibratable_params)
            ComponentVector(params=NamedTuple{Tuple(calibratable_params)}(ones(eltype(input), length(calibratable_params))))
        else
            # All params are fixed, create empty calibratable params
            ComponentVector()
        end

        # Extract bounds for calibratable parameters only
        lb_pas_full = get(kwargs, :lb_pas, nothing)
        ub_pas_full = get(kwargs, :ub_pas, nothing)

        if !isnothing(lb_pas_full) && !isnothing(ub_pas_full) && !isempty(calibratable_params)
            # Filter bounds for calibratable parameters
            # lb_pas_full and ub_pas_full are already filtered to only include calibratable params
            # So we just use them directly
            lb_pas = lb_pas_full
            ub_pas = ub_pas_full
        else
            lb_pas, ub_pas = nothing, nothing
        end

        pas_axes = getaxes(calibratable_pas)
        default_pas = calibratable_pas

        # Define objective function with fixed parameters
        objective = (p, c) -> begin
            full_pas = if !isempty(calibratable_params)
                cal_pas = ComponentVector(p, pas_axes)
                # Update complete parameters with calibratable values, then with fixed values
                update_ca(update_ca(default_pas, cal_pas), fixed_pas)
            else
                # All params are fixed, just use fixed params
                update_ca(default_pas, fixed_pas)
            end
            output = component(input, full_pas, c; timeidx=timeidx, initstates=default_initstates)
            loss_func(target[warm_up:end], output[end, warm_up:end])
        end
    else
        # Original behavior: calibrate all parameters
        pas_axes = getaxes(default_pas)
        lb_pas, ub_pas = get(kwargs, :lb_pas, nothing), get(kwargs, :ub_pas, nothing)

        # Define objective function without fixed parameters
        objective = (p, c) -> begin
            output = component(input, ComponentVector(p, pas_axes), c; timeidx=timeidx, initstates=default_initstates)
            loss_func(target[warm_up:end], output[end, warm_up:end])
        end
    end

    if isnothing(lb_pas) && isnothing(ub_pas)
        prob_kwargs = ()
    else
        prob_kwargs = (lb=lb_pas, ub=ub_pas)
    end

    # Create HydroConfig with only supported parameters
    hydro_config = HydroModels.HydroConfig(
        solver=solver,
        interpolator=interp,
        timeidx=timeidx
    )

    # Create config tuple with HydroConfig and optimization-specific parameters
    config = (hydro_config=hydro_config, solve_alg=solve_alg, sense_alg=sense_alg, solve_cb=solve_cb)

    if isnothing(adtype)
        optfunc = OptimizationFunction(objective)
    else
        optfunc = OptimizationFunction(objective, adtype)
    end
    prob = OptimizationProblem(optfunc, Vector(default_pas), config; prob_kwargs...)
    return prob
end

end
