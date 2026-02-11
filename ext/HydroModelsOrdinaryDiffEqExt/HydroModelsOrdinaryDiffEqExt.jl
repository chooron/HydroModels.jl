module HydroModelsOrdinaryDiffEqExt

using SciMLBase
using OrdinaryDiffEq
using SciMLSensitivity
using HydroModels

"""
    hydrosolve(::Val{HydroModels.ODESolver}, du_func, pas, initstates, timeidx, config)

Solve a continuous-time Ordinary Differential Equation (ODE) problem.

This method dispatches on `ODESolver` to set up and solve an `ODEProblem`
using the algorithms specified in the `config`.

# Arguments
- `du_func`: The derivative function `du/dt = f(u, p, t)`. (Note: Currently `ode_func!` is used internally).
- `pas`: Parameters for the ODE model.
- `initstates`: Initial state vector `u₀`.
- `timeidx`: Time points at which to save the solution.
- `config`: A dictionary-like object containing solver configurations, such as:
  - `:device`: The target device for the output array (e.g., `cpu` or `gpu`).
  - `:solve_alg`: The ODE solver algorithm (defaults to `Tsit5()`).
  - `:sense_alg`: The sensitivity analysis algorithm for automatic differentiation.
  - `:solve_cb`: Callbacks to be applied during solving.

# Returns
- An `Array` containing the solution at the specified time points, moved to the target device.
"""
function HydroModels.hydrosolve(::Val{HydroModels.ODESolver}, du_func, pas, initstates, timeidx, config)
    device = get(config, :device, identity)
    solve_alg = get(config, :solve_alg, Tsit5())
    sense_alg = get(config, :sense_alg, GaussAdjoint())
    solve_cb = get(config, :solve_cb, nothing)

    function ode_func!(du, u, p, t)
        du[:] = du_func(u, p, t)
        return nothing
    end

    prob = ODEProblem{true}(ode_func!, initstates, (timeidx[1], timeidx[end]), pas)
    sol = solve(
        prob, solve_alg;
        saveat=timeidx,
        sensealg=sense_alg,
    )
    return Array(sol) |> device
end

"""
    hydrosolve(::Val{HydroModels.DiscreteSolver}, du_func, pas, initstates, timeidx, config)

Solve a discrete-time problem (map or recurrence relation).

This method dispatches on `DiscreteSolver` to set up and solve a `DiscreteProblem`.
It is suitable for models that evolve in discrete steps rather than continuously.

# Arguments
- `du_func`: The function defining the discrete update rule `uₙ₊₁ = f(uₙ, p, tₙ)`. (Note: Currently `ode_func!` is used internally).
- `pas`: Parameters for the model.
- `initstates`: Initial state vector `u₀`.
- `timeidx`: Time points at which to evaluate and save the solution.
- `config`: A dictionary-like object for solver configurations, such as:
  - `:device`: The target device for the output array.
  - `:solve_alg`: The discrete solver algorithm (defaults to `FunctionMap`).
  - `:sense_alg`: The sensitivity analysis algorithm (defaults to `ReverseDiffAdjoint`).
  - `:solve_cb`: Callbacks to be applied during solving.

# Returns
- An `Array` containing the solution at the specified time points, moved to the target device.
"""
function HydroModels.hydrosolve(::Val{HydroModels.DiscreteSolver}, du_func, pas, initstates, timeidx, config)
    device = get(config, :device, identity)
    solve_alg = get(config, :solve_alg, FunctionMap{true}())
    sense_alg = get(config, :sense_alg, ReverseDiffAdjoint())
    # solve_cb = get(config, :solve_cb, nothing)

    function ode_func!(du, u, p, t)
        du[:] = du_func(u, p, t)
        return nothing
    end

    prob = DiscreteProblem(ode_func!, initstates, (timeidx[1], timeidx[end]), pas)
    sol = solve(
        prob, solve_alg;
        # callback=solve_cb,
        saveat=timeidx,
        sensealg=sense_alg,
    )
    return Array(sol) |> device
end

"""
    SciMLBase.ODEProblem(bucket::HydroModels.HydroBucket, input::AbstractArray{T,2}; kwargs...) where T

Construct an `ODEProblem` from a `HydroBucket` and input data.

# Arguments
- `bucket`: A `HydroBucket` instance containing the model's ODE and flux functions.
- `input`: A 2D array of input data (forcings), where rows are variables and columns are time steps.

# Keyword Arguments
- `params`: Parameter vector for the model.
- `interpolator`: The interpolation method (defaults to `Val(ConstantInterpolation)`).
- `timeidx`: Time points corresponding to columns of `input` (defaults to `1:size(input, 2)`).
- `initstates`: Initial states (defaults to zero vector).

# Returns
- A `Tuple{ODEProblem, SavedValues}` containing the configured problem and callback for saved flux values.
"""
function SciMLBase.ODEProblem(bucket::HydroModels.HydroBucket, input::AbstractArray{T,2}; kwargs...) where T
    # TODO: Enzyme.jl support for DE.jl+DataInterpolations.jl is not yet mature

    params = get(kwargs, :params, nothing)
    isnothing(params) && error("params keyword argument is required")

    interp = get(kwargs, :interpolator, Val(HydroModels.ConstantInterpolation))
    timeidx = get(kwargs, :timeidx, collect(1:size(input, 2)))
    initstates = get(kwargs, :initstates, zeros(T, length(HydroModels.get_state_names(bucket))))
    itpfuncs = HydroModels.hydrointerp(interp, input, timeidx)

    function ode_func!(du, u, p, t)
        du[:] = bucket.ode_func(itpfuncs(t), u, p)
    end

    saved_values = SavedValues(T, Vector{Vector{T}})
    cb = SavingCallback(
        (u, t, integrator) -> bucket.flux_func(itpfuncs(t), u, integrator.p),
        saved_values
    )
    prob = ODEProblem{true}(ode_func!, initstates, (timeidx[1], timeidx[end]), params; callback=cb)
    return prob, saved_values
end

end