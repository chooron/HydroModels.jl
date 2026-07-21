module HydroModelsDiffEqCallbacksExt

using HydroModels
using SciMLBase
using DiffEqCallbacks: SavedValues, SavingCallback

"""
    SciMLBase.ODEProblem(bucket::HydroModels.HydroBucket, input::AbstractArray{T,2}; kwargs...)

Construct a bucket ODE and a `SavingCallback` for flux values.  The callback
extension is intentionally gated only by `DiffEqCallbacks`.
"""
function SciMLBase.ODEProblem(
    bucket::HydroModels.HydroBucket,
    input::AbstractArray{T,2};
    kwargs...,
) where {T}
    params = get(kwargs, :params, nothing)
    isnothing(params) && error("params keyword argument is required")

    interp = get(kwargs, :interpolator, HydroModels.ConstantInterpolation)
    timeidx = get(kwargs, :timeidx, collect(1:size(input, 2)))
    initstates = get(
        kwargs,
        :initstates,
        zeros(T, length(HydroModels.get_state_names(bucket))),
    )
    itpfuncs = HydroModels.hydrointerp(interp, input, timeidx)

    function ode_func!(du, u, p, t)
        du .= bucket.ode_func(itpfuncs(t), u, p)
        return nothing
    end

    # HydroBucket flux functions return a freshly allocated Vector{T} for the
    # single-node path.  The copy in the save function also protects the
    # callback from user-defined bucket functions returning a view into `u`.
    # The solver promotes integer tspan values to the state/input element
    # type; keep SavedValues' time type aligned with that callback time.
    saved_values = SavedValues(T, Vector{T})
    save_func = (u, t, integrator) ->
        copy(bucket.flux_func(itpfuncs(t), u, integrator.p))
    cb = SavingCallback(save_func, saved_values)

    prob = ODEProblem{true}(
        ode_func!,
        initstates,
        (first(timeidx), last(timeidx)),
        params;
        callback=cb,
    )

    return prob, saved_values
end

end
