@doc """
    Compute the geometric mean difference.

    .. image:: /pictures/GMD.png

    **Range:**

    **Notes:** For the difference of geometric means, the geometric mean is computed for each of
    two samples then their difference is taken.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The geometric mean difference value.

    References
    ----------
 """
function g_mean_diff(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    sim_log = log1p(simulated_array)
    obs_log = log1p(observed_array)
    exp(gmean(sim_log) .- gmean(obs_log))
end