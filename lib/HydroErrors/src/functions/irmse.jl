@doc """
    Compute the inertial root mean square error (IRMSE) between the simulated and observed data.

    .. image:: /pictures/IRMSE.png

    **Range:** 0 ≤ IRMSE < inf, lower is better.

    **Notes:** This metric is the RMSE devided by by the standard deviation of the gradient of the
    observed timeseries data. This metric is meant to be help understand the ability of the model
    to predict changes in observation.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The inertial root mean square error.

    References
    ----------
    - Daga, M., Deo, M.C., 2009. Alternative data-driven methods to estimate wind from waves by
      inverse modeling. Natural Hazards 49(2) 293-310.
 """
function irmse(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    # Getting the gradient of the observed data
    obs_len = observed_array.size
    obs_grad = observed_array[1:obs_len] .- observed_array[0:obs_len .- 1]
    
    # Standard deviation of the gradient
    obs_grad_std = std(obs_grad, ddof=1)
    
    # Divide RMSE by the standard deviation of the gradient of the observed data
    rmse_value = sqrt(mean((simulated_array .- observed_array) .^ 2))
    rmse_value ./ obs_grad_std
end