@doc """
Compute the mean squared log error of the simulated and observed data.

    .. image:: /pictures/MSLE.png

    **Range:** 0 ≤ MSLE < inf, data units squared, smaller is better.

    **Notes** Same as the mean squared error (MSE) only use log ratios as the error term. Limits
    the impact of outliers, more evenly weights high and low values.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The mean squared log error value.

    References
    ----------
    - Törnqvist, Leo, Pentti Vartia, and Yrjö O. Vartia. “How Should Relative Changes Be Measured?”
      The American Statistician 39, no. 1 (1985): 43–46.
 """
function msle(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    sim_log = log1p(simulated_array)
    obs_log = log1p(observed_array)
    mean((sim_log .- obs_log) .^ 2)
end