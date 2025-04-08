@doc """
Compute the mean absolute log error of the simulated and observed data.

    .. image:: /pictures/MALE.png

    **Range:** 0 ≤ MALE < inf, data units squared, smaller is better.

    **Notes** Same as MAE only use log ratios as the error term. Limits the impact of outliers,
    more evenly weights high and low flows.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.
        
    Returns
    -------
    float
        The mean absolute log error value.
 
    References
    ----------
    - Törnqvist, Leo, Pentti Vartia, and Yrjö O. Vartia. “How Should Relative Changes Be Measured?”
      The American Statistician 39, no. 1 (1985): 43–46.
 """
function male(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    sim_log = log1p(simulated_array)
    obs_log = log1p(observed_array)
    mean(abs(sim_log .- obs_log))
end