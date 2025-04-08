@doc"""
    Compute the root mean square log error between the simulated and observed data.

    .. image:: /pictures/RMSLE.png

    **Range:** 0 ≤ RMSLE < inf. Smaller is better, and it does not indicate bias.

    **Notes:** Random errors do not cancel while using this metric. This metric limits the
    impact of outliers by more evenly weighting high and low values. To calculate the log values,
    each value in the observed and simulated array is increased by one unit in order to avoid
    run-time errors and nan values (function np.log1p).

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.
    
    Returns
    -------
    float
        The root mean square log error value.
    
    References
    ----------
    - Tornquist, L., Vartia, P., Vartia, Y.O., 1985. How Should Relative Changes be Measured?
      The American Statistician 43-46.
    - Willmott, C.J., Matsuura, K., 2005. Advantages of the mean absolute error (MAE) over the
      root mean square error (RMSE) in assessing average model performance.
      Climate Research 30(1) 79-82.
"""
function rmsle(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    sqrt(mean(power(log1p(simulated_array) .- log1p(observed_array), 2)))
end