@doc"""
    Compute the root mean square error between the simulated and observed data.

    .. image:: /pictures/RMSE.png

    **Range** 0 ≤ RMSE < inf, smaller is better.

    **Notes:** The standard deviation of the residuals. A lower spread indicates that the points
    are better concentrated around the line of best fit (linear). Random errors do not cancel.
    This metric will highlights larger errors.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.
    
    Returns
    -------
    float
        The root mean square error value.
    
    References
    ----------
    - Willmott, C.J., Matsuura, K., 2005. Advantages of the mean absolute error (MAE) over the
      root mean square error (RMSE) in assessing average model performance.
      Climate Research 30(1) 79-82.
    - Hyndman, R.J., Koehler, A.B., 2006. Another look at measures of forecast accuracy.
      International Journal of Forecasting 22(4) 679-688.
"""
function rmse(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    sqrt(mean((simulated_array .- observed_array) .^ 2))
end