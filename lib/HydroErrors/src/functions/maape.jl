@doc """
    Compute the the Mean Arctangent Absolute Percentage Error (MAAPE).

    .. image:: /pictures/MAAPE.png

    **Range:** 0 ≤ MAAPE < π/2, does not indicate bias, smaller is better.

    **Notes:** Represents the mean absolute error as a percentage of the observed values. Handles
    0s in the observed data. This metric is not as biased as MAPE by under-over predictions.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The mean arctangent absolute percentage error.

    References
    ----------
    - Kim, S., Kim, H., 2016. A new metric of absolute percentage error for intermittent demand
      forecasts. International Journal of Forecasting 32(3) 669-679.
 """
function maappe(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    a = simulated_array .- observed_array
    b = abs(a ./ observed_array)
    mean(arctan(b))
end