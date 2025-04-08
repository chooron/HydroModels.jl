@doc """
    Compute the the mean absolute percentage error (MAPE).

    .. image:: /pictures/MAPE.png

    **Range:** 0% ≤ MAPE ≤ inf. 0% indicates perfect correlation, a larger error indicates a
    larger percent error in the data.

    **Notes:**

    Parameters
    ----------
    simulated_array: one dimensional ndarray
    An array of simulated data from the time series.

    observed_array: one dimensional ndarray
    An array of observed data from the time series.
        
    Returns
    -------
    float
        The mean absolute percentage error.
 
    References
    ----------
 """
function mape(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    a = simulated_array .- observed_array
    b = abs(a ./ observed_array)
    c = 100 ./ simulated_array.size
    c .* sum(b)
end