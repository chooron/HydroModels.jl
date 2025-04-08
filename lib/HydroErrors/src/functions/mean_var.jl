@doc """
    Compute the mean variance.

    .. image:: /pictures/MV.png

    **Range:**

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
        The mean variance.
    
    References
    ----------
   
 """
function mean_var(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    var(log1p(observed_array) .- og1p(simulated_array))
end