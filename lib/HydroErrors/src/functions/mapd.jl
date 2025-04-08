@doc """
    Compute the the mean absolute percentage deviation (MAPD).

    .. image:: /pictures/MAPD.png

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
        The mean absolute percentage deviation.
 
    References
    ----------
 """
function mapd(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    a = sum(abs(simulated_array .- observed_array))
    b = sum(abs(observed_array))
    a ./ b
end