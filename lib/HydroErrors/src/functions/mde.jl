@doc """
 Compute the median error (MdE) between the simulated and observed data.

    .. image:: /pictures/MdE.png

    **Range** -inf < MdE < inf, closer to zero is better.

    **Notes** This metric indicates bias. It is similar to the mean error (ME), only it takes the
    median rather than the mean. Median measures reduces the impact of outliers.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.
        
    Returns
    -------
    float
        The median error value.
 """
function mde(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    median(simulated_array .- observed_array)
end