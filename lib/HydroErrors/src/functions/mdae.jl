@doc """
    Compute the median absolute error (MdAE) between the simulated and observed data.

    .. image:: /pictures/MdAE.png

    **Range** 0 ≤ MdAE < inf, closer to zero is better.

    **Notes** Random errors (noise) do not cancel. It is the same as the mean absolute error (MAE), only it takes the
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
        The median absolute error value.
 """
function mdae(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    median(abs(simulated_array .- observed_array))
end