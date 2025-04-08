@doc """
    Compute the median squared error (MdSE) between the simulated and observed data.

    .. image:: /pictures/MdSE.png

    **Range** 0 ≤ MdSE < inf, closer to zero is better.

    **Notes** Random errors (noise) do not cancel. It is the same as the mean squared error (MSE), only it takes the
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
        The median squared error value.
 """
function mdse(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    median((simulated_array .- observed_array) .^ 2)
end