@doc """
Compute the mean squared error of the simulated and observed data.

    .. image:: /pictures/MSE.png

    **Range:** 0 ≤ MSE < inf, data units squared, smaller is better.

    **Notes:** Random errors do not cancel, highlights larger errors, also referred to as a
    squared L2-norm.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The mean squared error value.

    References
    ----------
    - Wang, Zhou, and Alan C. Bovik. “Mean Squared Error: Love It or Leave It? A New Look at Signal
      Fidelity Measures.” IEEE Signal Processing Magazine 26, no. 1 (2009): 98–117.
 """
function mse(simulated_array::AbstractVector{T}, observed_array::AbstractVector{T}; kwargs...)::T where {T}
    mean((simulated_array .- observed_array) .^ 2)
end