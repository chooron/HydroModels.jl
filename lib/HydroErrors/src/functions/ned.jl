@doc """
    Compute the normalized Euclidian distance between the simulated and observed data in vector
    space.

    .. image:: /pictures/NED.png

    **Range** 0 ≤ NED < inf, smaller is better.

    **Notes** Also sometimes referred to as the squared L2-norm.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The normalized euclidean distance value.

    References
    ----------
    - Kennard, M. J., Mackay, S. J., Pusey, B. J., Olden, J. D., & Marsh, N. (2010). Quantifying
      uncertainty in estimation of hydrologic metrics for ecohydrological studies. River Research
      and Applications, 26(2), 137-156.
 """
function ned(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    a = observed_array ./ mean(observed_array)
    b = simulated_array ./ mean(simulated_array)
    linalg.norm(a .- b)
end