@doc """
 Compute the Euclidean distance between predicted and observed values in vector space.

    .. image:: /pictures/ED.png

    **Range** 0 ≤ ED < inf, smaller is better.
    **Notes** Also sometimes referred to as the L2-norm.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The euclidean distance error value.

    References
    ----------
    - Kennard, M. J., Mackay, S. J., Pusey, B. J., Olden, J. D., & Marsh, N. (2010). Quantifying
      uncertainty in estimation of hydrologic metrics for ecohydrological studies. River Research
      and Applications, 26(2), 137-156.
 """
function ed(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    linalg.norm(observed_array .- simulated_array)
end