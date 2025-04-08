@doc """
Compute the Spectral Angle (SA).

    .. image:: /pictures/SA.png

    **Range:** -π/2 ≤ SA < π/2, closer to 0 is better.

    **Notes:** The spectral angle metric measures the angle between the two vectors in hyperspace.
    It indicates how well the shape of the two series match – not magnitude.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The Spectral Angle value.

    References
    ----------
    - Robila, S.A., Gershman, A., 2005. Spectral matching accuracy in processing hyperspectral
      data, Signals, Circuits and Systems, 2005. ISSCS 2005. International Symposium on. IEEE,
      pp. 163-166.
 """
function sa(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    a = dot(simulated_array, observed_array)
    b = linalg.norm(simulated_array) .* linalg.norm(observed_array)
    arccos(a ./ b)
end