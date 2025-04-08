@doc """
    Compute the Spectral Correlation (SC).

    .. image:: /pictures/SC.png

    **Range:** -π/2 ≤ SA < π/2, closer to 0 is better.

    **Notes:** The spectral correlation metric measures the angle between the two vectors in
    hyperspace. It indicates how well the shape of the two series match – not magnitude.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The Spectral Correlation value.

    References
    ----------
    - Robila, S.A., Gershman, A., 2005. Spectral matching accuracy in processing hyperspectral
      data, Signals, Circuits and Systems, 2005. ISSCS 2005. International Symposium on. IEEE,
      pp. 163-166.
 """
function sc(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    a = dot(observed_array .- mean(observed_array), simulated_array .- mean(simulated_array))
    b = linalg.norm(observed_array .- mean(observed_array))
    c = linalg.norm(simulated_array .- mean(simulated_array))
    e = b .* c
    arccos(a ./ e)
end