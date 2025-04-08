@doc """
    Compute the Spectral Information Divergence (SID).

    .. image:: /pictures/SID.png

    **Range:** -π/2 ≤ SID < π/2, closer to 0 is better.

    **Notes:** The spectral information divergence measures the angle between the two vectors in
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
        The Spectral information divergence value.

    References
    ----------
    - Robila, S.A., Gershman, A., 2005. Spectral matching accuracy in processing hyperspectral
      data, Signals, Circuits and Systems, 2005. ISSCS 2005. International Symposium on. IEEE,
      pp. 163-166.
 """
function sid(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    first = (observed_array ./ mean(observed_array)) .- (simulated_array ./ mean(simulated_array))
    second1 = log10(observed_array) .- log10(mean(observed_array))
    second2 = log10(simulated_array) .- log10(mean(simulated_array))
    dot(first, second1 ,- second2)
end