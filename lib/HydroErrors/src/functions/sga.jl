@doc """
    Compute the Spectral Gradient Angle (SGA).

    .. image:: /pictures/SGA.png

    **Range:** -π/2 ≤ SID < π/2, closer to 0 is better.

    **Notes:** The spectral gradient angle measures the angle between the two vectors in
    hyperspace. It indicates how well the shape of the two series match – not magnitude.
    SG is the gradient of the simulated or observed time series.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The Spectral Gradient Angle.

    References
    ----------
    - Robila, S.A., Gershman, A., 2005. Spectral matching accuracy in processing hyperspectral
      data, Signals, Circuits and Systems, 2005. ISSCS 2005. International Symposium on. IEEE,
      pp. 163-166.
 """
function sga(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    sgx = observed_array[1:end] .- observed_array[:observed_array.size .- 1]
    sgy = simulated_array[1:end] .- simulated_array[:simulated_array.size .- 1]
    a = dot(sgx, sgy)
    b = linalg.norm(sgx) .* linalg.norm(sgy)
    arccos(a ./ b)
end