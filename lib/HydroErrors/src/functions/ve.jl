@doc """
    Compute the Volumetric Efficiency (VE).

    .. image:: /pictures/VE.png

    **Range:** 0 ≤ VE < 1 smaller is better, does not indicate bias.

    **Notes:** Represents the error as a percentage of flow.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The Volumetric Efficiency value.

    References
    ----------
    - Criss, R.E., Winston, W.E., 2008. Do Nash values have value? Discussion and alternate
      proposals. Hydrological Processes 22(14) 2723.
 """
function ve(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    a = sum(abs(simulated_array .- observed_array))
    b = sum(observed_array)
    1 .- (a ./ b)
end