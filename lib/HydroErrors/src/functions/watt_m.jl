@doc """
    Compute Watterson's M (M).

    .. image:: /pictures/M.png

    **Range:** -1 ≤ M < 1, does not indicate bias, larger is better.

    **Notes:**

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        Watterson's M value.

    References
    ----------
    - Watterson, I.G., 1996. Non‐dimensional measures of climate model performance. International
      Journal of Climatology 16(4) 379-391.
 """
function watt_m(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    a = 2 ./ pi
    b = mean((simulated_array .- observed_array) .^ 2)  # MSE
    c = std(observed_array, ddof=1) .^ 2 + std(simulated_array, ddof=1) .^ 2
    e = (mean(simulated_array) .- mean(observed_array)) .^ 2
    f = c .+ e
    a .* arcsin(1 .- (b ./ f))
end