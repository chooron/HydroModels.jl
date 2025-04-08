@doc """
    Compute the the modified index of agreement (dmod).

    .. image:: /pictures/dmod.png

    **Range:** 0 ≤ dmod < 1, does not indicate bias, larger is better.

    **Notes:** When j=1, this metric is the same as d1. As j becomes larger, outliers have a larger
    impact on the value.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    j: int or float
        Optional input indicating the j values desired. A higher j places more emphasis on
        outliers. j is 1 by default.

    Returns
    -------
    float
        The modified index of agreement.

    References
    ----------

    - Krause, P., Boyle, D., Bäse, F., 2005. Comparison of different efficiency criteria for
      hydrological model assessment. Advances in geosciences 5 89-97.
 """
function dmod(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    a = (abs(simulated_array .- observed_array)) .^ j
    b = abs(simulated_array .- mean(observed_array))
    c = abs(observed_array .- mean(observed_array))
    e = (b .+ c) .^ j
    1 .- (sum(a) ./ sum(e))
end