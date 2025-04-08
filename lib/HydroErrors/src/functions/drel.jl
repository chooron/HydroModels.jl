@doc """
    Compute the the relative index of agreement (drel).

    .. image:: /pictures/drel.png

    **Range:** 0 ≤ drel < 1, does not indicate bias, larger is better.

    **Notes:** Instead of absolute differences, this metric uses relative differences.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The relative index of agreement.

    References
    ----------
    - Krause, P., Boyle, D., Bäse, F., 2005. Comparison of different efficiency criteria for
      hydrological model assessment. Advances in geosciences 5 89-97.
 """
function drel(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    a = ((simulated_array .- observed_array) ./ observed_array) .^ 2
    b = abs(simulated_array .- mean(observed_array))
    c = abs(observed_array .- mean(observed_array))
    e = ((b .+ c) ./ mean(observed_array)) .^ 2
    1 .- (sum(a) ./ sum(e))
end