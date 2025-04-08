@doc """
    Compute the the index of agreement (d).

    .. image:: /pictures/d.png

    **Range:** 0 ≤ d < 1, does not indicate bias, larger is better.

    **Notes:** This metric is a modified approach to the Nash-Sutcliffe Efficiency metric.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The index of agreement (1).

    References
    ----------
    - Legates, D.R., McCabe Jr, G.J., 1999. Evaluating the use of “goodness‐of‐fit” Measures in
      hydrologic and hydroclimatic model validation. Water Resources Research 35(1) 233-241.
    - Willmott, C.J., Robeson, S.M., Matsuura, K., 2012. A refined index of model performance.
      International Journal of Climatology 32(13) 2088-2094.
 """
function d(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    a = (observed_array - simulated_array) ^ 2
    b = abs(simulated_array - mean(observed_array))
    c = abs(observed_array - mean(observed_array))
    1 .- (sum(a) ./ sum((b + c) .^ 2))
end