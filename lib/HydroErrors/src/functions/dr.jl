@doc """
    Compute the the refined index of agreement (dr).

    .. image:: /pictures/dr.png

    **Range:** -1 ≤ dr < 1, does not indicate bias, larger is better.

    **Notes:** Reformulation of Willmott’s index of agreement. This metric was created to address
    issues in the index of agreement and the Nash-Sutcliffe efficiency metric. Meant to be a
    flexible metric for use in climatology.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The refined index of agreement.

    References
    ----------
    - Willmott, C.J., Robeson, S.M., Matsuura, K., 2012. A refined index of model performance.
      International Journal of Climatology 32(13) 2088-2094.
 """
function dr(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    a = sum(abs(simulated_array .- observed_array))
    b = 2 .* sum(abs(observed_array .- observed_array.mean()))
    if a <= b
         1 .- (a ./ b)
    else
         (b ./ a) .- 1
    end
end