@doc """
    Compute the the index of agreement (d1).

    .. image:: /pictures/d1.png

    **Range:** 0 ≤ d < 1, does not indicate bias, larger is better.

    **Notes:** This metric is a modified approach to the Nash-Sutcliffe Efficiency metric. Compared
    to the other index of agreement (d) it has a reduced impact of outliers.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The index of agreement (d1).

    References
    ----------
    - Willmott, C.J., Robeson, S.M., Matsuura, K., 2012. A refined index of model performance.
      International Journal of Climatology 32(13) 2088-2094.
 """
function d1(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    obs_mean = mean(observed_array)
    a = sum(abs(simulated_array .- observed_array))
    b = abs(simulated_array .- obs_mean)
    c = abs(observed_array .- obs_mean)
    1 .- sum(a) ./ sum(b .+ c)
end

@doc """
    Compute the Legate-McCabe Index of Agreement.

    .. image:: /pictures/D1p.png

    **Range:** 0 ≤ d1' < 1, does not indicate bias, larger is better.

    **Notes:** The obs_bar_p argument represents a seasonal or other selected average.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    obs_bar_p: float
        Seasonal or other selected average. If None, the mean of the observed array will be used.

    Returns
    -------
    float
        The Legate-McCabe Efficiency index of agreement.

    References
    ----------
    - Legates, D.R., McCabe Jr, G.J., 1999. Evaluating the use of “goodness‐of‐fit” Measures in
      hydrologic and hydroclimatic model validation. Water Resources Research 35(1) 233-241.
      Lehmann, E.L., Casella, G., 1998. Springer Texts in Statistics. Springer-Verlag, New York.
 """
function d1_p(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    if obs_bar_p ! None
        a = abs(observed_array .- simulated_array)
        b = abs(simulated_array .- obs_bar_p) .+ abs(observed_array .- obs_bar_p)
        1 .- (sum(a) ./ sum(b))
    else
        mean_obs = mean(observed_array)
        a = abs(observed_array .- simulated_array)
        b = abs(simulated_array .- mean_obs) .+ abs(observed_array .- mean_obs)
        1 .- (sum(a) ./ sum(b))
    end
end