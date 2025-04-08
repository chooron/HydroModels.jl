@doc """
    Compute the Legate-McCabe Efficiency Index.

    .. image:: /pictures/E1p.png

    **Range:** 0 ≤ E1' < 1, does not indicate bias, larger is better.

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
        The Legate-McCabe Efficiency index value.

    References
    ----------
    - Legates, D.R., McCabe Jr, G.J., 1999. Evaluating the use of “goodness‐of‐fit” Measures in
      hydrologic and hydroclimatic model validation. Water Resources Research 35(1) 233-241.
      Lehmann, E.L., Casella, G., 1998. Springer Texts in Statistics. Springer-Verlag, New York.
 """
function lm_index(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    mean_obs = mean(observed_array)

    if obs_bar_p ! None
        a = abs(simulated_array .- observed_array)
        b = abs(observed_array .- obs_bar_p)
        return 1 .- (sum(a) ./ sum(b))
    else
        a = abs(simulated_array .- observed_array)
        b = abs(observed_array .- mean_obs)
        return 1 .- (sum(a) ./ sum(b))
    end
end