@doc """
    Compute the pearson correlation coefficient.

    .. image:: /pictures/R_pearson.png

    **Range:** -1 ≤ R (Pearson) ≤ 1. 1 indicates perfect postive correlation, 0 indicates
    complete randomness, -1 indicate perfect negative correlation.

    **Notes:** The pearson r coefficient measures linear correlation. It is sensitive to outliers.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The Pearson correlation coefficient.

    References
    ----------
    - Pearson, K. (1895). Note on regression and inheritance in the case of two parents.
      Proceedings of the Royal Society of London, 58, 240-242.
 """
function person_r(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    sim_mean = mean(simulated_array)
    obs_mean = mean(observed_array)

    top = sum((observed_array .- obs_mean) .* (simulated_array .- sim_mean))
    bot1 = sqrt(sum((observed_array .- obs_mean) .^ 2))
    bot2 = sqrt(sum((simulated_array .- sim_mean) .^ 2))

    top ./ (bot1 .* bot2)
end

@doc """
    Compute the spearman rank correlation coefficient.

    .. image:: /pictures/R_spearman.png

    **Range:** -1 ≤ R (Pearson) ≤ 1. 1 indicates perfect postive correlation, 0 indicates
    complete randomness, -1 indicate perfect negative correlation.

    **Notes:** The spearman r coefficient measures the monotonic relation between simulated and
    observed data. Because it uses a nonparametric measure of rank correlation, it is less sensitive
    to outliers compared to the Pearson correlation coefficient.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The Spearman rank correlation coefficient.

    References
    ----------
    - Spearman C (1904). "The proof and measurement of association between two things". American
      Journal of Psychology. 15: 72–101. doi:10.2307/1412159
 """
function spearman_r(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    rank_sim = rankdata(simulated_array)
    rank_obs = rankdata(observed_array)

    mean_rank_sim = mean(rank_sim)
    mean_rank_obs = mean(rank_obs)

    top = mean((rank_obs .- mean_rank_obs) .* (rank_sim .- mean_rank_sim))
    bot = sqrt(mean((rank_obs .- mean_rank_obs) .^ 2) .* mean((rank_sim - mean_rank_sim) .^ 2))

    top ./ bot
end

@doc """
    Compute the the Coefficient of Determination (r2).

    .. image:: /pictures/r2.png

    **Range:** 0 ≤ r2 ≤ 1. 1 indicates perfect correlation, 0 indicates complete randomness.

    **Notes:** The Coefficient of Determination measures the linear relation between simulated and
    observed data. Because it is the pearson correlation coefficient squared, it is more heavily
    affected by outliers than the pearson correlation coefficient.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The coefficient of determination (R^2).

    References
    ----------
 """
function r_squared(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    a = observed_array .- mean(observed_array)
    b = simulated_array .- mean(simulated_array)
    (sum(a .* b)) .^ 2 / (sum(a .^ 2) .* sum(b .^ 2))
end

@doc """
    Compute Mielke-Berry R value (MB R).

    .. image:: /pictures/MB_R.png

    **Range:** 0 ≤ MB R < 1, does not indicate bias, larger is better.

    **Notes:** Compares prediction to probability it arose by chance.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The Mielke-Berry R value.
    Notes
    -----
    If a more optimized version is desired, the `numba package <http://numba.pydata.org/doc.html>`_ can be implemented
    for a much more optimized performance when computing this metric. An example is given below.

    References
    ----------
    - Berry, K.J., Mielke, P.W., 1988. A Generalization of Cohen's Kappa Agreement Measure to
      Interval Measurement and Multiple Raters. Educational and Psychological Measurement 48(4)
      921-933.
    - Mielke, P.W., Berry, K.J., 2007. Permutation methods: a distance function approach.
      Springer Science & Business Media.
 """
function mb_r(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    n = simulated_array.size
    tot = 0.0
    for i in range(n)
        tot = tot .+ sum(abs(simulated_array .- observed_array[i]))
    mae_val = sum(abs(simulated_array .- observed_array)) ./ n
    mb = 1 .- ((n .^ 2) .* mae_val ./ tot)
    mb
    end
end