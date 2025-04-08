@doc"""
Compute the Kling-Gupta efficiency (2009).

    .. image:: /pictures/KGE_2009.png

    **Range:** -inf < KGE (2009) < 1, larger is better.

    **Notes:** Gupta et al. (2009) created this metric to demonstrate the relative importance of
    the three components of the NSE, which are correlation, bias and variability. This was done
    with hydrologic modeling as the context. This metric is meant to address issues with the NSE.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    s: tuple of length three
        Represents the scaling factors to be used for re-scaling the Pearson product-moment
        correlation coefficient (r), Alpha, and Beta, respectively.
    
    Returns
    -------
    float (tuple of float)
        The Kling-Gupta (2009) efficiency value, unless the return_all parameter is True.
    
        References
    ----------
    - Gupta, H. V., Kling, H., Yilmaz, K. K., & Martinez, G. F. (2009). Decomposition of the mean
      squared error and NSE performance criteria: Implications for improving hydrological modelling.
      Journal of Hydrology, 377(1-2), 80-91.
"""
function kge_2009(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    s = get(kwargs, :s, (1, 1, 1))

    # Means
    sim_mean = mean(simulated_array)
    obs_mean = mean(observed_array)

    # Standard Deviations
    sim_sigma = std(simulated_array, ddof=1)
    obs_sigma = std(observed_array, ddof=1)

    # Pearson R
    top_pr = sum((observed_array .- obs_mean) .* (simulated_array .- sim_mean))
    bot1_pr = sqrt(sum((observed_array .- obs_mean).^2))
    bot2_pr = sqrt(sum((simulated_array .- sim_mean).^2))
    pr = top_pr ./ (bot1_pr .* bot2_pr)

    # Ratio between mean of simulated and observed data
    if obs_mean != 0
        beta = sim_mean ./ obs_mean
    else
        beta = NaN
    end

    # Relative variability between simulated and observed values
    if obs_sigma != 0
        alpha = sim_sigma ./ obs_sigma
    else
        alpha = NaN
    end
    if !isnan(beta) & !isnan(alpha)
        kge = 1 .- sqrt(
            (s[1] .* (pr .- 1)) .^ 2 .+ (s[2] .* (alpha .- 1)) .^2 .+ (s[3] .* (beta .- 1)) .^2)
    else
        if obs_mean == 0
            @warn "Warning: The observed data mean is 0. Therefore, Beta is infinite and the KGE " *
                  "value cannot be computed."
        end
        if obs_sigma == 0
            @warn "Warning: The observed data standard deviation is 0. Therefore, Alpha is infinite " *
                  "and the KGE value cannot be computed."
        end
        kge = NaN
    end
    kge
end

@doc"""
    Compute the Kling-Gupta efficiency (2012).

    .. image:: /pictures/KGE_2012.png

    **Range:** -inf < KGE (2012) < 1, does not indicate bias, larger is better.

    **Notes:** The modified version of the KGE (2009). Kling proposed this version to avoid
    cross-correlation between bias and variability ratios.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    s: tuple of length three
        Represents the scaling factors to be used for re-scaling the Pearson product-moment
        correlation coefficient (r), gamma, and Beta, respectively.
    
    Returns
    -------
    float (tuple of float)
        The Kling-Gupta (2012) efficiency value, unless the return_all parameter is True.
    
    References
    ----------
    - Kling, H., Fuchs, M., & Paulin, M. (2012). Runoff conditions in the upper Danube basin under
      an ensemble of climate change scenarios. Journal of Hydrology, 424, 264-277.
"""
function kge_2012(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    s = get(kwargs, :s, (1, 1, 1))

    # Means
    sim_mean = mean(simulated_array)
    obs_mean = mean(observed_array)
    
    # Standard Deviations
    sim_sigma = std(simulated_array)
    obs_sigma = std(observed_array)
    
    # Pearson R
    top_pr = sum((observed_array .- obs_mean) .* (simulated_array .- sim_mean))
    bot1_pr = sqrt(sum((observed_array .- obs_mean) .^ 2))
    bot2_pr = sqrt(sum((simulated_array .- sim_mean) .^ 2))
    pr = top_pr ./ (bot1_pr .* bot2_pr)
    
    # Ratio between mean of simulated and observed data
    beta = sim_mean ./ obs_mean
    
    # CV is the coefficient of variation (standard deviation / mean)
    sim_cv = sim_sigma ./ sim_mean
    obs_cv = obs_sigma ./ obs_mean
    
    # Variability Ratio, or the ratio of simulated CV to observed CV
    gam = sim_cv ./ obs_cv
    
    if obs_mean != 0 & obs_sigma != 0 & sim_mean != 0
        kge = 1 .- qrt(
            (s[0] .* (pr .- 1)) .^ 2 .+ (s[1] .* (gam .- 1)) .^ 2 .+ (s[2] .* (beta .- 1)) .^ 2)
    else
        if obs_mean == 0
        @warn"Warning: The observed data mean is 0. Therefore, Beta is infinite and the KGE "*
            "value cannot be computed."
        end
        if obs_sigma == 0
        @warn"Warning: The observed data standard deviation is 0. Therefore, Gamma is infinite "*
            "and the KGE value cannot be computed."
        end
        if sim_mean == 0
        @warn"Warning: The simulated data mean is 0. Therefore, Gamma is infinite "
            "and the KGE value cannot be computed."
        end
        kge = NaN
    end
    kge
end