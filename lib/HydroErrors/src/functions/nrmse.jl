@doc """
   Compute the range normalized root mean square error between the simulated and observed data.

    .. image:: /pictures/NRMSE_Range.png

    **Range:** 0 ≤ NRMSE < inf.

    **Notes:** This metric is the RMSE normalized by the range of the observed time series (x).
    Normalizing allows comparison between data sets with different scales. The NRMSErange is the
    most sensitive to outliers of the three normalized rmse metrics.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The range normalized root mean square error value.

    References
    ----------
    - Pontius, R.G., Thontteh, O., Chen, H., 2008. Components of information for multiple
      resolution comparison between maps that share a real variable. Environmental and Ecological
      Statistics 15(2) 111-142.
 """
function nrmse_range(simulated_array::AbstractVector{T}, observed_array::AbstractVector{T}; kwargs...)::T where {T}
    rmse_value = sqrt(mean((simulated_array .- observed_array) .^ 2))
    obs_max = max(observed_array)
    obs_min = min(observed_array)
    rmse_value ./ (obs_max .- obs_min)
end

@doc """
    Compute the mean normalized root mean square error between the simulated and observed data.

    .. image:: /pictures/NRMSE_Mean.png

    **Range:** 0 ≤ NRMSE < inf.

    **Notes:** This metric is the RMSE normalized by the mean of the observed time series (x).
    Normalizing allows comparison between data sets with different scales.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The mean normalized root mean square error.

    References
    ----------
    - Pontius, R.G., Thontteh, O., Chen, H., 2008. Components of information for multiple
      resolution comparison between maps that share a real variable. Environmental and Ecological
      Statistics 15(2) 111-142.
 """
function nrmse_mean(simulated_array::AbstractVector{T}, observed_array::AbstractVector{T}; kwargs...)::T where {T} 
    rmse_value = sqrt(mean((simulated_array .- observed_array) .^ 2))
    obs_mean = mean(observed_array)
    rmse_value ./ obs_mean
end

@doc """
    Compute the IQR normalized root mean square error between the simulated and observed data.

    .. image:: /pictures/NRMSE_IQR.png

    **Range:** 0 ≤ NRMSE < inf.

    **Notes:** This metric is the RMSE normalized by the interquartile range of the observed time
    series (x). Normalizing allows comparison between data sets with different scales.
    The NRMSEquartile is the least sensitive to outliers of the three normalized rmse metrics.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The IQR normalized root mean square error.

    References
    ----------
    - Pontius, R.G., Thontteh, O., Chen, H., 2008. Components of information for multiple
      resolution comparison between maps that share a real variable. Environmental and Ecological
      Statistics 15(2) 111-142.
 """
function nrmse_iqr(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    rmse_value = sqrt(mean((simulated_array .- observed_array) .^ 2))
    q1 = percentile(observed_array, 25)
    q3 = percentile(observed_array, 75)
    iqr = q3 .- q1
    rmse_value ./ iqr
end