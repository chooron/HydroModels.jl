@doc """
    Compute the Nash-Sutcliffe Efficiency.

    .. image:: /pictures/NSE.png

    **Range:** -inf < NSE < 1, does not indicate bias, larger is better.

    **Notes:** The Nash-Sutcliffe efficiency metric compares prediction values to naive predictions
    (i.e. average value). One major flaw of this metric is that it punishes a higher variance in
    the observed values (denominator). This metric is analogous to the mean absolute error skill
    score (MAESS) using the mean flow as a benchmark.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The Nash-Sutcliffe Efficiency value.

    References
    ----------
    - Krause, P., Boyle, D., Bäse, F., 2005. Comparison of different efficiency criteria for
      hydrological model assessment. Advances in geosciences 5 89-97.
    - McCuen, R.H., Knight, Z., Cutter, A.G., 2006. Evaluation of the Nash-Sutcliffe Efficiency
      Index. Journal of Hydraulic Engineering.
    - Nash, J.E., Sutcliffe, J.V., 1970. River flow forecasting through conceptual models part
      I — A discussion of principles. Journal of Hydrology 282-290.
    - Willmott, C.J., Robeson, S.M., Matsuura, K., 2012. A refined index of model performance.
      International Journal of Climatology 32(13) 2088-2094.
 """
function nse(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    a = @. abs(simulated_array - observed_array) ^ 2
    b = @. abs(observed_array - mean(observed_array)) ^ 2
    1 .- (sum(a) ./ sum(b))
end

@doc """
    Compute the modified Nash-Sutcliffe efficiency (NSE mod).

    .. image:: /pictures/NSEmod.png

    **Range:** -inf < NSE (mod) < 1, does not indicate bias, larger is better.

    **Notes:** The modified Nash-Sutcliffe efficiency metric gives less weight to outliers if j=1,
    or more weight to outliers if j is higher. Generally, j=1.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    j: int or float
        If given, sets the value of j to the input. j is 1 by default. A higher j gives more
        emphasis to outliers

    Returns
    -------
    float
        The modified Nash-Sutcliffe efficiency value.

    References
    ----------
    - Krause, P., Boyle, D., Bäse, F., 2005. Comparison of different efficiency criteria for
      hydrological model assessment. Advances in geosciences 5 89-97.
 """
function nse_mod(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    a = (abs(simulated_array .- observed_array)) .^ j
    b = (abs(observed_array .- mean(observed_array))) .^ j
    1 .- (sum(a) ./ sum(b))
end

@doc """
    Compute the relative Nash-Sutcliffe efficiency (NSE rel).

    .. image:: /pictures/NSErel.png

    **Range:** -inf < NSE (rel) < 1, does not indicate bias, larger is better.

    **Notes:** The modified Nash-Sutcliffe efficiency metric gives less weight to outliers if j=1,
    or more weight to outliers if j is higher. Generally, j=1.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The relative Nash-Sutcliffe efficiency value.

    References
    ----------
    - Krause, P., Boyle, D., Bäse, F., 2005. Comparison of different efficiency criteria for
      hydrological model assessment. Advances in geosciences 5 89-97.
 """
function nse_rel(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    a = (abs((simulated_array .- observed_array) ./ observed_array)) .^ 2
    b = (abs((observed_array .- mean(observed_array)) ./ mean(observed_array))) .^ 2
    1 .- (sum(a) ./ sum(b))
end