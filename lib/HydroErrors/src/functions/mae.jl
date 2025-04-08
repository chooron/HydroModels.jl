@doc """
Compute the mean absolute error of the simulated and observed data.

    .. image:: /pictures/MAE.png

    **Range:** 0 ≤ MAE < inf, data units, smaller is better.

    **Notes:** The ME measures the absolute difference between the simulated data and the observed
    data. For the mean abolute error, a smaller number indicates a better fit to the original data.
    Also note that random errors do not cancel. Also referred to as an L1-norm.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The mean absolute error value.

    References
    ----------
    - Willmott, Cort J., and Kenji Matsuura. “Advantages of the Mean Absolute Error (MAE) over the
      Root Mean Square Error (RMSE) in Assessing Average Model Performance.” Climate Research 30,
      no. 1 (2005): 79–82.
    - Willmott, Cort J., and Kenji Matsuura. “On the Use of Dimensioned Measures of Error to
      Evaluate the Performance of Spatial Interpolators.” International Journal of Geographical
      Information Science 20, no. 1 (2006): 89–102.
 """
function mae(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    mean(absolute(simulated_array .- observed_array))
end