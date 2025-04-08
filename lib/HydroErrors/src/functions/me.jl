@doc """
Compute the mean error of the simulated and observed data.

    .. image:: /pictures/ME.png

    **Range:** -inf < MAE < inf, data units, closer to zero is better, indicates bias.

    **Notes:** The mean error (ME) measures the difference between the simulated data and the
    observed data. For the mean error, a smaller number indicates a better fit to the original
    data. Note that if the error is in the form of random noise, the mean error will be very small,
    which can skew the accuracy of this metric. ME is cumulative and will be small even if there
    are large positive and negative errors that balance.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The mean error value.

    Examples
    --------
    Note that in this example the random noise cancels, leaving a very small ME.
 
    References
    ----------
    - Fisher, R.A., 1920. A Mathematical Examination of the Methods of Determining the Accuracy of
      an Observation by the Mean Error, and by the Mean Square Error. Monthly Notices of the Royal
      Astronomical Society 80 758 - 770.
 """
function me(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    mean(simulated_array .- observed_array)
end