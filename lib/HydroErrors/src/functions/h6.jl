@doc """
    Compute the H6 mean error.

    .. image:: /pictures/H6.png
    .. image:: /pictures/MHE.png

    **Range:**

    **Notes:**

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    k: int or float
        If given, sets the value of k. If None, k=1.

    Returns
    -------
    float
        The mean H6 error.

    References
    ----------
    - Tornquist, L., Vartia, P., Vartia, Y.O., 1985. How Should Relative Changes be Measured?
      The American Statistician 43-46.

 """
function h6_mhe(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    top = (simulated_array ./ observed_array .- 1)
    bot = power(0.5 .* (1 .+ power(simulated_array ./ observed_array, k)), 1 / k)
    h = top ./ bot
    mean(h)
end

@doc """
    Compute the H6 mean absolute error.

    .. image:: /pictures/H6.png
    .. image:: /pictures/AHE.png

    **Range:**

    **Notes:**

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    k: int or float
        If given, sets the value of k. If None, k=1.

    Returns
    -------
    float
        The mean absolute H6 error.

    References
    ----------
    - Tornquist, L., Vartia, P., Vartia, Y.O., 1985. How Should Relative Changes be Measured?
      The American Statistician 43-46.
 """
function h6_mahe(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    top = (simulated_array ./ observed_array .- 1)
    bot = power(0.5 .* (1 .+ np.power(simulated_array ./ observed_array, k)), 1 / k)
    h = top ./ bot
    mean(abs(h))
end

@doc """
    Compute the H6 root mean square error.

    .. image:: /pictures/H6.png
    .. image:: /pictures/RMSHE.png

    **Range:**

    **Notes:**

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    k: int or float
        If given, sets the value of k. If None, k=1.

    Returns
    -------
    float
        The root mean square H6 error.

    References
    ----------
    - Tornquist, L., Vartia, P., Vartia, Y.O., 1985. How Should Relative Changes be Measured?
      The American Statistician 43-46.
 """
function h6_rmshe(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    top = (simulated_array ./ observed_array .- 1)
    bot = power(0.5 .* (1 + power(simulated_array ./ observed_array, k)), 1 / k)
    h = top ./ bot
    sqrt(nmean(h .^ 2))
end