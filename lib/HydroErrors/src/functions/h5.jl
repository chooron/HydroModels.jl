@doc """
    Compute the H5 mean error.

    .. image:: /pictures/H5.png
    .. image:: /pictures/MHE.png

    **Range:**

    **Notes:**

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The mean H5 error.

    References
    ----------
    - Tornquist, L., Vartia, P., Vartia, Y.O., 1985. How Should Relative Changes be Measured?
      The American Statistician 43-46.

 """
function h5_mhe(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    top = (simulated_array .- observed_array)
    bot = reciprocal(0.5 .* (reciprocal(observed_array) .+ reciprocal(simulated_array)))
    h = top ./ bot
    mean(h)
end

@doc """
    Compute the H5 mean absolute error.

    .. image:: /pictures/H5.png
    .. image:: /pictures/AHE.png

    **Range:**

    **Notes:**

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The mean absolute H5 error.

    References
    ----------
    - Tornquist, L., Vartia, P., Vartia, Y.O., 1985. How Should Relative Changes be Measured?
      The American Statistician 43-46.
 """
function h5_mahe(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    top = (simulated_array - observed_array)
    bot = reciprocal(0.5 .* (reciprocal(observed_array) .+ reciprocal(simulated_array)))
    h = top ./ bot
    mean(abs(h))
end

@doc """
    Compute the H5 root mean square error.

    .. image:: /pictures/H5.png
    .. image:: /pictures/RMSHE.png

    **Range:**

    **Notes:**

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The root mean square H5 error.

    References
    ----------
    - Tornquist, L., Vartia, P., Vartia, Y.O., 1985. How Should Relative Changes be Measured?
      The American Statistician 43-46.

 """
function h5_rmshe(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    top = (simulated_array .- observed_array)
    bot = reciprocal(0.5 .* (reciprocal(observed_array) .+ reciprocal(simulated_array)))
    h = top ./ bot
    sqrt(mean(h .^ 2))
end