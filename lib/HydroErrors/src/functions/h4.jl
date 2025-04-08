@doc """
    Compute the H4 mean error.

    .. image:: /pictures/H4.png
    .. image:: /pictures/MHE.png

    **Range:**

    **Notes:**

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.
        An array of observed data from the time series.

    Returns
    -------
    float
        The mean H4 error.

    References
    ----------
    - Tornquist, L., Vartia, P., Vartia, Y.O., 1985. How Should Relative Changes be Measured?
      The American Statistician 43-46.

 """
function h4_mhe(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    h = (simulated_array .- observed_array) ./ sqrt(simulated_array .* observed_array)
    mean(h)
end

@doc """
    Compute the H4 mean absolute error.

    .. image:: /pictures/H4.png
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
        The mean absolute H4 error.

    References
    ----------
    - Tornquist, L., Vartia, P., Vartia, Y.O., 1985. How Should Relative Changes be Measured?
      The American Statistician 43-46.
 """
function h4_mahe(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    h = (simulated_array .- observed_array) ./ sqrt(simulated_array .* observed_array)
    mean(abs(h))
end

@doc """
    Compute the H4 mean error.

    .. image:: /pictures/H4.png
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
        The root mean square H4 error.

    References
    ----------
    - Tornquist, L., Vartia, P., Vartia, Y.O., 1985. How Should Relative Changes be Measured?
      The American Statistician 43-46.
 """
function h4_rmshe(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    h = (simulated_array .- observed_array) ./ sqrt(simulated_array .* observed_array)
    sqrt(mean(h .^ 2))
end