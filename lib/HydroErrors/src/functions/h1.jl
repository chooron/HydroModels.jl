@doc """
    Compute the H1 mean error.

    .. image:: /pictures/H1.png
    .. image:: /pictures/MHE.png

    **Range:**

    **Notes:**

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray

    Returns
    -------
    float
        The mean H1 error.

    References
    ----------
    - Tornquist, L., Vartia, P., Vartia, Y.O., 1985. How Should Relative Changes be Measured?
      The American Statistician 43-46.
 """
function h1_mhe(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    h = (simulated_array .- observed_array) ./ observed_array
    mean(h)
end

@doc """
    Compute the H1 absolute error.

    .. image:: /pictures/H1.png
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
        The H1 absolute error.

    References
    ----------
    - Tornquist, L., Vartia, P., Vartia, Y.O., 1985. How Should Relative Changes be Measured?
      The American Statistician 43-46.
 """
function h1_mahe(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    h = (simulated_array .- observed_array) ./ observed_array
    mean(abs(h))
end

@doc """
    Compute the H1 root mean square error.

    .. image:: /pictures/H1.png
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
        The root mean squared H1 error.

    References
    ----------
    - Tornquist, L., Vartia, P., Vartia, Y.O., 1985. How Should Relative Changes be Measured?
      The American Statistician 43-46.
 """
function h1_rmshe(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    h = (simulated_array .- observed_array) ./ observed_array
    sqrt(mean(h .^ 2))
end