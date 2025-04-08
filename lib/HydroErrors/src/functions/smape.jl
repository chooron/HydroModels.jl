@doc """
    Compute the the Symmetric Mean Absolute Percentage Error (1) (SMAPE1).

    .. image:: /pictures/SMAPE1.png

    **Range:** 0 ≤ SMAPE1 < 100%, smaller is better, symmetrical.

    **Notes:** This metric is an adjusted version of the MAPE.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The symmetric mean absolute percentage error (1).

    References
    ----------
    - Flores, B.E., 1986. A pragmatic view of accuracy measurement in forecasting. Omega 14(2)
      93-98.
    - Goodwin, P., Lawton, R., 1999. On the asymmetry of the symmetric MAPE. International Journal
      of Forecasting 15(4) 405-408.
 """
function smape1(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    a = 100 ./ simulated_array.size
    b = abs(simulated_array .- observed_array)
    c = abs(simulated_array) .+ abs(observed_array)
    a .* sum(b ./ c)
end




@doc """
    Compute the the Symmetric Mean Absolute Percentage Error (2) (SMAPE2).

    .. image:: /pictures/SMAPE2.png

    **Range:** 0 ≤ SMAPE1 < 200%, does not indicate bias, smaller is better, symmetrical.

    **Notes:** This metric is an adjusted version of the MAPE with only positive metric values.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The symmetric mean absolute percentage error (2).

    References
    ----------
    - Flores, B.E., 1986. A pragmatic view of accuracy measurement in forecasting. Omega 14(2)
      93-98.
    - Goodwin, P., Lawton, R., 1999. On the asymmetry of the symmetric MAPE. International Journal
      of Forecasting 15(4) 405-408.
 """
function smape2(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    a = simulated_array .- observed_array
    b = (simulated_array .+ observed_array) ./ 2
    c = 100 ./ simulated_array.size
    c .* sum(abs(a ./ b))
end