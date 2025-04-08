@doc """
    Compute the the anomaly correlation coefficient (ACC).

    .. image:: /pictures/ACC.png

    **Range:** -1 ≤ ACC ≤ 1. -1 indicates perfect negative correlation of the variation
    pattern of the anomalies, 0 indicates complete randomness of the variation patterns of the
    anomalies, 1 indicates perfect correlation of the variation pattern of the anomalies.

    **Notes:** Common measure in the verification of spatial fields. Measures the correlation
    between the variation pattern of the simulated data compared to the observed data.

    Parameters
    ----------
    simulated_array: one dimensional ndarray
        An array of simulated data from the time series.

    observed_array: one dimensional ndarray
        An array of observed data from the time series.

    Returns
    -------
    float
        The anomaly correlation coefficient.

    References
    ----------
    - Langland, Rolf H., and Ryan N. Maue. “Recent Northern Hemisphere Mid-Latitude Medium-Range
      Deterministic Forecast Skill.” Tellus A: Dynamic Meteorology and Oceanography 64,
      no. 1 (2012): 17531.
    - Miyakoda, K., G. D. Hembree, R. F. Strickler, and I. Shulman. “Cumulative Results of Extended
      Forecast Experiments I. Model Performance for Winter Cases.” Monthly Weather Review 100, no.
      12(1972): 836–55.
    - Murphy, Allan H., and Edward S. Epstein. “Skill Scores and Correlation Coefficients in Model
      Verification.” Monthly Weather Review 117, no. 3 (1989): 572–82.
 """
function acc(simulated_array::AbstractVector, observed_array::AbstractVector; kwargs...)
    a = simulated_array .- mean(simulated_array)
    b = observed_array .- mean(observed_array)
    c = std(observed_array, ddof=1) .* std(simulated_array, ddof=1) .* simulated_array.size
    dot(a, b ./ c)
end