module Baseflow
using ..HydroModels
using ..HydroModels: @variables, @parameters, Number


"""
    Constant baseflow
"""
function BASE_CONSTANT(;
    baseflow::Number=first(@variables baseflow),
    max_baseflow_rate::Number=first(@parameters max_baseflow_rate [description = "Maximum baseflow", bounds = (0, 1), unit = "mm/d"]),
    flux_name::Symbol=:base_constant
)
    @hydroflux flux_name baseflow ~ max_baseflow_rate
end

"""
    Linear storage
"""
function BASE_LINEAR(;
    baseflow::Number=first(@variables baseflow),
    waterstorage::Number=first(@variables waterstorage),
    baseflow_coeff::Number=first(@parameters baseflow_coeff [description = "Baseflow coefficient", bounds = (0, 1), unit = "mm/d"]),
    flux_name::Symbol=:base_linear
)
    @hydroflux flux_name baseflow ~ baseflow_coeff * waterstorage
end

"""
    Linear storage analytic
"""
function BASE_LINEAR_ANALYTIC(;
    baseflow::Number=first(@variables baseflow),
    waterstorage::Number=first(@variables waterstorage),
    baseflow_coeff::Number=first(@parameters baseflow_coeff [description = "Baseflow coefficient", bounds = (0, 1), unit = "mm/d"]),
    flux_name::Symbol=:base_linear_analytic
)
    @hydroflux flux_name baseflow ~ waterstorage * (1 - exp(-baseflow_coeff))
end

"""
    Non-linear storage (Power function)
"""
function BASE_POWER_LAW(;
    baseflow::Number=first(@variables baseflow),
    waterstorage::Number=first(@variables waterstorage),
    baseflow_coeff::Number=first(@parameters baseflow_coeff [description = "Baseflow coefficient", bounds = (0, 1), unit = "mm/d"]),
    baseflow_n::Number=first(@parameters baseflow_coeff [description = "User-specified Soil Parameter", bounds = (0, 1), unit = "mm/d"]),
    flux_name::Symbol=:base_power_law
)
    @hydroflux flux_name baseflow ~ clamp(baseflow_coeff * max(0.0, waterstorage)^baseflow_n, 0.0, waterstorage)
end

"""
    VIC baseflow method
"""
function BASE_VIC(;
    baseflow::Number=first(@variables baseflow),
    waterstorage::Number=first(@variables waterstorage),
    max_baseflow_rate::Number=first(@parameters max_baseflow_rate [description = "Baseflow baseflow rate", bounds = (0, 1), unit = "mm/d"]),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Maximum waterstorage", bounds = (0, 1), unit = "mm"]),
    baseflow_n::Number=first(@parameters baseflow_n [description = "User-specified Soil Parameter", bounds = (0, 1), unit = "mm/d"]),
    flux_name::Symbol=:base_vic
)
    @hydroflux flux_name baseflow ~ clamp(max_baseflow_rate * max(0.0, waterstorage / max_waterstorage)^baseflow_n, 0.0, waterstorage)
end

"""
    GR4J baseflow method
"""
function BASE_GR4J(;
    baseflow::Number=first(@variables baseflow),
    waterstorage::Number=first(@variables waterstorage),
    reference_waterstorage::Number=first(@parameters reference_waterstorage [description = "Reference waterstorage", bounds = (0, 1), unit = "mm"]),
    flux_name::Symbol=:base_gr4j
)
    @hydroflux flux_name baseflow ~ waterstorage * (1 - (1 + (waterstorage / reference_waterstorage)^4)^(1 / 4))
end

"""
    Topmodel baseflow method
"""
function BASE_TOPMODEL(;
    baseflow::Number=first(@variables baseflow),
    waterstorage::Number=first(@variables waterstorage),
    max_baseflow_rate::Number=first(@parameters max_baseflow_rate [description = "Baseflow baseflow rate", bounds = (0, 1), unit = "mm/d"]),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Maximum waterstorage", bounds = (0, 1), unit = "mm"]),
    baseflow_n::Number=first(@parameters baseflow_n [description = "user-specified soil parameter BASEFLOW_N", bounds = (0, 1), unit = "mm"]),
    lambda::Number=first(@parameters lambda [description = "mean of the power-transformed topographic index (terrain parameter TOPMODEL_LAMBDA)", bounds = (0, 1), unit = "mm"]),
    flux_name::Symbol=:base_topmodel
)
    @hydroflux flux_name baseflow ~ clamp(max_baseflow_rate * max_waterstorage / baseflow_n * (1 / lambda^baseflow_n) *
                                          max(0.0, waterstorage / max_waterstorage)^baseflow_n, 0.0, max(0.0, waterstorage))
end

"""
    Threshold-based baseflow method
"""
function BASE_THRESH_POWER(;
    baseflow::Number=first(@variables baseflow),
    waterstorage::Number=first(@variables waterstorage),
    max_baseflow_rate::Number=first(@parameters max_baseflow_rate [description = "Baseflow baseflow rate", bounds = (0, 1), unit = "mm/d"]),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Maximum waterstorage", bounds = (0, 1), unit = "mm"]),
    storage_threshold::Number=first(@parameters storage_threshold [description = "The threshold saturation at which baseflow begins", bounds = (0, 1), unit = "mm"]),
    baseflow_n::Number=first(@parameters baseflow_n [description = "user-specified soil parameter BASEFLOW_N", bounds = (0, 1), unit = "mm"]),
    flux_name::Symbol=:base_thresh_power
)
    @hydroflux flux_name baseflow ~ max_baseflow_rate * (waterstorage / max_waterstorage - storage_threshold) / max(0.0, 1 - storage_threshold)^baseflow_n
end

"""
    Threshold-based baseflow method (storage)
"""
function BASE_THRESH_STOR(;
    baseflow::Number=first(@variables baseflow),
    waterstorage::Number=first(@variables waterstorage),
    baseflow_coeff::Number=first(@parameters max_baseflow_rate [description = "Baseflow coefficient", bounds = (0, 1), unit = "mm/d"]),
    storage_threshold::Number=first(@parameters storage_threshold [description = "The threshold saturation at which baseflow begins", bounds = (0, 1), unit = "mm"]),
    flux_name::Symbol=:base_thresh_stor
)
    @hydroflux flux_name baseflow ~ baseflow_coeff * max(waterstorage - storage_threshold, 0.0)
end

end