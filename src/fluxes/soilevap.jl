module SoilEvap
using ..HydroModels

function SOILEVAP_ALL(;
    soil_evaporation::Number=first(@variables soil_evaporation),
    potential_evaporation::Number=first(@variables potential_evaporation),
    waterstorage::Number=first(@variables waterstorage),
    pet_corr::Number=first(@variables pet_corr [description = "Potential evapotranspiration correction factor", bounds = (0, 1)]),
)
    @hydroflux soil_evaporation ~ clamp(potential_evaporation * pet_corr, 0.0, max(0.0, waterstorage))
end

function SOILEVAP_TOPMODEL(;
    soil_evaporation::Number=first(@variables soil_evaporation),
    waterstorage::Number=first(@variables waterstorage),
    potential_evaporation::Number=first(@variables potential_evaporation),
    storage_tension::Number=first(@parameters storage_tension),
    pet_corr::Number=first(@variables pet_corr [description = "Potential evapotranspiration correction factor", bounds = (0, 1)]),
)
    @hydroflux soil_evaporation ~ clamp(potential_evaporation * pet_corr * min(waterstorage / storage_tension, 1.0), 0.0, max(0.0, waterstorage))
end

function SOILEVAP_HBV(;
    soil_evaporation::Number=first(@variables soil_evaporation),
    snow_depth::Number=first(@variables snow_depth),
    waterstorage::Number=first(@variables waterstorage),
    potential_evaporation::Number=first(@variables potential_evaporation),
    storage_tension::Number=first(@variables storage_tension),    
    pet_corr::Number=first(@variables pet_corr [description = "Potential evapotranspiration correction factor", bounds = (0, 1)]),
)
    @hydroflux soil_evaporation ~ potential_evaporation * pet_corr * ifelse(snow_depth > 0.0, 0.0, min(waterstorage / storage_tension, 1.0))
end

end