module PotentialMelt

using ..HydroModels

function POTMELT_DEGREE_DAY(;
    potential_melt::Number=first(@variables potential_melt),
    temp::Number=first(@variables temp),
    pet_corr::Number=first(@variables pet_corr [description = "Potential evapotranspiration correction factor", bounds = (0, 1)]),
    tf::Number=first(@parameters tf [description = "Freeze/melt temperature", bounds = (0, 1), unit = "degree celsius"]),
    ma::Number=first(@parameters ma [description = "Melt factor", bounds = (-10, 10), unit = "mm/day/̧degree celsius"]),
    flux_name::Symbol=:potmelt_degree_day,
)
    @hydroflux flux_name potential_melt ~ ma * max(0, temp - tf) * pet_corr
end


end