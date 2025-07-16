module SnowFreeze
using ..HydroModels

function FREEZE_DEGREE_DAY(
    snowfreeze::Num=first(@variables snowfreeze),
    liquidwater::Num=first(@variables liquidwater),
    temp::Num=first(@variables temp),
    kf::Num=first(@parameters kf [description = "Degree-day factor for freezing", bounds = (0, 1), unit = "mm/degree celsius/day"]),
    tf::Num=first(@parameters tf [description = "Freezing temperature", bounds = (-10, 10), unit = "degree celsius"]),
    flux_name::Symbol=:freeze_degree_day,
)
    @hydroflux flux_name snowfreeze ~ min(max(0.0, liquidwater), kf * max(0, tf - temp))
end

function FREEZE_HMETS(
    snowfreeze::Num=first(@variables snowfreeze),
    liquidwater::Num=first(@variables liquidwater),
    mean_temp::Num=first(@variables mean_temp),
    min_temp::Num=first(@variables min_temp),
    kf::Num=first(@parameters kf [description = "Degree-day factor for freezing", bounds = (0, 1), unit = "mm/degree celsius/day"]),
    tf::Num=first(@parameters tf [description = "Freezing temperature", bounds = (-10, 10), unit = "degree celsius"]),
    re::Num=first(@parameters re [description = "land use parameter REFREEZE_EXP", bounds = (-10, 10), unit = "-"]),
    flux_name::Symbol=:freeze_degree_day,
)
    @hydroflux flux_name snowfreeze ~ min(max(0.0, liquidwater), kf * max(0, trf - (mean_temp + min_temp) / 2)^refreeze_exp)
end

end