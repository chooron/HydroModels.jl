module Infiltration
using ..HydroModels

"""
    INF_HMETS
"""
function INF_HMETS(;
    infiltration::Number=first(@variables infiltration),
    prcpeff::Number=first(@variables prcpeff),
    waterstorage::Number=first(@variables waterstorage),
    runoff_coeff::Number=first(@parameters runoff_coeff [description = "Runoff coefficient", bounds = (0, 1), unit = "mm/d"]),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Maximum waterstorage", bounds = (0, 1), unit = "mm"]),
    flux_name::Symbol=:inf_hmets,
)
    @hydroflux flux_name infiltration ~ prcpeff * max(0.0, 1 - runoff_coeff * max(0.0, waterstorage) / max_waterstorage)
end

"""
    INF_VIC_ARNO
"""
function INF_VIC_ARNO(;
    infiltration::Number=first(@variables infiltration),
    prcpeff::Number=first(@variables prcpeff),
    waterstorage::Number=first(@variables waterstorage),
    exp_coeff::Number=first(@parameters exp_coeff [description = "exp coefficient", bounds = (0, 1), unit = "mm/d"]),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Maximum waterstorage", bounds = (0, 1), unit = "mm"]),
    flux_name::Symbol=:inf_hmets,
)
    @hydroflux flux_name infiltration ~ prcpeff * max(0.0, 1 - max(0.0, 1 - max(0.0, waterstorage) / max_waterstorage)^exp_coeff)
end

"""
    INF_VIC_ARNO
"""
function INF_HBV(;
    infiltration::Number=first(@variables infiltration),
    prcpeff::Number=first(@variables prcpeff),
    waterstorage::Number=first(@variables waterstorage),
    hbv_beta::Number=first(@parameters hbv_beta [description = "parameter beta in HBV model", bounds = (0, 1), unit = "mm/d"]),
    max_waterstorage::Number=first(@parameters max_waterstorage [description = "Maximum waterstorage", bounds = (0, 1), unit = "mm"]),
    flux_name::Symbol=:inf_hmets,
)
    @hydroflux flux_name infiltration ~ prcpeff * max(0.0, 1 - max(0.0, max(0.0, waterstorage) / max_waterstorage)^hbv_beta)
end
end