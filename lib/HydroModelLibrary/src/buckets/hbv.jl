@reexport module HBV

using ..HydroModels

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5
@variables soilwater snowpack meltwater suz slz
@variables prcp pet temp

@variables rainfall snowfall
@variables melt refreeze infil
@variables wetness excess recharge evap
@variables q0 q1 q2 q perc

@parameters TT CFMAX CFR CWH LP FC BETA PPERC UZL k0 k1 k2 kp

#* snowfall and rainfall split flux
split_flux = HydroFlux([prcp, temp] => [snowfall, rainfall], [TT], exprs=[step_func(TT - temp) * prcp, step_func(temp - TT) * prcp])
#* snowpack bucket
snow_funcs = [
    HydroFlux([temp, snowpack] => [melt], [TT, CFMAX], exprs=[min(snowpack, max(0.0, temp - TT) * CFMAX)]),
    HydroFlux([temp, meltwater] => [refreeze], [TT, CFR, CFMAX],
        exprs=[min(max((TT - temp), 0.0) * CFR * CFMAX, meltwater)]),
    HydroFlux([meltwater] => [infil], [CWH], exprs=[max(0.0, meltwater - snowpack * CWH)])
]
snow_dfuncs = [StateFlux([snowfall, refreeze] => [melt], snowpack), StateFlux([melt] => [refreeze, infil], meltwater)]
snow_bucket = HydroBucket(name=:hbv_snow, funcs=snow_funcs, dfuncs=snow_dfuncs)

#* soilwater bucket
soil_funcs = [
    HydroFlux([infil, rainfall] => [recharge], [FC, BETA], exprs=[(rainfall + infil) * clamp((max(soilwater / FC,0.0))^BETA, 0, 1)]),
    HydroFlux([soilwater] => [excess], [FC], exprs=[max(soilwater - FC, 0.0)]),
    HydroFlux([soilwater, pet] => [evap], [LP, FC], exprs=[clamp(soilwater / (LP * FC), 0, 1) * pet]),
]
soil_dfuncs = [StateFlux([rainfall, infil] => [recharge, excess, evap], soilwater)]
soil_bucket = HydroBucket(name=:hbv_soil, funcs=soil_funcs, dfuncs=soil_dfuncs)

zone_funcs = [
    HydroFlux([suz] => [perc], [PPERC], exprs=[suz * PPERC]),
    HydroFlux([suz] => [q0], [UZL, k0], exprs=[max(0.0, suz - UZL) * k0]),
    HydroFlux([suz] => [q1], [k1], exprs=[suz * k1]),
    HydroFlux([slz] => [q2], [k2], exprs=[slz * k2]),
    HydroFlux([q0, q1, q2] => [q], exprs=[q0 + q1 + q2]),
]

zone_dfuncs = [
    StateFlux([recharge, excess] => [perc, q0, q1], suz),
    StateFlux([perc] => [q2], slz),
]

zone_bucket = HydroBucket(
    name=:hbv_zone,
    funcs=zone_funcs,
    dfuncs=zone_dfuncs,
)

hbv_model = HydroModel(name=:dpl_hbv, components=[split_flux, snow_bucket, soil_bucket, zone_bucket])

export split_flux, snow_bucket, soil_bucket, zone_bucket, hbv_model
end