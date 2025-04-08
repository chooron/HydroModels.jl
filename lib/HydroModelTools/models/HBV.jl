using ModelingToolkit

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
snow_fluxes = [
    HydroFlux([temp, snowpack] => [melt], [TT, CFMAX], exprs=[min(snowpack, max(0.0, temp - TT) * CFMAX)]),
    HydroFlux([temp, meltwater] => [refreeze], [TT, CFR, CFMAX],
        exprs=[min(max((TT - temp), 0.0) * CFR * CFMAX, meltwater)]),
    HydroFlux([meltwater] => [infil], [CWH], exprs=[max(0.0, meltwater - snowpack * CWH)])
]
snow_dfluxes = [StateFlux([snowfall, refreeze] => [melt], snowpack), StateFlux([melt] => [refreeze, infil], meltwater)]
snow_bucket = HydroBucket(name=:hbv_snow, fluxes=snow_fluxes, dfluxes=snow_dfluxes)

#* soilwater bucket
soil_fluxes = [
    HydroFlux([infil, rainfall] => [recharge], [FC, BETA], exprs=[(rainfall + infil) * clamp((max(soilwater / FC,0.0))^BETA, 0, 1)]),
    HydroFlux([soilwater] => [excess], [FC], exprs=[max(soilwater - FC, 0.0)]),
    HydroFlux([soilwater, pet] => [evap], [LP, FC], exprs=[clamp(soilwater / (LP * FC), 0, 1) * pet]),
]
soil_dfluxes = [StateFlux([rainfall, infil] => [recharge, excess, evap], soilwater)]
soil_bucket = HydroBucket(name=:hbv_soil, fluxes=soil_fluxes, dfluxes=soil_dfluxes)

zone_fluxes = [
    HydroFlux([suz] => [perc], [PPERC], exprs=[suz * PPERC]),
    HydroFlux([suz] => [q0], [UZL, k0], exprs=[max(0.0, suz - UZL) * k0]),
    HydroFlux([suz] => [q1], [k1], exprs=[suz * k1]),
    HydroFlux([slz] => [q2], [k2], exprs=[slz * k2]),
    HydroFlux([q0, q1, q2] => [q], exprs=[q0 + q1 + q2]),
]

zone_dfluxes = [
    StateFlux([recharge, excess] => [perc, q0, q1], suz),
    StateFlux([perc] => [q2], slz),
]

zone_bucket = HydroBucket(
    name=:hbv_zone,
    fluxes=zone_fluxes,
    dfluxes=zone_dfluxes,
)

hbv_model = HydroModel(name=:dpl_hbv, components=[split_flux, snow_bucket, soil_bucket, zone_bucket])