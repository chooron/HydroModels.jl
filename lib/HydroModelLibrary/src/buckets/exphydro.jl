@reexport module ExpHydro

using ..HydroModels

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5
# define variables and parameters

@variables temp lday pet prcp snowfall rainfall snowpack melt
@parameters Tmin Tmax Df Smax Qmax f

@variables soilwater pet evap baseflow surfaceflow flow rainfall

# define model components
fluxes_1 = [
    HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
]
dfluxes_1 = [StateFlux([snowfall] => [melt], snowpack),]
snow_bucket = HydroBucket(name=:exphydro_snow, funcs=fluxes_1, dfuncs=dfluxes_1)

fluxes_2 = [
    HydroFlux([soilwater, pet] => [evap], [Smax], exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
    HydroFlux([soilwater] => [baseflow], [Smax, Qmax, f], exprs=[step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))]),
    HydroFlux([soilwater] => [surfaceflow], [Smax], exprs=[max(0.0, soilwater - Smax)]),
    HydroFlux([baseflow, surfaceflow] => [flow], exprs=[baseflow + surfaceflow]),
]
dfluxes_2 = [StateFlux([rainfall, melt] => [evap, flow], soilwater)]
soil_bucket = HydroBucket(name=:exphydro_soil, funcs=fluxes_2, dfuncs=dfluxes_2)

exphydro_model = HydroModel(name=:exphydro, components=[snow_bucket, soil_bucket])

export snow_bucket, soil_bucket, exphydro_model

end
