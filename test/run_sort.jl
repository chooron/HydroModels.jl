@testset "test flux sort function" begin
    @variables a b c d e f
    @parameters p1 p2 p3 p4

    flux_1 = HydroModels.HydroFlux([a, b] => [c, d], [p1, p2], exprs=[a * p1 + p2, b * p2 + p1])
    flux_2 = HydroModels.HydroFlux([a, c] => [e], [p3], exprs=[a * p3 + c])
    flux_3 = HydroModels.HydroFlux([e, d] => [f], exprs=[e + d])
    sorted_fluxes = HydroModels.sort_fluxes([flux_2, flux_3, flux_1])
    @test sorted_fluxes == [flux_1, flux_2, flux_3]
end

@test "test sort components" begin
    @variables temp lday pet prcp snowfall rainfall snowpack melt
    @parameters Tmin Tmax Df
    name = :test
    snow_fluxes = [
        HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
        HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
        HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
    ]
    snow_dfluxes = [StateFlux([snowfall] => [melt], snowpack),]
    snow_bucket = HydroBucket(Symbol(name, :_surface), fluxes=snow_fluxes, dfluxes=snow_dfluxes)

    @variables soilwater pet evap baseflow surfaceflow flow rainfall melt
    @parameters Smax Qmax f
    soil_fluxes = [
        HydroFlux([soilwater, pet] => [evap], [Smax], exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
        HydroFlux([soilwater] => [baseflow], [Smax, Qmax, f], exprs=[step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))]),
        HydroFlux([soilwater] => [surfaceflow], [Smax], exprs=[max(0.0, soilwater - Smax)]),
    ]
    soil_dfluxes = [StateFlux([rainfall, melt] => [evap, baseflow, surfaceflow], soilwater)]
    soil_bucket = HydroBucket(Symbol(name, :_soil), fluxes=soil_fluxes, dfluxes=soil_dfluxes)
    flow_flux = HydroFlux([baseflow, surfaceflow] => [flow], exprs=[baseflow + surfaceflow])
    sorted_fluxes = HydroModels.sort_components([flow_flux, soil_bucket, snow_bucket])
end