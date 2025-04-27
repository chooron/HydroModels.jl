module alpine2
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables P [description = "precipitation", unit = "mm/d"]
@variables T [description = "temperature", unit = "celsius"]
@variables Ep [description = "evaporate at the potential rate", unit = "mm/d"]

@variables Sn [description = "current snow storage", unit = "mm"]
@variables S [description = "current soil moisture storage", unit = "mm"]
@variables Ps [description = "precipitation that falls as snow", unit = "mm/d"]

@variables Pr [description = "rainfall", unit = "mm/d"]
@variables Melt [description = "snow melt", unit = "mm/d"]
@variables Ea [description = "evaporate at the actual rate", unit = "mm/d"]
@variables Qse [description = "saturation excess runoff", unit = "mm/d"]
@variables Qin [description = "interflow", unit = "mm/d"]
@variables Qbf [description = "baseflow", unit = "mm/d"]
@variables Qt [description = "Total runoff", unit = "mm/d"]

# Model parameters
@parameters Tt [description = "Threshold temperature for snowfall and melt", bounds = (-3, 5), unit = "celsius"]
@parameters ddf [description = "Degree-day factor", bounds = (0, 20), unit = "mm/degree celsius/d"]
@parameters Smax [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters Sfc [description = "Field capacity", bounds = (0.05, 0.95), unit = "mm"]
@parameters tc_in [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]
@parameters tc_bf [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]

# Soil water component
bucket_1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Ps ~ P * step_func(Tt - T)
        @hydroflux Melt ~ min(Sn, ddf * (T - Tt) * step_func(T - Tt))
    end

    dfluxes = begin
        @stateflux Sn ~ Ps - Melt
    end
end

bucket_2 = @hydrobucket :bucket2 begin
    fluxes = begin
        @hydroflux Pr ~ P * step_func(T - Tt)
        @hydroflux Ea ~ Ep * step_func(S)
        @hydroflux Qse ~ (Pr + Melt) * step_func(S - Smax)
        @hydroflux Qin ~ tc_in * max(S - Sfc, 0.00)
        @hydroflux Qbf ~ tc_bf * S
        @hydroflux Qt ~ Qse + Qin + Qbf
    end

    dfluxes = begin
        @stateflux S ~ Pr + Melt - Ea - Qse - Qin - Qbf
    end
end


model = @hydromodel :alpine2 begin
    bucket_1
    bucket_2
end

end