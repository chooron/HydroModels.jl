module alpine1
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables P [description = "Incoming precipitation", unit = "[mm/d]"]
@variables T [description = "Temperature", unit = "[degree celsius]"]
@variables Ps [description = "Precipitation that falls as snow", unit = "[mm/d]"]
@variables Pr [description = "Precipitation that falls as rain", unit = "[mm/d]"]
@variables Sn [description = "Current snow storage", unit = "[mm]"]
@variables Sm [description = "Current soil moisture storage", unit = "[mm]"]
@variables Qt [description = "Total runoff", unit = "[mm/d]"]
@variables QN [description = "Snow melt", unit = "[mm/d]"]
@variables Qse [description = "saturation excess runof", unit = "[mm/d]"]
@variables Qss [description = "Subsurface flow", unit = "[mm/d]"]
@variables Ea [description = "Evaporate rate", unit = "[mm/d]"]
@variables Ep [description = "The potential rate", unit = "[mm/d]"]

# Model parameters
@parameters Tt [description = "Threshold temperature for snowfall and melt", bounds = (-3, 5), unit = "[degree celsius]"]
@parameters ddf [description = "Degree-day factor", bounds = (0, 20), unit = "[mm/degree celsius/day]"]
@parameters Smax [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "[mm]"]
@parameters tc [description = "Runoff coefficient", bounds = (0, 1), unit = "[d-1]"]


function build_upper_soil_bucket(;
    name::Optional{Symbol}=nothing,
    P=P, T=T, Tt=Tt,
    ddf=ddf, Sn=Sn,
)
    return @hydrobucket name begin
        fluxes = begin
            @hydroflux Ps ~ step_func(Tt - T) * P
            @hydroflux QN ~ min(step_func(T - Tt) * ddf * (T - Tt), Sn)
        end
        dfluxes = begin
            @stateflux Sn ~ Ps - QN
        end
    end
end
# Soil water component
function build_lower_soil_bucket(;
    name::Optional{Symbol}=nothing,
    P=P, T=T, Tt=Tt,
    ddf=ddf, Sn=Sn,
    Sm=Sm, Smax
)
    return @hydrobucket name begin
        fluxes = begin
            @hydroflux Pr ~ step_func(T - Tt) * P
            @hydroflux Ea ~ step_func(Sm) * Ep
            @hydroflux Qse ~ step_func(Sm - Smax) * (Pr + QN)
            @hydroflux Qss ~ tc * Sm
        end
        dfluxes = begin
            @stateflux Sm ~ Pr + QN - Ea - Qse - Qss
        end
    end
end

model = @hydromodel :alpine1 begin
    build_upper_soil_bucket()
    build_lower_soil_bucket()
    @hydroflux Qt ~ Qse + Qss
end

end