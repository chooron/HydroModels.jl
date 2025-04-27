module unitedstates
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables P [description = "precipitation", unit = "mm/d"]
@variables Ep [description = "Evaporation rate", unit = "mm/d"]
@variables Sus [description = "current storage in the unsaturated zone", unit = "mm"]
@variables Ssat [description = "current storage in the saturated zone", unit = "mm"]
@variables Eus_ei [description = "evaporation from interception", unit = "mm/d"]
@variables Eus_veg [description = "transpiration through vegetation", unit = "mm/d"]
@variables Eus_bs [description = "bare soil evaporation", unit = "mm/d"]
@variables rg [description = "drainage to the saturated zone", unit = "mm/d"]
@variables Susfc [description = "estimated field capacity", unit = "-"]
@variables Se [description = "storage excess", unit = "-"]
@variables Esat_veg [description = "transpiration through vegetation", unit = "mm/d"]
@variables Esat_bs [description = "bare soil evaporation", unit = "mm/d"]
@variables Qse [description = "saturation excess overland flow", unit = "mm/d"]
@variables Qss [description = "subsurface flow", unit = "mm/d"]
@variables Qt [description = "Total runoff", unit = "mm/d"]

# Model parameters
@parameters Alpha_ei [description = "Intercepted fraction of precipitation", bounds = (0, 1), unit = "-"]
@parameters M [description = "Forest fraction", bounds = (0.05, 0.95), unit = "-"]
@parameters Smax [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters fc [description = "Field capacity as fraction of Smax", bounds = (0.05, 0.95), unit = "-"]
@parameters Alpha_ss [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]

# Soil water component
bucket = @hydrobucket :bucket begin
    fluxes = begin
        @hydroflux Eus_ei ~ Alpha_ei * P
        @hydroflux Susfc ~ fc * (Smax - Ssat)
        @hydroflux Eus_veg ~ min(Sus / Susfc, 1.0) * Sus * M * Ep / (Sus + Ssat + 1e-6)
        @hydroflux Eus_bs ~ (Sus / (Sus + Ssat+ 1e-6)) * (1 - M) * (Sus / (Smax - Ssat)) * Ep
        @hydroflux rg ~ P * step_func(Sus - Susfc)
        @hydroflux Se ~ max(0.0, Sus - Susfc)

        @hydroflux Esat_veg ~ (Ssat / Smax) * M * Ep
        @hydroflux Esat_bs ~ (Ssat / Smax) * (1 - M) * Ep
        @hydroflux Qse ~ rg * step_func(Sus - Smax)
        @hydroflux Qss ~ Alpha_ss * Ssat
        @hydroflux Qt ~ Qse + Qss
    end

    dfluxes = begin
        @stateflux Sus ~ P - Eus_ei - Eus_veg - Eus_bs - rg
        @stateflux Ssat ~ rg - Esat_veg - Esat_bs - Qse - Qss
    end
end

model = @hydromodel :unitedstates begin
    bucket
end
end

