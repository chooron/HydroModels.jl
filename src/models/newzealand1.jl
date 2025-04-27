module newzealand1
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables S [description = "current soil moisture storage", unit = "mm"]
@variables P [description = "precipitation", unit = "mm/d"]
@variables Ep [description = "Evaporation rate", unit = "mm/d"]
@variables Eveg [description = "Evaporation through vegetation", unit = "mm/d"]
@variables Ebs [description = "bare soil evaporation", unit = "mm/d"]
@variables Qss [description = "subsurface runoff", unit = "mm/d"]
@variables Qse [description = "saturation excess runoff", unit = "mm/d"]
@variables Qbf [description = "baseflow", unit = "mm/d"]
@variables Qt [description = "Total runoff", unit = "mm/d"]

# Model parameters
@parameters Smax [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters Sfc [description = "Field capacity as fraction of maximum soil moisture", bounds = (0.05, 0.95), unit = "-"]
@parameters M [description = "Fraction forest", bounds = (0.05, 0.95), unit = "-"]
@parameters a [description = "Subsurface runoff coefficient", bounds = (0, 1), unit = "d-1"]
@parameters b [description = "Runoff non-linearity", bounds = (1, 5), unit = "-"]
@parameters tc_bf [description = "Baseflow runoff coefficient", bounds = (0, 1), unit = "d-1"]

# Soil water component
soil_bucket = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux Eveg ~ min(S / Sfc, 1.0) * M * Ep
        @hydroflux Ebs ~ S / Smax * (1 - M) * Ep
        @hydroflux Qse ~ P * step_func(S - Smax)
        @hydroflux Qss ~ (a * max(0.0, S - Sfc))^b
        @hydroflux Qbf ~ tc_bf * S
        @hydroflux Qt ~ Qse + Qss + Qbf
    end

    dfluxes = begin
        @stateflux S ~ P - Eveg - Ebs - Qse - Qss - Qbf
    end
end
model = @hydromodel :newzealand1 begin
    soil_bucket
end

end


