module newzealand2

using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables Si [description = "current interception storage", unit = "mm"]
@variables Sm [description = "current soil moisture storage", unit = "mm"]

@variables P [description = "daily precipitation", unit = "mm/d"]
@variables Ep [description = "evaporate rate", unit = "mm/d"]

@variables Eint [description = "Intercepted water", unit = "mm/d"]
@variables Qtf [description = "through-fall towards soil moisture", unit = "mm/d"]
@variables Eveg [description = "Evaporation through vegetation", unit = "mm/d"]
@variables Ebs [description = "bare soil evaporation", unit = "mm/d"]
@variables Qse [description = "saturation excess runoff", unit = "mm/d"]
@variables Qss [description = "subsurface runoff", unit = "mm/d"]
@variables Qbf [description = "baseflow", unit = "mm/d"]
@variables Qt [description = "Total runoff", unit = "mm/d"]

# Model parameters
@parameters Imax [description = "Maximum interception capacity", bounds = (0, 5), unit = "mm"]
@parameters Smax [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters Sfc [description = "Field capacity", bounds = (0.05, 0.95), unit = "-"]
@parameters M [description = "Forest fraction", bounds = (0.05, 0.95), unit = "-"]
@parameters a [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]
@parameters b [description = "Runoff nonlinearity", bounds = (1, 5), unit = "-"]
@parameters tc_bf [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]
@parameters d [description = "Unit Hydrograph time base", bounds = (1, 120), unit = "d"]

# Soil water component
bucket_1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Eint ~ Ep
        @hydroflux Qtf ~ P * step_func(Si - Imax)
    end

    dfluxes = begin
        @stateflux Si ~ P - Eint - Qtf
    end
end

bucket_2 = @hydrobucket :bucket2 begin
    fluxes = begin
        @hydroflux Eveg ~ M * Ep * min(Sm / Sfc, 1.0)
        @hydroflux Ebs ~ (Sm / Smax) * (1 - M) * Ep
        @hydroflux Qse ~ P * step_func(Sm - Smax)
        @hydroflux Qss ~ (a * max(Sm - Sfc, 0.00))^b
        @hydroflux Qbf ~ tc_bf * Sm
        @hydroflux Qt ~ Qse + Qss + Qbf
    end

    dfluxes = begin
        @stateflux Sm ~ Qtf - Eveg - Ebs - Qse - Qss - Qbf
    end
end

model = @hydromodel :newzealand2 begin
    bucket_1
    bucket_2
end

end


