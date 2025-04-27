module collie3
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables P [description = "the precipitation input", unit = "mm/d"]
@variables Ep [description = "potential evapotranspiration", unit = "mm/d"]

@variables S [description = "current storage in the soil moisture", unit = "mm"]
@variables G [description = "groundwater storage", unit = "mm"]

@variables Eb [description = "between bare soil evaporation", unit = "mm/d"]
@variables Ev [description = "transpiration through vegetation", unit = "mm/d"]
@variables Qse [description = "saturation excess overland flow", unit = "mm/d"]
@variables Qss [description = "non-linear subsurface flow regulated by runoff coefficients a and b", unit = "mm/d"]
@variables Qsg [description = "non-linear groundwater flow that relies on the same parameters as subsurface flow uses", unit = "mm"]
@variables Qt [description = "total runoff", unit = "mm"]

# Model parameters
@parameters Smax [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters Sfc [description = "Field capacity", bounds = (0.05, 0.95), unit = "mm"]
@parameters a [description = "Runoff coefficient", bounds = (0, 1), unit = "1/d"]
@parameters M [description = "Forest fraction", bounds = (0.05, 0.95), unit = "-"]
@parameters b [description = "Runoff nonlinearity", bounds = (1, 5), unit = "-"]
@parameters λ [description = "Fraction subsurface flow to groundwater", bounds = (0, 1), unit = "-"]

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Eb ~ S / Smax * (1 - M) * Ep
        @hydroflux Ev ~ min(1.0, S / Sfc) * M * Ep
        @hydroflux Qse ~ step_func(S - Smax) * P
        @hydroflux Qss ~ (a * max(0.0, S - Sfc))^b
        @hydroflux Qsg ~ (a * max(0.0, G))^b
        @hydroflux Qt ~ Qse + (1 - λ) * Qss + Qsg

    end
    dfluxes = begin
        @stateflux S ~ P - Eb - Ev - Qse - Qss
        @stateflux G ~ Qss * λ - Qsg
    end
end

model = @hydromodel :collie3 begin
    bucket1
end

end