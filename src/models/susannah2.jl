module susannah2
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables P [description = "Current precipitation", unit = "mm"]
@variables Ep [description = "Evapotranspiration rate", unit = "mm"]
@variables Esat [description = "Evaporation from the saturated zone", unit = "mm"]
@variables Eus [description = "Evaporation from the unsaturated zone", unit = "mm"]
@variables Sus [description = "Current storage in the unsaturated store", unit = "mm"]
@variables Susfc [description = "Variable field capacity", unit = "mm"]
@variables Ssat [description = "Current storage on the saturated zone", unit = "mm"]
@variables Se [description = "Storage excess", unit = "mm"]
@variables D [description = "Soil depth",  unit ="mm"]
@variables rg [description = "Drainage",  unit ="mm"]
@variables Qse [description = "saturation excess runoff", unit = "mm"]
@variables Qr [description = "Recharge of deep groundwater", unit = "mm"]
@variables Qss [description = "Subsurface flow", unit = "mm"]
@variables Qt [description = "Total runoff", unit = "mm"]

# Model parameters
@parameters Smax [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters phi [description = "Porosity", bounds = (0.05, 0.95), unit = "-"]
@parameters fc [description = "Field capacity as fraction of Smax", bounds = (0.05, 0.95), unit = "-"]
@parameters r [description = "Fraction of subsurface outflow to deep groundwater", bounds = (0, 1), unit = "-"]
@parameters c [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]
@parameters d [description = "Runoff nonlinearity", bounds = (1, 5), unit = "-"]

bucket = @hydrobucket :bucket2 begin
    fluxes = begin
        @hydroflux Eus ~ Sus / Smax * Ep
        @hydroflux Susfc ~ (Smax - Ssat) * fc / phi
        @hydroflux rg ~ step_func(Sus - Susfc) * P
        @hydroflux Se ~ max(Sus - Susfc, 0)

        @hydroflux Esat ~ Ssat / Smax * Ep
        @hydroflux Qse ~ step_func(Ssat - Smax) * (rg + Se)
        @hydroflux Qss ~ min(Ssat, (1 - r) * c * abs(Ssat) ^ d)
        @hydroflux Qr ~ min(Ssat, r * c * abs(Ssat) ^ d)
    end
    dfluxes = begin
        @stateflux Sus ~ P - Eus - rg - Se
        @stateflux Ssat ~ rg - Esat - Qse - Qss - Qr
    end
end

model = @hydromodel :susannah2 begin
    bucket
    @hydroflux Qt ~ Qse + Qss
end

end