module australia
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables P [description = "current precipitation", unit = "mm/d"]
@variables Ep [description = "Evaporation rate", unit = "mm/d"]

@variables Sus [description = "current storage", unit = "mm"]
@variables Ssat [description = "maximum soil moisture storage", unit = "mm"]
@variables Gw [description = "current groundwater storage", unit = "mm"]

@variables rg [description = "drainage from the unsaturated store to the saturated store", unit = "mm/d"]
@variables Se [description = "storage excess", unit = "mm/d"]
@variables Eus [description = "evaporation from the unsaturated store", unit = "mm"]
@variables Esat [description = "evaporation from the saturated store", unit = "mm"]
@variables QSE [description = "saturation excess", unit = "mm/d"]
@variables QSS [description = "subsurface flow", unit = "mm/d"]
@variables QR [description = "recharge of deep groundwater ", unit = "mm/d"]
@variables QBF [description = "baseflow", unit = "mm/d"]
@variables Qt [description = "Total runoff", unit = "mm/d"]

# Model parameters
@parameters Smax [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters φ [description = "Porosity", bounds = (0.05, 0.95), unit = "-"]
@parameters fc [description = "Field capacity", bounds = (0.01, 1.00), unit = "-"]
@parameters αss [description = "Runoff coefficient", bounds = (0, 1.00), unit = "1/d"]
@parameters βss [description = "Runoff nonlinearity", bounds = (1, 5), unit = "-"]
@parameters Kdeep [description = "Runoff coefficient", bounds = (0, 1.00), unit = "1/d"]
@parameters αbf [description = "Capillary rise rate", bounds = (0, 1.00), unit = "1/d"]
@parameters βbf [description = "Runoff nonlinearity", bounds = (1, 5), unit = "-"]

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Eus ~ Sus / Smax * Ep
        @hydroflux rg ~ step_func(Sus - (Smax - Ssat) * fc / φ) * P
        @hydroflux Se ~ max(Sus - (Smax - Ssat) * fc / φ, 0)
        @hydroflux Esat ~ Ssat / Smax * Ep
        @hydroflux QSE ~ step_func(Ssat - Smax) * (rg + Se)
        @hydroflux QSS ~ αss * max(Ssat, 0)^βss
        @hydroflux QR ~ Kdeep * Ssat
    end
    dfluxes = begin
        @stateflux Sus ~ P - Eus - rg - Se
        @stateflux Ssat ~ rg - Esat - QSE - QSS - QR
    end
end

bucket2 = @hydrobucket :bucket2 begin
    fluxes = begin
        @hydroflux QBF ~ αbf * max(Gw, 0)^βbf
    end
    dfluxes = begin
        @stateflux Gw ~ QR - QBF
    end
end

flux1 = @hydroflux Qt ~ QSE + QBF + QSS

model = @hydromodel :australia begin
    bucket1
    bucket2
    flux1
end
end