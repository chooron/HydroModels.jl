module collie2
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables S [description = "current storage in the soil moisture"]
@variables P [description = "Net precipitation intensity"]
@variables Eb [description = "bare soil evaporation"]
@variables Ev [description = "vegetation"]
@variables Ep [description = "Evapotranspiration rate"]
@variables Qse [description = "saturation excess overland flow"]
@variables Qss [description = "subsurface flow regulated by runoff coefficient a"]
@variables Qt [description = "Qse+Qss"]

# Model parameters
@parameters fc [description = "Field capacity, fc as fraction of Smax", bounds = (0.05, 0.95), unit = "mm"]
@parameters a [description = "Runoff coefficient", bounds = (0, 1), unit = "1/d"]
@parameters Smax [description = "maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters M [description = "forest fraction", bounds = (0.05, 0.95), unit = "-"]

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Ev ~ min(1.0, S / (fc * Smax)) * M * Ep
        @hydroflux Eb ~ S / Smax * (1 - M) * Ep
        @hydroflux Qse ~ step_func(S - Smax) * P
        @hydroflux Qss ~ step_func(S - (fc * Smax)) * a * (S - (fc * Smax))
    end
    dfluxes = begin
        @stateflux S ~ P - Eb - Ev - Qse - Qss
    end
end

model = @hydromodel :collie2 begin
    bucket1
    @hydroflux :qflux Qt ~ Qse + Qss
end

end