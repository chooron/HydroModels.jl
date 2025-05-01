module gsmsocont
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables His [description = "current storage in the snow pack", unit = "mm"]
@variables Pices [description = "precipitation-as-snow", unit = "mm/d"]
@variables Mis [description = "melt", unit = "mm/d"]
@variables T [description = "temperature", unit = "oC"]
@variables P [description = "precipitation", unit = "mm/d"]
@variables Ep [description = "precipitation", unit = "mm/d"]

@variables Pice [description = "falls on the ice-covered part of the catchment", unit = "-"]
@variables Mis [description="uses a degree-day-factor asnow to estimate snow melt", unit="-"]
@variables Ts [description = "r snow melt", unit = "oC"]

@variables Sis [description = "current storage in the snow-water routing reservoir", unit = "mm"]
@variables Mis [description = "snow melt", unit = "mm/d"]
@variables Qis [description = "runoff", unit = "-"]
@variables Pirs [description = "occurs only if the current snow pack storage is above zero", unit = "-"]
@variables Picer [description = "precipitation-as-rain", unit = "-"]
@variables Pices [description = "precipitation-as-rain", unit = "-"]
@variables Pnir [description = "precipitation-as-rain", unit = "-"]
@variables Piri [description = "precipitation-as-rain", unit = "-"]

@variables Sii [description = "current storage in the ice-water routing reservoir", unit = "mm"]
@variables Mii [description = "glacier melt", unit = "mm/d"]
@variables Picei [description="rain-on-ice", unit="mm/d"]
@variables Qii [description = "runoff", unit = "mm/d"]
@variables His [description = "snow pack", unit = "-"]

@variables Hnis [description = "current snow pack storage", unit = "mm"]
@variables Pnis [description = "increases through snowfall", unit = "mm/d"]
@variables Mnis [description = "increases through snowfall", unit = "mm/d"]

@variables Snis [description = "current storage in soil moisture", unit = "mm"]
@variables Pinf [description = "infiltrated precipitation", unit = "mm/d"]
@variables ET [description = "evapotranspiration", unit = "mm/d"]
@variables Qsl [description = "slow flow", unit = "mm/d"]
@variables Pinf [description = "depends on the effective precipitation", unit = "-"]
@variables Peff [description = "precipitation", unit = "-"]
@variables Peq [description = "total of snow melt", unit = "-"]
@variables Mnis [description = "snow melt", unit = "-"]

@variables Sniq [description = "current storage in the direct runoff reservoir", unit = "mm"]
@variables Peff [description = "effective precipitation", unit = "mm/d"]
@variables Qqu [description = "quick flow", unit = "mm/d"]
@variables Pnonice [description = "", unit = ""]
@variables Qt [description = "", unit = ""]

# Model parameters
@parameters fice [description = "Fraction of catchment covered by glacier", bounds = (0, 1), unit = "-"]
@parameters t0 [description = "Threshold temperature for snowfall", bounds = (-3, 5), unit = "oC"]
@parameters asnow [description = "Degree-day factor for snow melt", bounds = (0, 20), unit = "mm/oC/d"]
@parameters tm [description = "Threshold temperature for snow melt", bounds = (-3, 3), unit = "oC"]
@parameters ks [description = "Shape parameter for evaporation", bounds = (0, 1), unit = "1/d"]
@parameters aice [description = "Degree-day factor for ice melt", bounds = (0, 20), unit = "mm/oC/d"]
@parameters ki [description = "Runoff coeficient for ice melt on glacier", bounds = (0, 1), unit = "1/d"]
@parameters a [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters x [description = "Evaporation non-linearity", bounds = (0, 10), unit = "-"]
@parameters y [description = "Infiltration non-linearity", bounds = (0, 5), unit = "-"]
@parameters ksl [description = "Routing delay  ", bounds = (0, 1), unit = "1/d"]
@parameters beta [description = "Runoff coefficient for quick flow", bounds = (0, 1), unit = "mm^(4/3)/d"]

ice_split_flux = @hydroflux begin
    Pice ~ fice * P
    Pnonice ~ (1 - fice) * P
end

pice_split_flux = @hydroflux begin
    Pices ~ step_func(T - t0) * Pice
    Picer ~ step_func(t0 - T) * Pice
end

pnonice_split_flux = @hydroflux begin
    Pnis ~ step_func(T - t0) * Pnonice
    Pnir ~ step_func(t0 - T) * Pnonice
end


bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Mis ~ min(His, asnow * max(0.0, T - tm))
        @hydroflux Mii ~ (1 - step_func(His)) * aice * max(0.0, T - t0)
        @hydroflux Pirs ~ step_func(His) * Picer
        @hydroflux Piri ~ Picer - Pirs
        @hydroflux Qis ~ ks * Sis
        @hydroflux Qii ~ ki * Sii
    end
    dfluxes = begin
        @stateflux His ~ Pices - Mis
        @stateflux Sis ~ Mis + Pirs - Qis
        @stateflux Sii ~ Mii + Piri - Qii
    end
end

bucket2 = @hydrobucket :bucket2 begin
    fluxes = begin
        @hydroflux Mnis ~ min(Hnis, asnow * max(0.0, T - tm))
        @hydroflux Peq ~ Mnis + Pnir
        @hydroflux Peff ~ Peq * max(0.0, Snis / a)^y
        @hydroflux Pinf ~ max(0.0, Peq - Peff)
        @hydroflux ET ~ Ep * max(0.0, Snis / a)^x
        @hydroflux Qsl ~ ksl * Snis
        @hydroflux Qqu ~ beta * max(0.0, Sniq)^(5 / 3)
    end
    dfluxes = begin
        @stateflux Hnis ~ Pnis - Mnis
        @stateflux Snis ~ Pinf - ET - Qsl
        @stateflux Sniq ~ Peff - Qqu
    end
end

flux1 = @hydroflux Qt ~ Qqu + Qsl + Qis + Qii

model = @hydromodel :gsmsocont begin
    ice_split_flux
    pice_split_flux
    pnonice_split_flux
    bucket1
    bucket2
    flux1
end

end