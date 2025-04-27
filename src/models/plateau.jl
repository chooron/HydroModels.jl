module plateau
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables P [description = "Incoming precipitation", unit = "mm/d"]
@variables Ep [description = "potential rate", unit = "mm/d"]

@variables Su [description = "above the wilting point Swp", unit = "mm"]
@variables R [description = "Storage excess", unit = "mm"]

@variables Pe [description = "further divided into infiltration", unit = "mm/d"]
@variables Pi [description = "based on the maximum infiltration rate Fmax", unit = "mm/d"]
@variables Pie [description = "infiltration excess", unit = "mm/d"]
@variables Pie_routed [description = "infiltration excess routed", unit = "mm/d"]
@variables Cap [description = "Capillary rise", unit = "mm/d"]
@variables Et [description = "Evaporation from soil moisture", unit = "mm/d"]
@variables Spgw [description = "current groundwater storage", unit = "mm"]
@variables Qpgw [description = "Groundwater flow ", unit = "mm/d"]
@variables Qt [description = "Total runoff", unit = "mm/d"]
@variables t
# Model parameters
@parameters Fmax [description = "Maximum infiltration rate", bounds = (0, 200), unit = "mm/d"]
@parameters Dp [description = "Interception evaporation", bounds = (0, 5), unit = "mm/d"]
@parameters Sumax [description = "Maximum soil moisture storage", bounds = (1, 200), unit = "mm"]
@parameters lp [description = "Wilting point as fraction of Smax", bounds = (0.05, 0.95), unit = "-"]
@parameters erf [description = "Evaporation reduction factor", bounds = (0, 1), unit = "-"]
@parameters Tp [description = "Unit Hydrograph time base", bounds = (1, 120), unit = "d"]
@parameters C [description = "Capillary rise rate", bounds = (0, 4), unit = "mm/d"]
@parameters Kp [description = "Runoff coefficient", bounds = (0, 1), unit = "1/d"]

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Pe ~ max(P - Dp, 0)
        @hydroflux Pi ~ min(Pe, Fmax)
        @hydroflux Pie ~ Pe - Pi
        @hydroflux Et ~ min(Su, Ep * max(erf * (Su - lp * Sumax) / (Sumax - lp * Sumax), 0))
        @hydroflux R ~ (Pi + C) * step_func(Su - Sumax)
        @hydroflux Qpgw ~ Kp * Spgw
    end
    dfluxes = begin
        @stateflux Su ~ Pi + C - Et - R
        @stateflux Spgw ~ R - C - Qpgw
    end
end

uh_1 = @unithydro begin
    uh_func = begin
        ceil(Tp) => 1 / (0.5 * Tp^2) * (0.5 * Tp^2 - 0.5 * (t - 1)^2)
        Tp => 1 / (0.5 * Tp^2) * (0.5 * t^2 - 0.5 * (t - 1)^2)
    end
    uh_vars = [Pie]
    configs = (solvetype=:SPARSE, suffix=:_routed)
end

flux1 = @hydroflux Qt ~ Pie_routed + Qpgw

model = @hydromodel :plateau begin
    bucket1
    uh_1
    flux1
end

end