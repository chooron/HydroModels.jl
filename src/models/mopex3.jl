module mopex3
using ..HydroModels
using ..HydroModels: step_func

# Model variables

@variables Sn [description = "current snow pack", unit = "mm"]
@variables S1 [description = "current storage in soil moisture", unit = "mm"]
@variables S2 [description = "current groundwater storage", unit = "mm"]
@variables Sc1 [description = "current storage in the fast flow routing reservoir", unit = "mm"]
@variables Sc2 [description = "current storage in the slow flow routing reservoir", unit = "mm"]

@variables P [description = "precipitation", unit = "mm/d"]
@variables T [description = "current temperature", unit = "oC"]
@variables Ep [description = "potential evapotranspiration", unit = "mm/d"]

@variables Ps [description = "snowfall", unit = "mm/d"]
@variables Qn [description = "Snowmelt", unit = "mm/d"]
@variables dd [description = "surface runoff", unit = "mm/oC/d"]
@variables Pr [description = "precipitation", unit = "mm/d"]
@variables ET1 [description = "evaporation", unit = "mm/d"]
@variables Q1f [description = "Saturation excess flow", unit = "mm/d"]
@variables Qw [description = "deeper groundwater", unit = "mm/d"]

@variables ET2 [description = "evaporation", unit = "mm/d"]
@variables Q2u [description = "slow runoff store", unit = "mm/d"]
@variables Q2f [description = "excess flow", unit = "mm/d"]
@variables Qu [description = "Routed flow", unit = "mm/d"]
@variables Qf [description = "Routed flow", unit = "mm/d"]
@variables Qt [description = "Total simulated flow", unit = "mm/d"]

# Model parameters
@parameters tcrit [description = "Snowfall & snowmelt temperature", bounds = (-3, 3), unit = "oC"]
@parameters ddf [description = "Degree-day factor for snowmelt", bounds = (0, 20), unit = "mm/oC/d"]
@parameters Smax1 [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters tw [description = "Groundwater leakage time", bounds = (0, 1), unit = "1/d"]
@parameters tu [description = "Slow flow routing response time", bounds = (0, 1), unit = "1/d"]
@parameters se [description = "Root zone storage capacity as fraction of Smax2", bounds = (0.05, 0.95), unit = "-"]
@parameters Smax2 [description = "Root zone storage capacity", bounds = (1, 2000), unit = "mm"]
@parameters tc [description = "Mean residence time", bounds = (0, 1), unit = "1/d"]

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Ps ~ step_func(tcrit - T) * P
        @hydroflux Qn ~ min(Sn, ddf * max(0.0, T - tcrit))
    end
    dfluxes = begin
        @stateflux Sn ~ Ps - Qn
    end
end

bucket2 = @hydrobucket :bucket2 begin
    fluxes = begin
        @hydroflux Pr ~ step_func(T - tcrit) * P
        @hydroflux ET1 ~ (S1 / Smax1) * Ep
        @hydroflux Q1f ~ step_func(S1 - Smax1) * (Pr + Qn)
        @hydroflux Qw ~ tw * S1
    end
    dfluxes = begin
        @stateflux S1 ~ Pr - ET1 - Q1f - Qw
    end
end

bucket3 = @hydrobucket :bucket3 begin
    fluxes = begin
        @hydroflux ET2 ~ (S2 / (se * Smax2)) * Ep
        @hydroflux Q2u ~ tu * S2
        @hydroflux Q2f ~ step_func(S2 - Smax2) * Qw
    end
    dfluxes = begin
        @stateflux S2 ~ Qw - ET2 - Q2u - Q2f
    end
end

bucket4 = @hydrobucket :bucket4 begin
    fluxes = begin
        @hydroflux Qf ~ tc * Sc1
    end
    dfluxes = begin
        @stateflux Sc1 ~ Q1f + Q2f - Qf
    end
end

bucket5 = @hydrobucket :bucket5 begin
    fluxes = begin
        @hydroflux Qu ~ tc * Sc2
    end
    dfluxes = begin
        @stateflux Sc2 ~ Q2u - Qu
    end
end

flux1 = @hydroflux Qt ~ Qf + Qu

model = @hydromodel :mopex3 begin
    bucket1
    bucket2
    bucket3
    bucket4
    bucket5
    flux1
end

end