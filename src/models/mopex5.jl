module mopex5
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
@variables Ep [description = "Evapotranspiration rate", unit = "mm"]
@variables Epc [description = "Evapotranspiration capacity", unit = "mm"]
@variables Qn [description = "Snowmelt", unit = "mm/d"]
@variables dd [description = "surface runoff", unit = "mm/oC/d"]
@variables Ps [description = "snowfall", unit = "mm/d"]
@variables Pr [description = "precipitation as rain", unit = "mm/d"]
@variables ET1 [description = "evaporation", unit = "mm/d"]
@variables ET2 [description = "Evaporation", unit = "mm/d"]
@variables GSI [description = "growing season index based on parameters", unit = "-"]
@variables I [description = "Interception", unit = "mm/d"]
@variables Q1f [description = "Saturation excess flow", unit = "mm/d"]
@variables Qw [description = "deeper groundwater", unit = "mm/d"]
@variables Q2u [description = "slow runoff store", unit = "mm/d"]
@variables Q2f [description = "excess flow", unit = "mm/d"]
@variables Qf [description = "flow from fast flow routing reservoir", unit = "mm/d"]
@variables Qu [description = "flow from slow flow routing reservoir", unit = "mm/d"]
@variables Qt [description = "Total simulated flow", unit = "mm/d"]

# Model parameters
@parameters tcrit [description = "Snowfall & snowmelt temperature", bounds = (-3, 3), unit = "oC"]
@parameters ddf [description = "Degree-day factor for snowmelt", bounds = (0, 20), unit = "mm/oC/d"]
@parameters sb1 [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters tw [description = "Groundwater leakage time", bounds = (0, 1), unit = "1/d"]
@parameters i_alpha [description = "Intercepted fraction of Pr ", bounds = (0, 1), unit = "-"]
@parameters i_s [description = "Maximum Leaf Area Index timing ", bounds = (1, 365), unit = "d"]
@parameters tmin [description = "Growing Season Index minimum temperature ", bounds = (-10, 0), unit = "oC"]
@parameters trange [description = "Growing Season Index temperature range ", bounds = (1, 20), unit = "oC"]
@parameters tu [description = "Slow flow routing response time", bounds = (0, 1), unit = "1/d"]
@parameters se [description = "Root zone storage capacity as fraction of Sb2", bounds = (0.05, 0.95), unit = "-"]
@parameters sb2 [description = "Root zone storage capacity", bounds = (1, 2000), unit = "mm"]
@parameters tc [description = "Mean residence time", bounds = (0, 1), unit = "1/d"]

tmax = tmin + trange

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Ps ~ step_func(tcrit - T) * P
        @hydroflux Qn ~ ddf * min(0.0, T - tcrit)
        @hydroflux Epc ~ Ep * clamp((T - tmin) / (tmax - tmin), 0, 1)
    end
    dfluxes = begin
        @stateflux Sn ~ Ps - Qn
    end
end

bucket2 = @hydrobucket :bucket2 begin
    fluxes = begin
        @hydroflux Pr ~ step_func(T - tcrit) * P
        @hydroflux ET1 ~ S1 / sb1 * Epc
        # todo 加入时间？
        # @hydroflux I ~ max(0, i_alpha + (1 - i_alpha) * sin(2π * (t + I_s / (365 / d))))
        @hydroflux Q1f ~ step_func(S1 - sb1) * P
        @hydroflux Qw ~ tw * S1
    end
    dfluxes = begin
        @stateflux S1 ~ Pr - ET1 - Q1f - Qw
    end
end

bucket3 = @hydrobucket :bucket3 begin
    fluxes = begin
        @hydroflux ET2 ~ S2 / se * Epc
        @hydroflux Q2u ~ tu * S2
        @hydroflux Q2f ~ step_func(S2 - sb2) * Qw
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


model = @hydromodel :plateau begin
    bucket1
    bucket2
    bucket3
    bucket4
    bucket5
    flux1
end

end