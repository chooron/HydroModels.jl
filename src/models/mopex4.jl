
module mopex4
# Model variables
@variables Sn [description = "current snow pack", unit = "mm"]
@variables P [description = "precipitation", unit = "mm/d"]
@variables Ps [description = "Precipitation occurs as snowfall", unit = "mm/d"]
@variables T [description = "current temperature", unit = "oC"]
@variables QN [description = "Snowmelt", unit = "mm/d"]
@variables S1 [description = "current storage in soil moisture", unit = "mm"]
@variables Pr [description = "precipitation as rain", unit = "mm/d"]
@variables ET1 [description = "Evaporation", unit = "mm/d"]
@variables Ep [description = "potential evapotranspiration", unit = "mm/d"]
@variables I [description = "Interception", unit = "mm/d"]
@variables tmax [description = "the length of the seasonal cycle", unit = "d"]
@variables Q1f [description = "Saturation excess flow", unit = "mm/d"]
@variables Qw [description = "infiltration to deeper groundwater", unit = "mm/d"]
@variables S2 [description = "current groundwater storage", unit = "mm"]
@variables ET2 [description = "Evaporation", unit = "mm/d"]
@variables Q2u [description = "Leakage to the slow runoff store", unit = "mm/d"]
@variables Q2f [description = "excess flow", unit = "mm/d"]
@variables Sc1 [description = "current storage in the fast flow routing reservoir", unit = "mm"]
@variables Qf [description = "Routed flow", unit = "mm/d"]
@variables Sc2 [description = "current storage in the slow flow routing reservoir", unit = "mm"]
@variables Qu [description = "Routed flow", unit = "mm/d"]
@variables Qt [description = "Total simulated flow", unit = "mm/d"]
# Model parameters
@parameters Tcrit [description = "Threshold temperature for snowfall and melt", bounds = (-3,3), unit = "oC"]
@parameters ddf [description = "Degree-day factor", bounds = (0,20), unit = "mm/oC/d"]
@parameters Sb1 [description = "Maximum soil moisture storage", bounds = (1,2000), unit = "mm"]
@parameters t [description = "Runoff coefficient", bounds = (0,1), unit = "d-1"]
@parameters tw [description = "Runoff coefficient", bounds = (0,1), unit = "d-1"]
@parameters I_alpha [description = "Mean intercepted fraction of precipitation", bounds = (0,1), unit = "-"]
@parameters I_s [description = "Timing of peak interception capacity", bounds = (1,365), unit = "d"]
@parameters Sb2 [description = "Maximum deep storage", bounds = (1,2000), unit = "mm"]
@parameters tu [description = "Runoff coefficient", bounds = (0,1), unit = "d-1"]
@parameters Se [description = "Maximum groundwater storage capacity", bounds = (0.05,0.95), unit = "mm"]
@parameters tc [description = "Runoff coefficient", bounds = (0,1), unit = "d-1"]

soil_bucket_1 = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux Ps ~ P*step_func(Tcrit-T)
        @hydroflux QN ~ ddf*max(T-Tcrit,0.00)

    end

    dfluxes = begin
        @stateflux Sn ~ Ps-QN
    end
end

soil_bucket_2 = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux Pr ~ ifelse(T > Tcrit, P, 0)
        @hydroflux ET1 ~ (S1 / Sb1) * Ep
        @hydroflux I ~ max(0, I_alpha + (1 - I_alpha) * cos(2π * (t - I_s) / tmax)) * Pr
        @hydroflux Q1f ~ ifelse(S1 >= Sb1, P, 0)
        @hydroflux Qw ~ tw * S1
    end

    dfluxes = begin
        @stateflux S1 ~ Pr - ET1 - I - Q1f - Qw
    end
end

soil_bucket_3 = @hydrobucket :soil2 begin
    fluxes = begin
        @hydroflux ET2 ~ (S2 / Se) * Ep
        @hydroflux Q2u ~ tu * S2
        @hydroflux Q2f ~ ifelse(S2 >= Sb2, Qw, 0)
    end

    dfluxes = begin
        @stateflux S2 ~ Qw - ET2 - Q2u - Q2f
    end
end

soil_bucket_4 = @hydrobucket :soil2 begin
    fluxes = begin
        @hydroflux Qf ~ tc*Sc1
      
    end

    dfluxes = begin
        @stateflux Sc1 ~ Q1f+Q2f-Qf
    end
end

soil_bucket_5 = @hydrobucket :soil2 begin
    fluxes = begin
        @hydroflux Qu ~ tc*Sc2
        @hydroflux Qt ~ Qf+Qu
      
    end

    dfluxes = begin
        @stateflux Sc2 ~ Q2u-Qu
    end
end

model = @hydromodel :mopex4 begin
    soil_bucket
end

end