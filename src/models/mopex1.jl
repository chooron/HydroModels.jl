module mopex1
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables S1 [description = "current storage in soil moisture", unit = "mm"]
@variables P [description = "precipitation", unit = "mm/d"]
@variables ET1 [description = "Evaporation", unit = "mm/d"]
@variables Ep [description = "potential evapotransporation", unit = "mm/d"]
@variables Q1f [description = "Saturation excess flow", unit = "mm/d"]
@variables Qw [description = "Infiltration to deeper groundwater", unit = "mm/d"]
@variables S2 [description = "current groundwater storage", unit = "mm"]
@variables ET2 [description = "Evaporation", unit = "mm/d"]
@variables Se [description = "groundwater storage capacity", unit = "mm"]
@variables Q2u [description = "Leakage to the slow runoff store", unit = "mm/d"]
@variables Sc1 [description = "current storage in the fast flow routing reservoir", unit = "mm"]
@variables Qf [description = "Routed flow", unit = "mm/d"]
@variables Sc2 [description = "current storage in the slow flow routing reservoir,", unit = "mm"]
@variables Qu [description = "Routed flow", unit = "mm/d"]
@variables Qt [description = "Total simulated flow", unit = "mm/d"]

# Model parameters
@parameters Sb1 [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters tw [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]
@parameters tu [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]
@parameters Se [description = "Maximum groundwater storage capacity", bounds = (1, 2000), unit = "mm"]
@parameters tc [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]

# Soil water component
bucket_1 = @hydrobucket :bucket_1 begin
    fluxes = begin
        @hydroflux ET1 ~ (S1 / Sb1) * Ep
        @hydroflux Q1f ~ step_func(S1 - Sb1) * P
        @hydroflux Qw ~ tw * S1
    end

    dfluxes = begin
        @stateflux S1 ~ P - ET1 - Q1f - Qw
    end
end

bucket_2 = @hydrobucket :bucket_2 begin
    fluxes = begin
        @hydroflux ET2 ~ (S2 / Se) * Ep
        @hydroflux Q2u ~ tu * S2
    end

    dfluxes = begin
        @stateflux S2 ~ Qw - ET2 - Q2u
    end
end

bucket_3 = @hydrobucket :bucket_3 begin
    fluxes = begin
        @hydroflux Qf ~ tc * Sc1
    end

    dfluxes = begin
        @stateflux Sc1 ~ Q1f - Qf
    end
end

bucket_4 = @hydrobucket :bucket_4 begin
    fluxes = begin
        @hydroflux Qu ~ tc * Sc2
        @hydroflux Qt ~ Qf + Qu
    end

    dfluxes = begin
        @stateflux Sc2 ~ Q2u - Qu
    end
end

model = @hydromodel :mopex1 begin
    bucket_1
    bucket_2
    bucket_3
    bucket_4
end

end