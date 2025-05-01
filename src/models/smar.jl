
module smar
# Model variables
@variables sum_Sn [description = "total storage in the upper soil layer", unit = "mm"]
@variables S1 [description = "current storage in the upper soil layer", unit = "mm"]
@variables I [description = "infiltration", unit = "mm/d"]
@variables P [description = "rainfall", unit = "mm/d"]
@variables Ep [description = "Evaporation", unit = "mm/d"]
@variables P_star [description = "effective precipitation", unit = "mm/d"]
@variables R1 [description = "direct runoff", unit = "mm/d"]
@variables R2 [description = "infiltration excess runoff", unit = "mm/d"]
@variables E1 [description = "evaporation", unit = "mm/d"]
@variables q1 [description = "flow towards deeper soil layers", unit = "mm/d"]
@variables Ep_star [description = "Evaporation from this soil layer occurs at the effect potential rate", unit = "mm/d"]
@variables S2 [description = "current storage in the second soil layer", unit = "mm"]
@variables E2 [description = "the evaporation scaled by parameter", unit = "mm/d"]
@variables q2 [description = "overflow into the next layer", unit = "mm/d"]
@variables S3 [description = "current storage in the second soil layer", unit = "mm"]
@variables E3 [description = "the evaporation scaled by parameter", unit = "mm/d"]
@variables q3 [description = "overflow into the next layer", unit = "mm/d"]
@variables S4 [description = "current storage in the second soil layer", unit = "mm"]
@variables E4 [description = "evaporation scaled by parameter", unit = "mm/d"]
@variables q4 [description = "overflow into the next layer", unit = "mm/d"]
@variables S5 [description = "current storage in the second soil layer", unit = "mm"]
@variables E5 [description = "the evaporation scaled by parameter", unit = "mm/d"]
@variables R3 [description = "overflow towards groundwater", unit = "mm/d"]
@variables Gw [description = "the current groundwater storage", unit = "mm"]
@variables Qg [description = "groundwater flow", unit = "mm/d"]
@variables Qr [description = "flow routed per time step", unit = "mm/d"]

# Model parameters
@parameters H [description = "Fraction of effective precipitation that is direct runoff", bounds = (0,1), unit = "-"]
@parameters Y [description = "Infiltration rate", bounds = (0,200), unit = "mm/d"]
@parameters Smax [description = "Maximum soil moisture storage", bounds = (1,2000), unit = "mm"]
@parameters C [description = "Evaporation reduction parameter", bounds = (0,1), unit = "-"]
@parameters G [description = "Fraction of subsurface flow to groundwater", bounds = (0,1), unit = "-"]
@parameters KG [description = "Runoff coefficient", bounds = (0,1), unit = "d-1"]
@parameters N [description = "Gamma function parameter", bounds = (1,10), unit = "-"]
@parameters K [description = "Routing time parameter", bounds = (1,120), unit = "d"]

# Soil water component
soil_bucket_1 = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux P_star  ~ ifelse(P > Ep, P - Ep, 0)
        @hydroflux R1 ~ P_star * H * (sum_Sn / Smax)
        @hydroflux I ~ ifelse(P_star - R1 >= Y, Y, P_star - R1)
        @hydroflux R2 ~ (P_star - R1) - I
        @hydroflux Ep_star ~ ifelse(Ep > P, Ep - P, 0)
        @hydroflux E1 ~ Ep_star   
        @hydroflux q1 ~ ifelse(S1 >= Smax / 5, P_star - R1 - R2, 0)
    end

    dfluxes = begin
        @stateflux S1 ~ I - E1 - q1
    end
end

soil_bucket_2 = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux E2 ~ ifelse(S1 == 0, C * Ep, 0)
        @hydroflux q2 ~ ifelse(S2 >= Smax / 5, q1, 0)
    end

    dfluxes = begin
        @stateflux S2 ~ q1 - E2 - q2
    end
end

soil_bucket_3 = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux E3 ~ ifelse(S2 == 0, C^2 * Ep, 0)
        @hydroflux q3 ~ ifelse(S3 >= Smax / 5, q2, 0)
    end

    dfluxes = begin
        @stateflux S3 ~ q2 - E3 - q3
    end
end

soil_bucket_4 = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux E4 ~ ifelse(S3 == 0, C^3 * Ep, 0)
        @hydroflux q4 ~ ifelse(S4 >= Smax / 5, q3, 0)
    end

    dfluxes = begin
        @stateflux S4 ~ q3 - E4 - q4
    end
end

soil_bucket_5= @hydrobucket :soil begin
    fluxes = begin
        @hydroflux E5 ~ ifelse(S4 == 0, C^4 * Ep, 0)
        @hydroflux R3 ~ ifelse(S5 >= Smax / 5, q4, 0)
    end

    dfluxes = begin
        @stateflux S5 ~ q4 - E5 - R3
    end
end

soil_bucket_6= @hydrobucket :soil begin
    fluxes = begin
        @hydroflux Rg ~ G * R3
        @hydroflux Qg ~ KG * Gw
        @hydroflux Qr ~ (R1 + R2 + R3) * (1 / (K * gamma(N))) * (t / K)^(N - 1) * exp(-t / K)
        @hydroflux Qt ~ Qr + Qg
    end

    dfluxes = begin
        @stateflux Gw ~ Rg - Qg
    end
end

model = @hydromodel :smar begin
    bucket
end

end