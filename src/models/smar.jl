module smar
using ..HydroModels
using ..HydroModels: step_func
# Model variables
@variables S1 [description = "current storage in the upper soil layer", unit = "mm"]
@variables I [description = "infiltration", unit = "mm/d"]
@variables P [description = "rainfall", unit = "mm/d"]
@variables Ep [description = "Evaporation", unit = "mm/d"]
@variables Pe [description = "effective precipitation", unit = "mm/d"]
@variables R1 [description = "direct runoff", unit = "mm/d"]
@variables R2 [description = "infiltration excess runoff", unit = "mm/d"]
@variables E1 [description = "evaporation", unit = "mm/d"]
@variables q1 [description = "flow towards deeper soil layers", unit = "mm/d"]
@variables EPe [description = "Evaporation from this soil layer occurs at the effect potential rate", unit = "mm/d"]
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
@variables Rg [description = "groundwater recharge", unit = "mm/d"]
@variables Qg [description = "groundwater flow", unit = "mm/d"]
@variables Qr [description = "flow routed per time step", unit = "mm/d"]
@variables Qt [description = "total flow", unit = "mm/d"]

# Model parameters
@parameters H [description = "Fraction of effective precipitation that is direct runoff", bounds = (0, 1), unit = "-"]
@parameters Y [description = "Infiltration rate", bounds = (0, 200), unit = "mm/d"]
@parameters Smax [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters C [description = "Evaporation reduction parameter", bounds = (0, 1), unit = "-"]
@parameters G [description = "Fraction of subsurface flow to groundwater", bounds = (0, 1), unit = "-"]
@parameters KG [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]
@parameters N [description = "Gamma function parameter", bounds = (1, 10), unit = "-"]
@parameters K [description = "Routing time parameter", bounds = (1, 120), unit = "d"]

# Soil water component
bucket_1 = @hydrobucket :bucket_1 begin
    fluxes = begin
        @hydroflux Pe ~ max(P - Ep, 0)
        @hydroflux EPe ~ max(Ep - P, 0)
        @hydroflux R1 ~ Pe * H * ((S1 + S2 + S3 + S4 + S5) / Smax)
        @hydroflux I ~ min(Y, Pe - R1)
        @hydroflux R2 ~ (Pe - R1) - I
        @hydroflux E1 ~ EPe
        @hydroflux q1 ~ step_func(S1 - Smax / 5) * (Pe - R1 - R2)

        @hydroflux E2 ~ min(S2, C * Ep * (1- step_func(S1)))
        @hydroflux q2 ~ step_func(S2 - Smax / 5) * q1

        @hydroflux E3 ~ min(S3, C^2 * Ep * (1- step_func(S2)))
        @hydroflux q3 ~ step_func(S3 - Smax / 5) * q2

        @hydroflux E4 ~ min(S4, C^3 * Ep * (1- step_func(S3)))
        @hydroflux q4 ~ step_func(S4 - Smax / 5) * q3

        @hydroflux E5 ~ min(S5, C^4 * Ep * (1- step_func(S4)))
        @hydroflux R3 ~ step_func(S5 - Smax / 5) * q4
    end

    dfluxes = begin
        @stateflux S1 ~ I - E1 - q1
        @stateflux S2 ~ q1 - E2 - q2
        @stateflux S3 ~ q2 - E3 - q3
        @stateflux S4 ~ q3 - E4 - q4
        @stateflux S5 ~ q4 - E5 - R3
    end
end

bucket_2 = @hydrobucket :bucket_2 begin
    fluxes = begin
        @hydroflux Rg ~ G * R3
        @hydroflux Qg ~ KG * Gw
        # todo: add routing
        @hydroflux Qr ~ (R1 + R2 + R3) #  * (1 / (K * gamma(N))) * (t / K)^(N - 1) * exp(-t / K)
        @hydroflux Qt ~ Qr + Qg
    end

    dfluxes = begin
        @stateflux Gw ~ Rg - Qg
    end
end

model = @hydromodel :smar begin
    bucket_1
    bucket_2
end

end