module xinanjiang

using ..HydroModels
using ..HydroModels: step_func
# Model variables
@variables W [description = "current tension water storage", unit = "mm"]
@variables Pi [description = "infiltration", unit = "mm/d"]
@variables E [description = "evaporation", unit = "mm/d"]
@variables R [description = "runoff", unit = "mm/d"]
@variables P [description = "precipitation", unit = "mm/d"]
@variables Aim [description = "precipitation", unit = "-"]
@variables Wmax [description = "maximum tension water storage", unit = "mm"]
@variables Ep [description = "Evaporation occurs at potential rate", unit = "mm/d"]
@variables LM [description = "storage threshold", unit = "mm"]
@variables S [description = "current storage of free water", unit = "mm"]
@variables RS [description = "surface runoff", unit = "mm/d"]
@variables RI [description = "interflow", unit = "mm/d"]
@variables RG [description = "baseflow", unit = "mm/d"]
@variables Rprop [description = "proportion of surface runoff (1 - S/Smax)^Ex", unit = "-"]
@variables Smax [description = "maximum free water storage", unit = "mm"]
@variables SI [description = "current storage in the interflow routing reservoir", unit = "mm"]
@variables QI [description = "delayed interflow", unit = "mm/d"]
@variables SG [description = "current storage in the baseflow routing reservoir", unit = "mm"]
@variables QG [description = "delayed baseflow", unit = "mm/d"]
@variables RB [description = "direct rainoff", unit = "mm/d"]
@variables QS [description = "surface runoff", unit = "mm/d"]
@variables Qt [description = "Total runoff", unit = "mm/d"]
# Model parameters
@parameters Aim [description = "Fraction impervious area", bounds = (0, 1), unit = "-"]
@parameters a [description = "Contributing area curve inflection point", bounds = (-0.49, 0.49), unit = "-"]
@parameters b [description = "Contributing area curve shape parameter", bounds = (0, 10), unit = "-"]
@parameters Wmax [description = "Maximum tension water storage", bounds = (0.01, 0.99), unit = "-"]
@parameters LM [description = "Threshold for evaporation behaviour change", bounds = (0.01, 0.99), unit = "-"]
@parameters c [description = "Threshold and evaporation reduction factor", bounds = (0.01, 0.99), unit = "-"]
@parameters Smax [description = "Maximum free water storage", bounds = (1, 2000), unit = "mm"]
@parameters Ex [description = "Contributing area curve shape parameter", bounds = (0, 10), unit = "-"]
@parameters kI [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]
@parameters kG [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]
@parameters cI [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]
@parameters cG [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]

# Soil water component
bucket_1 = @hydrobucket :bucket_1 begin
    fluxes = begin
        @hydroflux Pi ~ (1 - Aim) * P
        @hydroflux R ~ max(0.0, step_func((W / Wmax) - (0.5 - a)) *
                                (Pi * (0.5 - a)^(1 - b) * max(0.0, (W / Wmax))^b) +
                                step_func((0.5 - a) - (W / Wmax)) *
                                (Pi * (1 - (0.5 + a)^(1 - b) * max(0.0, (1 - W / Wmax))^b)))
        @hydroflux E ~ clamp(W / LM, c, 1.0) * Ep
    end
    dfluxes = begin
        @stateflux W ~ Pi - E - R
    end
end

bucket_2 = @hydrobucket :bucket_2 begin
    fluxes = begin
        @hydroflux Rprop ~ max(0.0, 1 - S / Smax)^Ex
        @hydroflux RS ~ R * Rprop
        @hydroflux RI ~ kI * S * Rprop
        @hydroflux RG ~ kG * S * Rprop
    end

    dfluxes = begin
        @stateflux S ~ R - RS - RI - RG
    end
end


bucket_3 = @hydrobucket :bucket_3 begin
    fluxes = begin
        @hydroflux QI ~ cI * SI
    end
    dfluxes = begin
        @stateflux SI ~ RI - QI
    end
end

bucket_4 = @hydrobucket :bucket_4 begin
    fluxes = begin
        @hydroflux QG ~ cG * SG
        @hydroflux RB ~ Aim * P
        @hydroflux QS ~ RS + RB
        @hydroflux Qt ~ QS + QI + QG
    end

    dfluxes = begin
        @stateflux SG ~ RG - QG
    end
end

model = @hydromodel :xinanjiang begin
    bucket_1
    bucket_2
    bucket_3
    bucket_4
end

end
