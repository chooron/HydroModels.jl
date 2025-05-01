
module xinanjiang
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
@variables Smax [description = "maximum free water storage", unit = "mm"]
@variables SI [description = "current storage in the interflow routing reservoir", unit = "mm"]
@variables QI [description = "delayed interflow", unit = "mm/d"]
@variables SG [description = "current storage in the baseflow routing reservoir", unit = "mm"]
@variables QG [description = "delayed baseflow", unit = "mm/d"]
@variables RB [description = "direct rainoff", unit = "mm/d"]
@variables QS [description = "surface runoff", unit = "mm/d"]
@variables Qt [description = "Total runoff", unit = "mm/d"]
# Model parameters
@parameters Aim [description = "Fraction impervious area", bounds = (0,1), unit = "-"]
@parameters a [description = "Contributing area curve inflection point", bounds = (-0.49,0.49), unit = "-"]
@parameters b [description = "Contributing area curve shape parameter", bounds = (0,10), unit = "-"]
@parameters Wmax [description = "Maximum tension water storage", bounds = (0.01,0.99), unit = "-"]
@parameters LM [description = "Threshold for evaporation behaviour change", bounds = (0.01,0.99), unit = "-"]
@parameters c [description = "Threshold and evaporation reduction factor", bounds = (0.01,0.99), unit = "-"]
@parameters Smax [description = "Maximum free water storage", bounds = (1,2000), unit = "mm"]
@parameters Ex [description = "Contributing area curve shape parameter", bounds = (0,10), unit = "-"]
@parameters kI [description = "Runoff coefficient", bounds = (0,1), unit = "d-1"]
@parameters kG [description = "Runoff coefficient", bounds = (0,1), unit = "d-1"]
@parameters cI [description = "Runoff coefficient", bounds = (0,1), unit = "d-1"]
@parameters cG [description = "Runoff coefficient", bounds = (0,1), unit = "d-1"]

# Soil water component
soil_bucket_1 = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux Pi ~ (1-Aim)*P
        @hydroflux R ~ ifelse(W/Wmax <= 0.5 - a,
        Pi * (0.5 - a)^(1 - b) * (W / Wmax)^b,
        Pi * (1 - (0.5 + a)^(1 - b) * (1 - W / Wmax)^b)
    )
        @hydroflux E ~ ifelse(W > LM,
        Ep,
        ifelse(W >= c * LM,
            (W / LM) * Ep,
            c * Ep
        )
    )
  
    end
    dfluxes = begin
        @stateflux W ~ Pi-E-R
    end
end

soil_bucket_2 = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux RS ~ R * (1 - (1 - S / Smax)^Ex)
        @hydroflux RI ~ kI * S * (1 - (1 - S / Smax)^Ex)
        @hydroflux RG ~ kG * S * (1 - (1 - S / Smax)^Ex)
    end

    dfluxes = begin
        @stateflux S ~ R - RS - RI - RG
    end
end


soil_bucket_3 = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux QI ~ cI*SI
        
    end

    dfluxes = begin
        @stateflux SI ~ RI-QI
    end
end

soil_bucket_4 = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux QG ~ cG*SG
        @hydroflux Qt ~ QS+QI+QG
        @hydroflux QS ~ RS+RB
        @hydroflux RB ~ Aim*P
    end

    dfluxes = begin
        @stateflux SG ~ RG-QG
    end
end






model = @hydromodel :xinanjiang begin
    soil_bucket
end

end
