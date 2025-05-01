module modhydrolog
# Model variables
@variables I [description = "current interception storage", unit = "mm"]
@variables P [description = "rainfall", unit = "mm/d"]
@variables Ei [description = "evaporation from the interception store", unit = "mm/d"]
@variables EXC [description = "excess rainfall", unit = "mm/d"]
@variables Ep [description = "Evaporation", unit = "mm/d"]
@variables SMS [description = "current storage in the soil moisture store", unit = "mm"]
@variables SMF [description = "infiltration", unit = "mm/d"]
@variables DINF [description = "delayed infiltration", unit = "mm/d"]
@variables INF [description = "total infiltration", unit = "mm/d"]
@variables INT [description = "interflow and saturation excess flow", unit = "mm/d"]
@variables REC [description = "preferential recharge of groundwater", unit = "mm/d"]
@variables ET [description = "evaporation from the soil moisture that occurs at the potential rate when possible", unit = "mm/d"]
@variables GWF [description = "flow to the groundwater store", unit = "mm/d"] 
@variables TRAP [description = "the part of overland flow captured in the depression store", unit = "mm/d"] 
@variables ED [description = "the evaporation from the depression store", unit = "mm/d"] 
@variables SEEP [description = "the exchange with a deeper aquifer", unit = "mm/d"]
@variables FLOW [description = "the exchange with the channel", unit = "mm/d"]
@variables Qt [description = "total runoff", unit = "mm/d"]

# Model parameters
@parameters INSC [description = "Maximum interception capacity", bounds = (0,5), unit = "mm"]
@parameters COEFF [description = "Maximum infiltration loss", bounds = (0,600), unit = "mm/d"]
@parameters SQ [description = "Infiltration loss exponent", bounds = (0,15), unit = "-"]
@parameters SMSC [description = "Maximum soil moisture storage", bounds = (1,2000), unit = "mm"]
@parameters SUB [description = "Proportionality constant", bounds = (0,1), unit = "-"]
@parameters CRAK [description = "Proportionality constant", bounds = (0,1), unit = "-"]
@parameters EM [description = "Proportionality constant", bounds = (0,20), unit = "mm/d"]
@parameters DSC [description = "Maximum depression storage", bounds = (0,50), unit = "mm"]
@parameters ADS [description = "Fraction of area functioning as depression store", bounds = (0,1), unit = "-"]
@parameters MD [description = "Depression store shape parameter", bounds = (0.99,1), unit = "-"]
@parameters VCOND [description = "Runoff coefficient", bounds = (0,0.5), unit = "mm/d"]
@parameters DLEV [description = "Datum of groundwater store", bounds = (-10,10), unit = "mm"]
@parameters k1 [description = "Runoff coefficient", bounds = (0,1), unit = "d-1"]
@parameters k2 [description = "Runoff coefficient", bounds = (0,1), unit = "d-1"]
@parameters k3 [description = "Runoff coefficient", bounds = (0,100), unit = "d-1"]

# Soil water component
soil_bucket_1 = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux Ei ~ ifelse(I > 0, Ep, 0)
        @hydroflux EXC ~ ifelse(I == INSC, P, 0)
    end

    dfluxes = begin
        @stateflux I ~ P - Ei - EXC
    end
end

soil_bucket_2 = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux INF ~ min(
            COEFF * exp(-SQ * SMS / SMSC),
            EXC
        )
        @hydroflux INT ~ SUB * (SMS / SMSC) * INF
        @hydroflux REC ~ CRAK * (SMS / SMSC) * (INF - INT)
        @hydroflux SMF ~ INF - INT - REC
        @hydroflux ET ~ min(
            EM * (SMS / SMSC),
            PET
        )
        @hydroflux GWF ~ ifelse(SMS == SMSC, SMF, 0)
    end

    dfluxes = begin
        @stateflux SMS ~ SMF + DINF - ET - GWF
    end
end

soil_bucket_3= @hydrobucket :soil begin
    fluxes = begin
        @hydroflux RUN ~ EXC - INF
        @hydroflux RATE ~ COEFF * exp(-SQ * SMS / SMSC) - INF - INT - REC
        @hydroflux TRAP ~ ADS * exp(-MD * D / (DSC - D)) * RUN
        @hydroflux ED ~ ifelse(D > 0, ADS * Ep, 0)
        @hydroflux DINF ~ ifelse(D > 0, ADS * RATE, 0)
    end

    dfluxes = begin
        @stateflux D ~ TRAP - ED - DINF
    end
end

soil_bucket_4 = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux SEEP ~ VCOND * (GW - DLEV)

        @hydroflux FLOW ~ ifelse(GW >= 0,
            k1 * abs(GW) + k2 * (1 - exp(-k3 * abs(GW))),
            - (k1 * abs(GW) + k2 * (1 - exp(-k3 * abs(GW))))
        )
    end

    dfluxes = begin
        @stateflux GW ~ REC + GWF - SEEP - FLOW
    end
end

soil_bucket_5= @hydrobucket :soil begin
    fluxes = begin
        @hydroflux SRUN ~ RUN - TRAP
        @hydroflux Q ~ ifelse(CH > 0, CH, 0)
    end

    dfluxes = begin
        @stateflux CH ~ SRUN + INT + FLOW - Q
    end
end

model = @hydromodel :modhydrolog begin
    soil_bucket
end

end