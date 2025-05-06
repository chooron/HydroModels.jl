module modhydrolog

using ..HydroModels
using ..HydroModels: step_func
# Model variables
@variables I [description = "current interception storage", unit = "mm"]
@variables P [description = "rainfall", unit = "mm/d"]
@variables Ei [description = "evaporation from the interception store", unit = "mm/d"]
@variables EXC [description = "excess rainfall", unit = "mm/d"]
@variables Ep [description = "Evaporation", unit = "mm/d"]
@variables SMS [description = "current storage in the soil moisture store", unit = "mm"]
@variables SMSProp [description = "current storage in the soil moisture store as a proportion of the maximum", unit = "-"]
@variables SMF [description = "infiltration", unit = "mm/d"]
@variables DINF [description = "delayed infiltration", unit = "mm/d"]
@variables INF [description = "total infiltration", unit = "mm/d"]
@variables INT [description = "interflow and saturation excess flow", unit = "mm/d"]
@variables REC [description = "preferential recharge of groundwater", unit = "mm/d"]
@variables ET [description = "evaporation from the soil moisture that occurs at the potential rate when possible", unit = "mm/d"]
@variables GWF [description = "flow to the groundwater store", unit = "mm/d"]
@variables RUN [description = "surface runoff", unit = "mm/d"]
@variables SRUN [description = "surface runoff", unit = "mm/d"]
@variables RATE [description = "rate of infiltration", unit = "mm/d"]
@variables TRAP [description = "the part of overland flow captured in the depression store", unit = "mm/d"]
@variables ED [description = "the evaporation from the depression store", unit = "mm/d"]
@variables SEEP [description = "the exchange with a deeper aquifer", unit = "mm/d"]
@variables FLOW [description = "the exchange with the channel", unit = "mm/d"]
@variables Qt [description = "total runoff", unit = "mm/d"]
@variables DH D GW CH

# Model parameters
@parameters INSC [description = "Maximum interception capacity", bounds = (0, 5), unit = "mm"]
@parameters COEFF [description = "Maximum infiltration loss", bounds = (0, 600), unit = "mm/d"]
@parameters SQ [description = "Infiltration loss exponent", bounds = (0, 15), unit = "-"]
@parameters SMSC [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters SUB [description = "Proportionality constant", bounds = (0, 1), unit = "-"]
@parameters CRAK [description = "Proportionality constant", bounds = (0, 1), unit = "-"]
@parameters EM [description = "Proportionality constant", bounds = (0, 20), unit = "mm/d"]
@parameters DSC [description = "Maximum depression storage", bounds = (0, 50), unit = "mm"]
@parameters ADS [description = "Fraction of area functioning as depression store", bounds = (0, 1), unit = "-"]
@parameters MD [description = "Depression store shape parameter", bounds = (0.99, 1), unit = "-"]
@parameters VCOND [description = "Runoff coefficient", bounds = (0, 0.5), unit = "mm/d"]
@parameters DLEV [description = "Datum of groundwater store", bounds = (-10, 10), unit = "mm"]
@parameters k1 [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]
@parameters k2 [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]
@parameters k3 [description = "Runoff coefficient", bounds = (0, 100), unit = "d-1"]

# Soil water component
bucket_1 = @hydrobucket :bucket_1 begin
    fluxes = begin
        @hydroflux Ei ~ min(I, Ep)
        @hydroflux EXC ~ P * step_func(I - INSC)
    end

    dfluxes = begin
        @stateflux I ~ P - Ei - EXC
    end
end

bucket_2 = @hydrobucket :bucket_2 begin
    fluxes = begin
        @hydroflux SMSProp ~ SMS / SMSC
        @hydroflux INF ~ min(COEFF * exp(-SQ * SMSProp), EXC)
        @hydroflux INT ~ SUB * SMSProp * INF
        @hydroflux REC ~ CRAK * SMSProp * (INF - INT)
        @hydroflux SMF ~ INF - INT - REC
        @hydroflux ET ~ min(EM * SMSProp, Ep)
        @hydroflux GWF ~ step_func(SMS - SMSC) * SMF

        @hydroflux RUN ~ EXC - INF
        @hydroflux RATE ~ COEFF * exp(-SQ * SMSProp) - INF - INT - REC
        @hydroflux TRAP ~ ADS * exp(-MD * D / (DSC - D)) * RUN
        @hydroflux ED ~ step_func(D) * ADS * Ep
        @hydroflux DINF ~ step_func(D) * ADS * RATE
    end

    dfluxes = begin
        @stateflux SMS ~ SMF + DINF - ET - GWF
        @stateflux D ~ TRAP - ED - DINF
    end
end

sigmoid(x) = 1 / (1 + exp(-x))

bucket_3 = @hydrobucket :bucket_3 begin
    fluxes = begin
        @hydroflux SEEP ~ VCOND * (GW - DLEV)
        @hydroflux SRUN ~ RUN - TRAP
        @hydroflux FLOW ~ max(-SRUN, sigmoid(GW) * (k1 * GW * sigmoid(GW) + k2 * (1 - exp(-k3 * GW * sigmoid(GW)))))
        @hydroflux Qt ~ max(CH, 0)
    end

    dfluxes = begin
        @stateflux GW ~ REC + GWF - SEEP - FLOW
        @stateflux CH ~ SRUN + INT + FLOW - Qt
    end
end


model = @hydromodel :modhydrolog begin
    bucket_1
    bucket_2
    bucket_3
end

end