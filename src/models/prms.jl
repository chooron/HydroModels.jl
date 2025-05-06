module prms
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables Sn [description = "thecurrent snow storage", unit = "mm"]
@variables Ps [description = "the rain that falls as snow", unit = "mm"]
@variables M [description = "the snowmelt", unit = "mm"]
@variables P [description = "precipitation", unit = "mm/d"]
@variables T [description = "temperature", unit = "oC"]
@variables XIN [description = " the current storage in the interception reservoir", unit = "mm"]
@variables Pin [description = " intercepted rainfall", unit = "mm/d"]
@variables Ein [description = "drained by evaporation", unit = "mm/d"]
@variables Ptf [description = "throughfall", unit = "mm/d"]
@variables Psm [description = "rainfall on non-impervious area", unit = "mm/d"]
@variables Pr [description = "rainfall", unit = "mm/d"]
@variables Ep [description = "the potential rate", unit = "mm/d"]
@variables Rstor
@variables Pim [description = "currentdepression storage", unit = "mm"]
@variables Mim [description = " given as the fraction 1 of snowmelt M", unit = ""]
@variables SAS [description = "surface runof", unit = "mm/d"]
@variables RECHR [description = " the current storage in the upper soil moisture zone", unit = "mm"]
@variables INF [description = "infiltration", unit = "mm/d"]
@variables Ea [description = "evaporation", unit = "mm/d"]
@variables Pc [description = "percolation", unit = "mm/d"]
@variables Msm [description = " incoming snowmelt ", unit = "mm/d"]
@variables Pby [description = " interception bypass ", unit = "mm/d"]
@variables SRO [description = " saturated area", unit = "mm/d"]
@variables Ei
@variables Eim [description = "evaporation", unit = "mm/d"]
@variables SMAV [description = "thecurrentstorage in the lower soil moisture zone", unit = "mm"]
@variables Et
@variables EXCS [description = " soil moisture excess ", unit = "mm/d"]
@variables RES [description = "the current storage in the runoff reservoir", unit = "mm"]
@variables QRES
@variables GAD [description = "groundwater drainage", unit = "mm/d"]
@variables RAS [description = "interflow component", unit = "mm/d"]
@variables SEP [description = "constant groundwater recharge", unit = "mm/d"]
@variables GW [description = "the current groundwater storage", unit = "mm"]
@variables BAS [description = "baseflow", unit = "mm/d"]
@variables SNK [description = "flow to deeper groundwater", unit = "mm/d"]
@variables Qt [description = "Total flow", unit = "mm/d"]

# Model parameters
@parameters tt [description = "Threshold temperature for snowfall and snowmelt", bounds = (-3, 5), unit = "oC"]
@parameters ddf [description = "Degree-day factor for snowmelt", bounds = (0, 20), unit = "mm/oC/d"]
@parameters alpha [description = "Fraction of precipitation on soil that is intercepted", bounds = (0, 1), unit = "-"]
@parameters beta [description = "Fraction of precipitation that falls on soil", bounds = (0, 1), unit = "-"]
@parameters stor [description = "Maximum interception storage", bounds = (0, 5), unit = "mm"]
@parameters retip [description = "Maximum depression storage", bounds = (0, 50), unit = "mm"]
@parameters fscn [description = "Minimum contributing area to surface runoff", bounds = (0, 1), unit = "-"]
@parameters scx [description = "Maximum contributing area to surface runoff", bounds = (0, 1), unit = "-"]
@parameters flz [description = "Maximum contributing area to surface runoff", bounds = (0.005, 0.995), unit = "-"]
@parameters stot [description = "Total soil moisture storage", bounds = (1.0, 2000.0), unit = "mm"]
@parameters cgw [description = "Maximum groundwater recharge rate", bounds = (0.0, 20.0), unit = "mm/d"]
@parameters resmax [description = "Maximum runoff routing storage", bounds = (1.0, 300.0), unit = "mm"]
@parameters k1 [description = "Runoff coefficient", bounds = (0.0, 1.0), unit = "1/d"]
@parameters k2 [description = "Runoff non-linearity", bounds = (1.0, 5.0), unit = "-"]
@parameters k3 [description = "Runoff coefficient", bounds = (0.0, 1.0), unit = "1/d"]
@parameters k4 [description = "Runoff coefficient", bounds = (0.0, 1.0), unit = "1/(mm*d)"]
@parameters k5 [description = "Runoff coefficient", bounds = (0.0, 1.0), unit = "1/d"]
@parameters k6 [description = "Runoff coefficient", bounds = (0.0, 1.0), unit = "1/d"]

scn = fscn * scx # Minimum contributing fraction area to saturation excess flow
remx = (1 - flz) * stot # Maximum upper soil moisture storage [mm]
smax = flz * stot # Maximum lower soil moisture storage [mm]

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Ps ~ step_func(tt - T) * P
        @hydroflux M ~ min(Sn, ddf * max(0.0, T - tt))
    end
    dfluxes = begin
        @stateflux Sn ~ Ps - M
    end
end

bucket2 = @hydrobucket :bucket2 begin
    fluxes = begin
        @hydroflux Pr ~ step_func(T - tt) * P
        @hydroflux Psm ~ beta * Pr
        @hydroflux Pin ~ alpha * Psm
        @hydroflux Ein ~ min(XIN, step_func(XIN) * beta * Ep)
        @hydroflux Ptf ~ min(XIN, step_func(XIN - stor) * Pin)
    end
    dfluxes = begin
        @stateflux XIN ~ Pin - Ein - Ptf
    end
end

bucket3 = @hydrobucket :bucket3 begin
    fluxes = begin
        @hydroflux Pim ~ (1 - beta) * Pr
        @hydroflux Mim ~ (1 - beta) * M
        @hydroflux Eim ~ min(Rstor, step_func(Rstor) * (1 - beta) * Ep)
        @hydroflux SAS ~ min(Rstor, step_func(Rstor - retip) * (Pim + Mim))
    end
    dfluxes = begin
        @stateflux Rstor ~ Pim + Mim - Eim - SAS
    end
end

bucket4 = @hydrobucket :bucket4 begin
    fluxes = begin
        @hydroflux Msm ~ beta * M
        @hydroflux Pby ~ (1 - alpha) * Psm
        @hydroflux SRO ~ (scn + ((scx - scn) * (RECHR / remx))) * (Msm + Ptf + Pby)
        @hydroflux INF ~ Msm + Ptf + Pby - SRO
        @hydroflux Ea ~ RECHR / remx * (Ep - Ein - Eim)
        @hydroflux Pc ~ step_func(RECHR - remx) * INF
    end
    dfluxes = begin
        @stateflux RECHR ~ INF - Ea - Pc
    end
end

bucket5 = @hydrobucket :bucket5 begin
    fluxes = begin
        @hydroflux Et ~ step_func(Ep - Ein - Eim - RECHR) * (SMAV / smax * (Ep - Ein - Eim - Ea))
        @hydroflux EXCS ~ step_func(SMAV - (smax - remx)) * Pc
        @hydroflux SEP ~ min(cgw, EXCS)
        @hydroflux QRES ~ min(EXCS - SEP, 0)
        @hydroflux GAD ~ k1 * ((RES / resmax)^k2)
        @hydroflux RAS ~ k3 * RES + k4 * (RES^2)
        @hydroflux BAS ~ k5 * GW
        @hydroflux SNK ~ k6 * GW
    end
    dfluxes = begin
        @stateflux SMAV ~ Pc - Et - EXCS
        @stateflux RES ~ QRES - GAD - RAS
        @stateflux GW ~ SEP + GAD - BAS - SNK
    end
end

flux1 = @hydroflux Qt ~ SAS + SRO + RAS + BAS

model = @hydromodel :prms begin
    bucket1
    bucket2
    bucket3
    bucket4
    bucket5
    flux1
end

end