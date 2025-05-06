module echo
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables I [description = "current interception storage", unit = "mm"]
@variables T [description = "temperature", unit = "oC"]
@variables P [description = "precipitation", unit = "mm/d"]
@variables Ei [description = "evaporation", unit = "mm/d"]
@variables Pn [description = "net precipitation", unit = "mm/d"]
@variables Ep [description = "the potential rate", unit = "mm/d"]
@variables Hs [description = "current storage in the snow pack", unit = "mm"]
@variables Ps [description = "precipitation-as-snow", unit = "mm/d"]
@variables Fs [description = "refreezing of melted snow", unit = "mm/d"]
@variables Ms [description = "snowmelt", unit = "mm/d"]
@variables Gs [description = "ground-heat flux", unit = "mm/d"]
@variables Hw [description = "current storage of liquid water in the snow pack", unit = "mm"]
@variables Pr [description = "precipitation-as-rain", unit = "mm/d"]
@variables Mw [description = "outflow of melt water", unit = "mm/d"]
@variables S [description = "current storage in the soil moisture zone", unit = "mm"]
@variables Fi [description = "infiltration", unit = "mm/d"]
@variables RD [description = "Dunne-type runoff", unit = "mm/d"]
@variables Et_pot [description = "potential evapotranspiration", unit = "mm/d"]
@variables Et [description = "evapotranspiration", unit = "mm/d"]
@variables L [description = "leakage", unit = "mm/d"]
@variables Peq [description = "equivalent precipitation", unit = "mm/d"]
@variables RH [description = "Horton-type runoff", unit = "mm/d"]
@variables Sfast [description = "current storage in the fast runoff reservoir", unit = "mm"]
@variables Lf [description = "leakage-to-fast-flow", unit = "mm/d"]
@variables Rf [description = "fast runoff", unit = "mm/d"]
@variables Ls [description = "leakage-to-slow-flow", unit = "mm/d"]
@variables Lmax [description = "maximum leakage rate", unit = "mm/d"]
@variables Sslow [description = "current storage in the slow runoff reservoir", unit = "mm"]
@variables Ls [description = "leakage-to-slow-flow", unit = "mm/d"]
@variables Rs [description = "slow runoff", unit = "mm/d"]
@variables Qt [description = "Total flow", unit = "mm/d"]

# Model parameters
@parameters rho [description = "Maximum interception storage", bounds = (0, 5), unit = "mm"]
@parameters Ts [description = "Threshold temperature for snowfall", bounds = (-3, 5), unit = "oC"]
@parameters Tm [description = "Threshold temperature for snowmelt", bounds = (-3, 3), unit = "oC"]
@parameters an [description = "Degree-day factor for snowmelt", bounds = (0, 20), unit = "mm/oC/d"]
@parameters af [description = "Degree-day factor reduction factor for refreezing", bounds = (0, 1), unit = "-"]
@parameters Gmax [description = "Snow melt through ground heat flux rate", bounds = (0, 2), unit = "mm/d"]
@parameters theta [description = "Water holding capacity as fraction of current snow pack", bounds = (0, 1), unit = "-"]
@parameters phi [description = "Maximum Horton type flow rate", bounds = (0, 200), unit = "mm/d"]
@parameters Smax [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters sw [description = "Wilting point", bounds = (0.05, 0.95), unit = "-"]
@parameters sm [description = "Plant stress point", bounds = (0.05, 0.95), unit = "-"]
@parameters Ksat [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]
@parameters c [description = "Runoff non-linearity", bounds = (0, 5), unit = "-"]
@parameters Lmax [description = "Maximum leakage rate", bounds = (0, 20), unit = "mm/d"]
@parameters kf [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]
@parameters ks [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]

# Soil water component
bucket_1 = @hydrobucket :bucket_1 begin
    fluxes = begin
        @hydroflux Ei ~ step_func(I) * Ep
        @hydroflux Pn ~ step_func(I - rho) * P
    end

    dfluxes = begin
        @stateflux I ~ P - Ei - Pn
    end
end

bucket_2 = @hydrobucket :bucket_2 begin
    fluxes = begin
        @hydroflux Ps ~ Pn * step_func(Ts - T)
        @hydroflux Pr ~ Pn * step_func(T - Ts)
        @hydroflux Ms ~ clamp(an * (T - Tm), 0.0, Hs)
        @hydroflux Mw ~ step_func(Hw - theta * Hs) * (Pr + Ms)
        @hydroflux Fs ~ clamp(af * an * (Tm - T), 0.0, Hw)
        @hydroflux Gs ~ min(Hs, Gmax)
    end

    dfluxes = begin
        @stateflux Hs ~ Ps + Fs - Ms - Gs
        @stateflux Hw ~ Pr + Ms - Fs - Mw
    end
end

bucket_3 = @hydrobucket :bucket_3 begin
    fluxes = begin
        @hydroflux Peq ~ Mw + Gs
        @hydroflux RH ~ step_func(S - Smax) * max(Peq - phi, 0)
        @hydroflux Fi ~ Peq - RH
        @hydroflux RD ~ step_func(S - Smax) * Peq
        @hydroflux Et ~ (Ep - Ei) * clamp((S - sw) / (sm - sw), 0, 1)
        @hydroflux L ~ Ksat * max(0.0, S)^c
    end

    dfluxes = begin
        @stateflux S ~ Fi - RD - Et - L
    end
end

bucket_4 = @hydrobucket :bucket_4 begin
    fluxes = begin
        @hydroflux Ls ~ min(L, Lmax)
        @hydroflux Lf ~ L - Ls
        @hydroflux Rf ~ kf * Sfast
    end

    dfluxes = begin
        @stateflux Sfast ~ Lf - Rf
    end
end

bucket_5 = @hydrobucket :bucket_5 begin
    fluxes = begin
        @hydroflux Rs ~ ks * Sslow
        @hydroflux Qt ~ RH + RD + Rf + Rs
    end

    dfluxes = begin
        @stateflux Sslow ~ Ls - Rs
    end
end
model = @hydromodel :echo begin
    bucket_1
    bucket_2
    bucket_3
    bucket_4
    bucket_5
end

end