module hymod
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables Sm [description = "moisture storage", unit = "mm"]
@variables P [description = " the precipitation input", unit = "mm/d"]
@variables Ea [description = "evaporation", unit = "mm/d"]
@variables Pe [description = " a the fraction of Pe that flows into the fast stores", unit = "mm/d"]
@variables Ep [description = "evaporation", unit = "mm/d"]
@variables S [description = "the current storage in store S", unit = "mm"]
@variables F1 [description = "the current storage in store F1 ", unit = "mm"]
@variables Pf [description = "fast flow", unit = "mm/d"]
@variables Ps [description = "slow flow", unit = "mm/d"]
@variables Qf1 [description = "fast flow in store Sf1", unit = "mm/d"]
@variables Sf1 [description = "fast flow store 1", unit = "mm"]
@variables Qf2 [description = "fast flow in store Sf2", unit = "mm/d"]
@variables Sf2 [description = "fast flow store 2", unit = "mm"]
@variables Qf3 [description = "fast flow in store Sf3", unit = "mm/d"]
@variables Sf3 [description = "fast flow store 3", unit = "mm"]
@variables Qs [description = "slow flow", unit = "mm/d"]
@variables Qt [description = "total flow", unit = "mm/d"]

# Model parameters
@parameters Smax [description = "Maximumsoil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters b [description = "Contributing area curve shape parameter", bounds = (0, 10), unit = "-"]
@parameters a [description = "Fraction of effective precipitation that is fast flow", bounds = (0, 1), unit = "-"]
@parameters kf [description = "Runoffcoefficient", bounds = (0, 1), unit = "1/d"]
@parameters ks [description = "Runoffcoefficient", bounds = (0, 1), unit = "1/d"]

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Ea ~ min(Sm / Smax * Ep, Sm)
        @hydroflux Pe ~ (1 - max(0.0, 1 - Sm / Smax)^b) * P
        @hydroflux Pf ~ a * Pe
        @hydroflux Ps ~ (1 - a) * Pe
    end
    dfluxes = begin
        @stateflux Sm ~ P - Ea - Pe
    end
end

bucket2 = @hydrobucket :bucket2 begin
    fluxes = begin
        @hydroflux Qf1 ~ kf * Sf1
        @hydroflux Qf2 ~ kf * Sf2
        @hydroflux Qf3 ~ kf * Sf3
    end
    dfluxes = begin
        @stateflux Sf1 ~ Pf - Qf1
        @stateflux Sf2 ~ Qf1 - Qf2
        @stateflux Sf3 ~ Qf2 - Qf3
    end
end

bucket3 = @hydrobucket :bucket3 begin
    fluxes = begin
        @hydroflux Qs ~ ks * S
        @hydroflux Qt ~ Qs + Qf3
    end
    dfluxes = begin
        @stateflux S ~ Ps - Qs
    end
end

model = @hydromodel :hymod begin
    bucket1
    bucket2
    bucket3
end

end