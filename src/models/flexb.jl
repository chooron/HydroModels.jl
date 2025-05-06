module flexb
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables UR [description = "the current storage in the unsaturated zone", unit = "mm"]
@variables Ru [description = " the inflow into UR based on its current storage", unit = "mm"]
@variables Eur [description = " the evaporation", unit = "mm/d"]
@variables Rp [description = "the percolation from UR to the slow reservoir SR", unit = "mm/d"]
@variables Cr [description = "the current storage in the slow reservoir SR", unit = "mm"]
@variables P [description = "routed towards the unsaturated zone based on Cr", unit = "mm/d"]
@variables Ep [description = "potential evapotranspiration", unit = "mm/d"]
@variables Ps [description = "the percolation from UR to the slow reservoir SR", unit = "mm/d"]
@variables Rs [description = "the remainder being divided into preferential recharge", unit = "mm/d"]
@variables Rf [description = " fast runoff ", unit = "mm/d"]
@variables FR [description = " thecurrent storage ", unit = "mm"]
@variables SR [description = " the current storage", unit = "mm"]
@variables Rfl [description = "Rf routed through unit hydrograph", unit = "mm/d"]
@variables Rsl [description = "Rs routed through unit hydrograph", unit = "mm/d"]
@variables Qf [description = "Outflow from fast reservoir", unit = "mm/d"]
@variables Qs [description = "Outflow from slow reservoir", unit = "mm/d"]
@variables Qt [description = "Total outflow from fast and slow reservoirs", unit = "mm/d"]
@variables t

# Model parameters
@parameters URmax [description = " Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters beta [description = "Shape parameter", bounds = (0, 10), unit = "-"]
@parameters D [description = " Fraction effective precipitation to slow store", bounds = (0, 1), unit = "-"]
@parameters Percmax [description = "Maximumpercolation rate", bounds = (0, 20), unit = "mm/d"]
@parameters Lp [description = "Wilting point as fraction of URmax", bounds = (0.05, 0.95), unit = "-"]
@parameters Nlagf [description = "Unit Hydrograph time base", bounds = (1, 5), unit = "d"]
@parameters Nlags [description = "Unit Hydrograph time base", bounds = (1, 15), unit = "d"]
@parameters Kf [description = "Runoff coefficient", bounds = (0, 1), unit = "1/d"]
@parameters Ks [description = "Runoff coefficient", bounds = (0, 1), unit = "1/d"]

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Cr ~ 1 / (1 + exp((-UR / URmax + 1 / 2) / beta))
        @hydroflux Ru ~ (1 - Cr) * P
        @hydroflux Eur ~ Ep * min(1, UR / (URmax * Lp))
        @hydroflux Ps ~ Percmax * UR / URmax
    end
    dfluxes = begin
        @stateflux UR ~ Ru - Eur - Ru
    end
end

flux1 = @hydroflux begin
    Rs ~ (P - Ru) * D + Ps
    Rf ~ (P - Ru) * (1 - D)
end

uh1 = @unithydro begin
    uh_func = begin
        Nlagf => 1 / (0.5 * Nlagf^2) * (0.5 * min(t, Nlagf)^2 + 0.5 * (t - 1)^2)
    end
    uh_vars = [Rf]
    configs = (solvetype=:SPARSE, outputs=[Rfl])
end

uh2 = @unithydro begin
    uh_func = begin
        Nlags => 1 / (0.5 * Nlags^2) * (0.5 * min(t, Nlags)^2 + 0.5 * (t - 1)^2)
    end
    uh_vars = [Rs]
    configs = (solvetype=:SPARSE, outputs=[Rsl])
end

bucket2 = @hydrobucket :bucket2 begin
    fluxes = begin
        @hydroflux Qf ~ Kf * FR
        @hydroflux Qs ~ Ks * SR
        @hydroflux Qt ~ Qf + Qs
    end
    dfluxes = begin
        @stateflux FR ~ Rfl - Qf
        @stateflux SR ~ Rsl - Qs
    end
end

model = @hydromodel :flexb begin
    bucket1
    flux1
    uh1
    uh2
    bucket2
end

end