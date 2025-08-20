module gr4j
using ..HydroModels

# Model variables
@variables P [description = "Precipitation input", unit = "mm/d"]
@variables Ep [description = "Potential evapotranspiration input", unit = "mm/d"]
@variables ps [description = "Net rainfall that enters the production store", unit = "mm/d"]
@variables pn [description = "Net rainfall (precipitation minus evapotranspiration when positive)", unit = "mm/d"]
@variables es [description = "Actual evaporation from the production store", unit = "mm/d"]
@variables en [description = "Net evaporation (evapotranspiration minus precipitation when positive)", unit = "mm/d"]
@variables perc [description = "Percolation from production store to routing", unit = "mm/d"]
@variables Q9 [description = "Slow flow component", unit = "mm/d"]
@variables Q1 [description = "Fast flow component", unit = "mm/d"]
@variables Q9_routed [description = "Slow flow component routed through unit hydrograph 1", unit = "mm/d"]
@variables Q1_routed [description = "Fast flow component routed through unit hydrograph 2", unit = "mm/d"]
@variables Qroute [description = "Outflow from routing store", unit = "mm/d"]
@variables Qt [description = "Total runoff", unit = "mm/d"]
@variables exch [description = "Water exchange between groundwater and surface water", unit = "mm/d"]
@variables t 

@variables S [description = "Production store level", unit = "mm"]
@variables R [description = "Routing store level", unit = "mm"]
# Model parameters
@parameters x1 [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters x2 [description = "Subsurface water exchange", bounds = (-20, 20), unit = "mm/d"]
@parameters x3 [description = "Routing store depth", bounds = (1, 300), unit = "mm"]
@parameters x4 [description = "Unit Hydrograph time base = (5, 15)", bounds = (5, 15), unit = "d"]

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux pn ~ max(0.0, P - Ep)
        @hydroflux en ~ max(0.0, Ep - P)
        @hydroflux ps ~ pn * (1 - (S / x1)^2)
        @hydroflux es ~ min(S, en * (2S / x1 - (S / x1)^2))
        @hydroflux perc ~ min(S, ((x1)^-4) / 4 * ((4 / 9)^-4) * (S^5))
    end
    dfluxes = begin
        @stateflux S ~ ps - es - perc
    end
end

split_flux = @hydroflux begin
    Q9 ~ (perc + pn - ps) * 0.9
    Q1 ~ (perc + pn - ps) * 0.1
end

uh_1 = @unithydro begin
    uh_func = begin
        x4 => (t / x4)^2.5
    end
    uh_vars = [Q9]
    configs = (solvetype=:SPARSE, suffix=:_routed)
end

uh_2 = @unithydro begin
    uh_func = begin
        2x4 => (1 - 0.5 * (2 - t / x4)^2.5)
        x4 => (0.5 * (t / x4)^2.5)
    end
    uh_vars = [Q1]
    configs = (solvetype=:SPARSE, suffix=:_routed)
end

bucket2 = @hydrobucket :bucket2 begin
    fluxes = begin
        @hydroflux exch ~ x2 * max(0.0, R / x3)^3.5
        @hydroflux Qroute ~ min(R, x3^(-4) / 4 * (max(0.0, R))^5)
        @hydroflux Qt ~ Qroute + max(Q1_routed + exch, 0.0)
    end
    dfluxes = begin
        @stateflux R ~ Q9_routed + exch - Qroute
    end
end

model = @hydromodel :gr4j begin
    bucket1
    split_flux
    uh_1
    uh_2
    bucket2
end

export model
end