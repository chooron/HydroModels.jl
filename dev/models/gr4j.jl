using HydroModels
using ModelingToolkit: @parameters, @variables

# Parameters in the GR4J model
@parameters x1 x2 x3 x4 lag area_coef

# Variables in different components
# Production store
@variables prcp ep soilwater pn en ps es perc pr slowflow fastflow
# Unit hydrograph
@variables slowflow_routed fastflow_routed
# Routing store
@variables routingstore exch routedflow flow q q_routed s_river

# Production store bucket
prod_bucket = @hydrobucket :gr4j_prod begin
    fluxes = [
        @hydroflux begin
            pn ~ prcp - min(prcp, ep)
            en ~ ep - min(prcp, ep)
        end,
        @hydroflux ps ~ max(0.0, pn * (1 - (soilwater / x1)^2)),
        @hydroflux es ~ en * (2 * soilwater / x1 - (soilwater / x1)^2),
        @hydroflux perc ~ ((x1)^(-4)) / 4 * ((4 / 9)^(4)) * (soilwater^5),
        @hydroflux pr ~ pn - ps + perc,
        @hydroflux begin
            slowflow ~ 0.9 * pr
            fastflow ~ 0.1 * pr
        end
    ]
    dfluxes = [
        @stateflux soilwater ~ ps - (es + perc)
    ]
end

# Routing store bucket
routing_bucket = @hydrobucket :gr4j_rst begin
    fluxes = [
        @hydroflux exch ~ x2 * abs(routingstore / x3)^3.5,
        @hydroflux routedflow ~ x3^(-4) / 4 * (routingstore + slowflow + exch)^5,
        @hydroflux flow ~ routedflow + max(fastflow_routed + exch, 0.0)
    ]
    dfluxes = [
        @stateflux routingstore ~ slowflow + exch - routedflow
    ]
end

# Unit hydrograph components
uh_flux_1 = HydroModels.UnitHydroFlux(slowflow, slowflow_routed, x4, 
    uhfunc=HydroModels.UHFunction(:UH_1_HALF), solvetype=:SPARSE)
uh_flux_2 = HydroModels.UnitHydroFlux(fastflow, fastflow_routed, x4, 
    uhfunc=HydroModels.UHFunction(:UH_2_FULL), solvetype=:SPARSE)

# Create the complete model
gr4j_model = @hydromodel :gr4j begin
    prod_bucket
    uh_flux_1
    uh_flux_2
    routing_bucket
end

export gr4j_model, prod_bucket, routing_bucket, uh_flux_1, uh_flux_2