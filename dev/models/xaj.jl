using HydroModels

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

# Model variables
## Precipitation variables
@variables prcp [description = "Observed precipitation intensity", unit="mm/d"]
@variables pn [description = "Net precipitation intensity"]
@variables iu [description = "Upper to lower soil layer percolation intensity"]
@variables il [description = "Lower to deep soil layer percolation intensity"]
@variables fw [description = "Effective infiltration coefficient"]

## Evapotranspiration variables
@variables pet [description = "Potential evapotranspiration"]
@variables en [description = "Net evapotranspiration intensity"]
@variables eu [description = "Upper soil layer evapotranspiration intensity"]
@variables el [description = "Lower soil layer evapotranspiration intensity"]
@variables ed [description = "Deep soil layer evapotranspiration intensity"]
@variables et [description = "Actual evapotranspiration intensity"]

## Runoff variables
@variables r [description = "Total runoff generation intensity"]
@variables rs [description = "Surface runoff intensity (permeable area)"]
@variables ri [description = "Interflow intensity"]
@variables rg [description = "Groundwater runoff intensity"]
@variables rt [description = "Surface runoff intensity (total basin)"]

## Streamflow variables
@variables qi [description = "Interflow linear reservoir outflow intensity"]
@variables qg [description = "Groundwater linear reservoir outflow intensity"]
@variables qt [description = "Total channel network inflow intensity"]
@variables q1 [description = "Nash unit hydrograph 1st reservoir outflow intensity"]
@variables q2 [description = "Nash unit hydrograph 2nd reservoir outflow intensity"]
@variables q3 [description = "Nash unit hydrograph 3rd reservoir outflow intensity"]
@variables q [description = "Channel network outflow intensity"]

# Model parameters
## Evapotranspiration parameters
@parameters Ke [description = "Evapotranspiration potential coefficient", bounds = (0.6, 1.5)]
@parameters c [description = "Deep evapotranspiration coefficient", bounds = (0.01, 0.2)]

## Soil layer parameters
@parameters Wum [description = "Upper layer soil average potential water capacity", bounds = (5.0, 30.0)]
@parameters Wlm [description = "Lower layer soil average potential water capacity", bounds = (60.0, 90.0)]
@parameters Wdm [description = "Deep soil average potential water capacity", bounds = (15.0, 60.0)]

## Surface parameters
@parameters Aimp [description = "Impervious area ratio", bounds = (0.01, 0.2)]
@parameters b [description = "Potential water capacity curve index", bounds = (0.1, 0.4)]
@parameters Smax [description = "Average free water capacity", bounds = (10.0, 50.0)]
@parameters ex [description = "Free water capacity curve index", bounds = (1.0, 1.5)]

## Discharge parameters
@parameters Ki [description = "Sediment discharge coefficient", bounds = (0.1, 0.55)]
@parameters Kg [description = "Groundwater discharge coefficient", bounds = (0.1, 0.55)]
@parameters ci [description = "Subsurface runoff dissipation coefficient", bounds = (0.5, 0.9)]
@parameters cg [description = "Groundwater runoff dissipation coefficient", bounds = (0.98, 0.998)]
@parameters Kf [description = "Subsurface runoff discharge coefficient", bounds = (0.1, 10.0)]

## State variables
@variables wu [description = "Upper layer soil average potential water capacity"]
@variables wl [description = "Lower layer soil average potential water capacity"]
@variables wd [description = "Deep soil average potential water capacity"]
@variables s0 [description = "Free water capacity"]
@variables oi [description = "Interflow"]
@variables og [description = "Groundwater"]
@variables F1 [description = "Linear reservoir 1"]
@variables F2 [description = "Linear reservoir 2"]
@variables F3 [description = "Linear reservoir 3"]


eu ~ step_func(wu) * en


e_us_ei ~ alpha_ei * p

# Soil water component
soil_bucket = @hydrobucket :soil begin
    fluxes = [
        @hydroflux begin
            pn ~ max(0.0, prcp - Ke * pet)
            en ~ max(0.0, Ke * pet - prcp)
        end,
        @hydroflux eu ~ step_func(wu) * en,
        @hydroflux el ~ step_func(wl) * max(c, wl / Wlm) * (en - eu),
        @hydroflux ed ~ step_func(wd) * max(c * (en - eu) - el, 0.0),
        @hydroflux et ~ eu + el + ed,
        @hydroflux fw ~ (1 - Aimp) * (1 - (1 - (wu + wl + wd) / (Wum + Wlm + Wdm))^(b / (1 + b))),
        @hydroflux r ~ pn * (fw + Aimp),
        @hydroflux iu ~ step_func(wu - Wum) * (pn - r - eu),
        @hydroflux il ~ step_func(wl - Wlm) * (iu - el)
    ]
    dfluxes = [
        @stateflux wu ~ pn - (r + eu + iu),
        @stateflux wl ~ iu - (el + il),
        @stateflux wd ~ il - ed
    ]
end

# Free water component
free_water_bucket = @hydrobucket :zone begin
    fluxes = [
        @hydroflux rs ~ pn * fw * (1 - (1 - s0 / Smax)^(ex / (1 + ex))),
        @hydroflux begin
            ri ~ fw * s0 * (-Ki * log(1 - Ki - Kg) / (Ki + Kg))
            rg ~ fw * s0 * (-Kg * log(1 - Ki - Kg) / (Ki + Kg))
        end,
        @hydroflux rt ~ pn * Aimp + rs
    ]
    dfluxes = [
        @stateflux s0 ~ pn - (rs + ri + rg) / fw
    ]
end

# Land routing component
land_routing_bucket = @hydrobucket :landroute begin
    fluxes = [
        @hydroflux qi ~ -oi * log(ci),
        @hydroflux qg ~ -og * log(cg),
        @hydroflux qt ~ rt + qi + qg
    ]
    dfluxes = [
        @stateflux oi ~ ri - qi,
        @stateflux og ~ rg - qg
    ]
end

# River routing component
river_routing_bucket = @hydrobucket :riverroute begin
    fluxes = [
        @hydroflux q1 ~ F1 * Kf,
        @hydroflux q2 ~ F2 * Kf,
        @hydroflux q3 ~ F3 * Kf,
        @hydroflux q ~ q3
    ]
    dfluxes = [
        @stateflux F1 ~ qt - q1,
        @stateflux F2 ~ q1 - q2,
        @stateflux F3 ~ q2 - q3
    ]
end

# Complete model
xaj_model = @hydromodel :xaj begin
    soil_bucket
    free_water_bucket
    land_routing_bucket
    river_routing_bucket
end

export xaj_model, soil_bucket, free_water_bucket, land_routing_bucket, river_routing_bucket