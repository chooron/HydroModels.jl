module hbv_edu
using ..HydroModels
using ..HydroModels: step_func

# Model variables by component
# State variables
@variables soilwater snowpack meltwater suz slz
# Input variables
@variables P Ep T
# Snow component variables
@variables rainfall snowfall melt refreeze infil
# Soil component variables
@variables wetness excess recharge evap
# Response component variables
@variables q0 q1 q2 Qt perc

# Model parameters
@parameters TT [description = "Terature threshold for snowmelt", bounds = (-1.5, 1.2), unit = "C"]
@parameters CFMAX [description = "Maximum snowmelt capacity", bounds = (1.0, 8.0), unit = "mm"]
@parameters CWH [description = "Water holding capacity", bounds = (0.0, 0.2), unit = "mm"]
@parameters CFR [description = "Refreeze capacity", bounds = (0.0, 0.1), unit = "mm"]
@parameters FC [description = "Field capacity", bounds = (50.0, 500.0), unit = "mm"]
@parameters LP [description = "Evaporation capacity", bounds = (0.3, 1.0), unit = "mm"]
@parameters BETA [description = "Evaporation exponent", bounds = (1.0, 6.0), unit = "-"]
@parameters PPERC [description = "Percolation capacity", bounds = (0.0, 3.0), unit = "mm"]
@parameters UZL [description = "Upland zone length", bounds = (0, 70.0), unit = "d"]
@parameters k0 [description = "Routing constant", bounds = (0.05, 0.5), unit = "1/d"]
@parameters k1 [description = "Routing constant", bounds = (0.01, 0.3), unit = "1/d"]
@parameters k2 [description = "Routing constant", bounds = (0.001, 0.15), unit = "1/d"]

# Precipitation separation into rainfall and snowfall
split = @hydroflux begin
    snowfall ~ step_func(TT - T) * P
    rainfall ~ step_func(T - TT) * P
end

# Snow routine
snow_bucket = @hydrobucket :hbv_snow begin
    fluxes = begin
        @hydroflux melt ~ min(snowpack, max(0.0, T - TT) * CFMAX)
        @hydroflux refreeze ~ min(max((TT - T), 0.0) * CFR * CFMAX, meltwater)
        @hydroflux infil ~ max(0.0, meltwater - snowpack * CWH)
    end
    dfluxes = begin
        @stateflux snowpack ~ snowfall + refreeze - melt
        @stateflux meltwater ~ melt - refreeze - infil
    end
end

# Soil moisture routine
soil_bucket = @hydrobucket :hbv_soil begin
    fluxes = begin
        @hydroflux recharge ~ (rainfall + infil) * clamp((max(soilwater / FC, 0.0))^BETA, 0, 1)
        @hydroflux excess ~ max(soilwater - FC, 0.0)
        @hydroflux evap ~ clamp(soilwater / (LP * FC), 0, 1) * Ep
    end
    dfluxes = begin
        @stateflux soilwater ~ rainfall + infil - (recharge + excess + evap)
    end
end

# Response routine
zone_bucket = @hydrobucket :hbv_zone begin
    fluxes = begin
        @hydroflux perc ~ suz * PPERC
        @hydroflux q0 ~ max(0.0, suz - UZL) * k0
        @hydroflux q1 ~ suz * k1
        @hydroflux q2 ~ slz * k2
        @hydroflux Qt ~ q0 + q1 + q2
    end
    dfluxes = begin
        @stateflux suz ~ recharge + excess - (perc + q0 + q1)
        @stateflux slz ~ perc - q2
    end
end

route_uh = @unithydro 

# Complete model
model = @hydromodel :hbv begin
    split
    snow_bucket
    soil_bucket
    zone_bucket
end

end