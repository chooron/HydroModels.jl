using ModelingToolkit

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

# Model variables by component
# State variables
@variables soilwater snowpack meltwater suz slz
# Input variables
@variables prcp pet temp
# Snow component variables
@variables rainfall snowfall melt refreeze infil
# Soil component variables
@variables wetness excess recharge evap
# Response component variables
@variables q0 q1 q2 q perc

# Model parameters
@parameters TT CFMAX CFR CWH     # Snow routine parameters
@parameters LP FC BETA           # Soil routine parameters
@parameters PPERC UZL           # Response routine parameters
@parameters k0 k1 k2 kp         # Routing parameters

# Precipitation separation into rainfall and snowfall
split = @hydroflux begin
    snowfall ~ step_func(TT - temp) * prcp
    rainfall ~ step_func(temp - TT) * prcp
end

# Snow routine
snow_bucket = @hydrobucket :hbv_snow begin
    fluxes = [
        @hydroflux melt ~ min(snowpack, max(0.0, temp - TT) * CFMAX),
        @hydroflux refreeze ~ min(max((TT - temp), 0.0) * CFR * CFMAX, meltwater),
        @hydroflux infil ~ max(0.0, meltwater - snowpack * CWH)
    ]
    dfluxes = [
        @stateflux snowpack ~ snowfall + refreeze - melt,
        @stateflux meltwater ~ melt - refreeze - infil
    ]
end

# Soil moisture routine
soil_bucket = @hydrobucket :hbv_soil begin
    fluxes = [
        @hydroflux recharge ~ (rainfall + infil) * clamp((max(soilwater / FC, 0.0))^BETA, 0, 1),
        @hydroflux excess ~ max(soilwater - FC, 0.0),
        @hydroflux evap ~ clamp(soilwater / (LP * FC), 0, 1) * pet
    ]
    dfluxes = [
        @stateflux soilwater ~ rainfall + infil - (recharge + excess + evap)
    ]
end

# Response routine
zone_bucket = @hydrobucket :hbv_zone begin
    fluxes = [
        @hydroflux perc ~ suz * PPERC,
        @hydroflux q0 ~ max(0.0, suz - UZL) * k0,
        @hydroflux q1 ~ suz * k1,
        @hydroflux q2 ~ slz * k2,
        @hydroflux q ~ q0 + q1 + q2
    ]
    dfluxes = [
        @stateflux suz ~ recharge + excess - (perc + q0 + q1),
        @stateflux slz ~ perc - q2
    ]
end

# Complete model
hbv_model = @hydromodel :hbv begin
    split
    snow_bucket
    soil_bucket
    zone_bucket
end

export hbv_model, snow_bucket, soil_bucket, zone_bucket