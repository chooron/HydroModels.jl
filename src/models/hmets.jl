module HMets

using ..HydroModels

#--------------------------------------------------------------------------------
# Parameters
#--------------------------------------------------------------------------------

# Snow parameters
@parameters ddfmin [description = "Minimum degree-day-factor in mm/°C/day"]
@parameters ddfplus [description = "Maximum degree-day-factor in mm/°C/day (ddfmin + ddfplus = ddfmax)"]
@parameters Tbm [description = "Base melting temperature in °C"]
@parameters Kcum [description = "Empirical parameter for the calculation of the degree-day-factor in mm⁻¹"]
@parameters fcmin [description = "Minimum fraction for the snowpack water retention capacity"]
@parameters fcplus [description = "Maximum fraction of the snowpack water retention capacity (fcmin + fcplus = fcmax)"]
@parameters Ccum [description = "Parameter for the calculation of water retention capacity in mm⁻¹"]
@parameters Tbf [description = "Base refreezing temperature in °C"]
@parameters Kf [description = "Degree-day factor for refreezing in mm/°C/day"]
@parameters Fe [description = "Empirical exponent for the freezing equation"]

# Real evapotranspiration parameter
@parameters ETeff [description = "Fraction of the potential evapotranspiration"]

# Subsurface parameters
@parameters cr [description = "Fraction of the water for surface and delayed runoff"]
@parameters cvp [description = "Fraction of the water for groundwater recharge"]
@parameters cv [description = "Fraction of the water for hypodermic flow"]
@parameters cp [description = "Fraction of the water for groundwater flow"]
@parameters LVmax [description = "Maximum level of the vadose zone in mm"]
@parameters LPmax [description = "Maximum level of the phreatic zone in mm"]

# Unit hydrograph parameters
@parameters α1 [description = "Shape parameter for the gamma distribution used on the surface unit hydrograph"]
@parameters β1 [description = "Rate parameter for the gamma distribution used on the surface unit hydrograph"]
@parameters α2 [description = "Shape parameter for the gamma distribution used on the delayed unit hydrograph"]
@parameters β2 [description = "Rate parameter for the gamma distribution used on the delayed unit hydrograph"]


#--------------------------------------------------------------------------------
# Variables
#--------------------------------------------------------------------------------
# Forcings
@variables Tmean [description = "Mean daily temperature in °C"]
@variables Tmin [description = "Minimum daily temperature in °C"]
@variables SNF [description = "Snowfall in mm/day"]
@variables PET [description = "Potential Evapotranspiration in mm/day"]

# State variables
@variables SNW [description = "Snow water equivalent in mm"]
@variables LV [description = "Vadose zone level in mm"]
@variables LP [description = "Phreatic zone level in mm"]

# Fluxes and other internal variables
@variables TD [description = "Meandiurnal temperature in °C"]
@variables POR [description = "Potential amount of overnight refreezing in mm"]
@variables DDF [description = "Actual degree day factor in mm/°C/day"]
@variables PSM [description = "Potential snowmelt in mm"]
@variables WRF [description = "Water retention fraction"]
@variables WAR [description = "Water available for runoff and infiltration in mm"]
@variables RET [description = "Real Evapotranspiration in mm"]
@variables INF [description = "Infiltration in mm"]
@variables GR [description = "Groundwater recharge in mm"]
@variables H1 [description = "Surface runoff component from rainfall in mm"]
@variables H2 [description = "Surface runoff component from snowmelt in mm"]
@variables H3 [description = "Hypodermic flow in mm"]
@variables H4 [description = "Groundwater flow in mm"]
@variables Q [description = "Total runoff in mm"]


#--------------------------------------------------------------------------------
# Model definition
#--------------------------------------------------------------------------------

# Parameter-derived values
ddfmax = ddfmin + ddfplus
fcmax = fcmin + fcplus

snow_bucket = @hydrobucket begin
    fluxes = begin
        @hydroflux Tmean ~ (Tmax + Tmin) / 2
        @hydroflux TD ~ (Tmean + Tmin) / 2
        @hydroflux POR ~ max(0, Kf * (Tbf - TD)^Fe)
        @hydroflux DDF ~ ddfmin * (1 + Kcum * SNW)
        @hydroflux PSM ~ max(0, DDF * (Tmean - Tbm))
        @hydroflux WRF ~ max(fcmin, fcmax * (1 - Ccum * SNW))
        @hydroflux WAR ~ max(0, SNW + P - WRF)
    end

    dfluxes = begin
        @stateflux SNW ~ POR + SNF - PSM
    end
end

soil_bucket = @hydrobucket begin
    fluxes = begin
        @hydroflux RET ~ ETeff * PET
        @hydroflux H1 ~ cr * (LV / LVmax) * WAR
        @hydroflux INF ~ WAR - H1 - RET
        @hydroflux GR ~ cvp * LV
        @hydroflux H2 ~ cr * INF * (LV / LVmax)
        @hydroflux H3 ~ cv * LV
        @hydroflux H4 ~ cp * LP
        @hydroflux Q ~ H2 + max(0, LV - LVmax) + max(0, LP - LPmax) + H3 + H4 + H1
    end

    dfluxes = begin
        @stateflux LV ~ INF - RET - H1 - H2 - H3 - GR
        @stateflux LP ~ GR - H4
    end
end

model = @hydromodel begin
    snow_bucket
    soil_bucket
end

# todo unit hydrograph
end