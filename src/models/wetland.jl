module wetland
using HydroModels
using ..HydroModels: step_func

# Model variables
@variables P [description = "Incoming precipitation", unit = "[mm/d]"]
@variables Ep [description = "Evapotranspiration rate", unit = "[mm/d]"]
@variables Pe [description = "Efficiency precipitation", unit = "[mm/d]"]
@variables Ew [description = "Evaporation from soil moisture", unit = "[mm/d]"]
@variables Sw [description = "The current soil water storage", unit = "[mm]"]
@variables Qw_sof [description = "Saturation excess surface runoff", unit = "[mm/d]"]
@variables Qw_gw [description = "Groundwater flow", unit = "[mm/d]"]
@variables Qt [description = "Total flow", unit = "[mm/d]"]

# Model parameters
@parameters Dw [description = " Interception evaporation ", bounds = (0, 5), unit = "[mm/d]"]
@parameters Betaw [description = "Non-linearity parameter for contributing area", bounds = (0, 10), unit = "[-]"]
@parameters Sw_max [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "[mm]"]
@parameters Kw [description = "Runoff coefficient", bounds = (0, 1), unit = "[d-1]"]

# Soil water component
soil_bucket = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux Pe ~ max(P - Dw, 0)
        @hydroflux Ew ~ step_func(Sw) * Ep
        @hydroflux Qw_sof ~ (1 - (1 - Sw / Sw_max)^Betaw) * Pe
        @hydroflux Qw_gw ~ Kw * Sw
        @hydroflux Qt ~ Qw_gw + Qw_sof
    end

    dfluxes = begin
        @stateflux Sw ~ Pe - Ew - Qw_sof - Qw_gw
    end
end

model = @hydromodel :wetland begin
    soil_bucket
end

end