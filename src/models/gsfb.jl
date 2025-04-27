
module gsfb
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables P [description = "precipitation", unit = "mm/d"]
@variables Ep [description = "potential rate", unit = "mm/d"]

@variables S [description = "current storage in the upper zone", unit = "mm"]
@variables SS [description = "current storage in the subsurface store", unit = "mm"]
@variables DS [description = "current deep storage", unit = "mm"]

@variables Qdr [description = "deep groundwater", unit = "mm/d"]
@variables Ea [description = "evaporation", unit = "mm/d"]
@variables Qs [description = "surface runoff", unit = "mm/d"]
@variables F [description = "infiltration", unit = "mm/d"]
@variables Qb [description = "baseflow", unit = "mm/d"]
@variables Dp [description = "deep percolation", unit = "mm/d"]
@variables Qt [description = "Total flow", unit = "mm/d"]

# Model parameters
@parameters C [description = "Recharge coefficient", bounds = (0, 1), unit = "d-1"]
@parameters NDC [description = "Recharge coefficient", bounds = (0.05, 0.95), unit = "-"]
@parameters Smax [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters Emax [description = "Maximum evaporation rate", bounds = (0, 20), unit = "mm/d"]
@parameters Frate [description = "Recharge rate", bounds = (0, 200), unit = "mm/d"]
@parameters B [description = "Fraction subsurface flow to stream", bounds = (0, 1), unit = "-"]
@parameters DPF [description = "Runoff coefficient", bounds = (0, 1), unit = "d-1"]
@parameters SDRmax [description = "Threshold for subsurface flow generation", bounds = (1, 20), unit = "mm"]

# Soil water component
bucket = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Qdr ~ C * DS * (1.0 - min(1.0, S / (NDC * Smax)))
        @hydroflux Ea ~ min(Ep, Emax * min(S / (NDC * Smax), 1.00))
        @hydroflux Qs ~ step_func(S - Smax) * P
        @hydroflux F ~ Frate * step_func(S - NDC * Smax)

        @hydroflux Qb ~ B * DPF * max(0.0, SS - SDRmax)
        @hydroflux Dp ~ (1 - B) * DPF * SS

        @hydroflux Qt ~ Qs + Qb
    end

    dfluxes = begin
        @stateflux S ~ P + Qdr - Ea - Qs - F
        @stateflux SS ~ F - Qb - Dp
        @stateflux DS ~ Dp - Qdr
    end
end

model = @hydromodel :gsfb begin
    bucket
end

end