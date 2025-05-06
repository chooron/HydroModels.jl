module hillslope
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables Sw [description = "current soil water storage", unit = "mm"]
@variables Pe [description = "effective precipitation", unit = "mm/d"]
@variables Ea [description = "evaporation", unit = "mm/d"]
@variables Qse [description = " Storage excess surface runof", unit = "mm/d"]
@variables P [description = "Incoming precipitation", unit = "mm/d"]
@variables Dh [description = "interception", unit = "mm/d"]
@variables cap [description = "capillary rise", unit = "mm/d"]
@variables Ep [description = "the potential rate", unit = "mm/d"]
@variables Sh [description = "soil moisture storage", unit = "mm"]
@variables Shgw [description = "current groundwater storage", unit = "mm"]
@variables Qses [description = "surface fraction of storage excess flow", unit = "mm/d"]
@variables Qhsrf [description = "surface fraction of storage excess flow routed through unit hydrograph", unit = "mm/d"]
@variables Qseg [description = "groundwater fraction of storage excess flow", unit = "mm/d"]
@variables Qhgw [description = "Groundwater flow", unit = "mm/d"]
@variables Qt [description = "Total runoff", unit = "mm/d"]
@variables t
# Model parameters
@parameters Dw [description = "Interception evaporation", bounds = (0, 5), unit = "mm/d"]
@parameters Shmax [description = "Maximum soil moisture storage", bounds = (0, 10), unit = "mm"]
@parameters Betah [description = "Non-linearity parameter for contributing area", bounds = (1, 2000), unit = "-"]
@parameters a [description = "Fraction saturation excess to groundwater", bounds = (0, 1), unit = "-"]
@parameters Th [description = "Unit Hydrograph time base", bounds = (1, 120), unit = "d"]
@parameters C [description = " Capillary rise", bounds = (0, 4), unit = "mm/d"]
@parameters Kh [description = " Runoff coefficient", bounds = (0, 1), unit = "1/d"]

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Pe ~ max(P - Dw, 0.0)
        @hydroflux Ea ~ min(step_func(Sw) * Ep, Sw)
        @hydroflux Qse ~ (1 - max(0.0, 1 - Sw / Shmax)^Betah) * Pe
        @hydroflux cap ~ min(C, Shgw)

        @hydroflux Qses ~ a * Qse
        @hydroflux Qseg ~ (1 - a) * Qse
        @hydroflux Qhgw ~ Kh * Shgw
    end
    dfluxes = begin
        @stateflux Sw ~ Pe + cap - Ea - Qse
        @stateflux Shgw ~ Qseg - cap - Qhgw
    end
end

uh1 = @unithydro begin
    uh_func = begin
        Th => 1 / (0.5 * Th^2) * (0.5 * min(t, Th)^2 + 0.5 * (t - 1)^2)
    end
    uh_vars = [Qses]
    configs = (solvetype=:SPARSE, outputs=[Qhsrf])
end

flux1 = @hydroflux Qt ~ Qhsrf + Qhgw

model = @hydromodel :hillslope begin
    bucket1
    uh1
    flux1
end

end