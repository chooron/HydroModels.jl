module hbv
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables P [description = "Incoming precipitation", unit = "mm/d"]
@variables T [description = "Temperature", unit = "°C"]
@variables Ep [description = "Evaporation potential", unit = "mm/d"]

@variables SP [description = "current snow storage", unit = "mm"]
@variables WC [description = "current liquid storage", unit = "mm"]
@variables SM [description = "current soil moisture storage", unit = "mm"]
@variables UZ [description = "current upper zone storage", unit = "mm"]
@variables LZ [description = "current low zone storage", unit = "mm"]

@variables sf [description = "precipitation that occurs as snowfall", unit = "mm/d"]
@variables rf [description = "precipitation that occurs as rainfall", unit = "mm/d"]
@variables refr [description = "refreezing of liquid snow", unit = "mm/d"]
@variables melt [description = "snowmelt", unit = "mm/d"]
@variables infi [description = "infiltration to soil moisture", unit = "mm/d"]
@variables excess [description = "excess runoff", unit = "mm/d"]
@variables evap [description = "evaporation from soil moisture", unit = "mm/d"]
@variables cf [description = "capillary", unit = "mm/d"]
@variables r [description = "flow to the upper zone", unit = "mm/d"]
@variables perc [description = "percolation to the lower zone", unit = "mm/d"]
@variables q0 [description = "Outflow from the upper zone", unit = "mm/d"]
@variables q1 [description = "Outflow from the lower zone", unit = "mm/d"]
@variables q [description = "total outflow", unit = "mm/d"]
@variables q_routed [description = "Routed outflow", unit = "mm/d"]
@variables t
# Model parameters
@parameters TT [description = "Threshold temperature for snowfall", bounds = (-3, 5), unit = "°C"]
@parameters TTI [description = "Interval length of rain-snow spectrum", bounds = (0, 17), unit = "°C"]
@parameters TTM [description = "Threshold temperature for snowmelt", bounds = (-3, 3), unit = "°C"]
@parameters CFR [description = "Coefficient of refreezing of melted snow", bounds = (0, 1), unit = "-"]
@parameters CFMAX [description = "Degree-day factor of snowmelt and refreezing", bounds = (0, 20), unit = "mm/°C/d"]
@parameters WHC [description = "Maximum water holding content of snow pack", bounds = (0, 1), unit = "-"]
@parameters CFLUX [description = "Maximum rate of capillary rise", bounds = (0, 4), unit = "mm/d"]
@parameters FC [description = "Maximum soil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters LP [description = "Wilting point as fraction of FC", bounds = (0.05, 0.95), unit = "-"]
@parameters BETA [description = "Non-linearity coefficient of upper zone recharge", bounds = (0, 10), unit = "-"]
@parameters K0 [description = "Runoff coefficient from upper zone", bounds = (0, 1), unit = "d⁻¹"]
@parameters ALPHA [description = "Non-linearity coefficient of runoff from upper zone", bounds = (0, 4), unit = "-"]
@parameters PPERC [description = "Maximum rate of percolation to lower zone", bounds = (0, 20), unit = "mm/d"]
@parameters K1 [description = "Runoff coefficient from lower zone", bounds = (0, 1), unit = "d⁻¹"]
@parameters MAXBAS [description = "Flow routing delay", bounds = (1, 120), unit = "d"]

bucket_01 = @hydrobucket :snow begin
    fluxes = begin
        @hydroflux sf ~ P * clamp((TT + TTI / 2) - T, 0, 1)
        @hydroflux rf ~ P * clamp(T - (TT + TTI / 2), 0, 1)
        @hydroflux refr ~ clamp(CFR * CFMAX * (TTM - T), 0, WC)
        @hydroflux melt ~ clamp(CFMAX * (T - TTM), 0, SP)
        @hydroflux infi ~ step_func(WC - WHC * SP) * (rf + melt)
        @hydroflux excess ~ max(0.0, WC - WHC * SP)
    end
    dfluxes = begin
        @stateflux SP ~ sf + refr - melt
        @stateflux WC ~ rf + melt - refr - infi - excess
    end
end

bucket_02 = @hydrobucket :snow begin
    fluxes = begin
        @hydroflux cf ~ min(UZ, CFLUX * (1 - SM / FC))
        @hydroflux evap ~ min(Ep * min(1.0, SM / (FC * LP)), SM)
        @hydroflux r ~ (infi + excess) * max(0.0, SM / FC)^BETA

        @hydroflux q0 ~ K0 * (max(0.0, UZ)^(1 + ALPHA))
        @hydroflux perc ~ min(PPERC, UZ)
        @hydroflux q1 ~ K1 * LZ

        @hydroflux q ~ q0 + q1
    end
    dfluxes = begin
        @stateflux SM ~ infi + excess + cf - evap - r
        @stateflux UZ ~ r - cf - q0 - perc
        @stateflux LZ ~ perc - q1
    end
end

uh_1 = @unithydro begin
    uh_func = begin
        ceil(MAXBAS) => 1 / (0.5 * MAXBAS^2) * (0.5 * MAXBAS^2 - 0.5 * (t - 1)^2)
        MAXBAS => 1 / (0.5 * MAXBAS^2) * (0.5 * t^2 - 0.5 * (t - 1)^2)
    end
    uh_vars = [q]
    configs = (solvetype=:SPARSE, suffix=:_routed)
end

model = @hydromodel :hbv begin
    bucket_01
    bucket_02
    uh_1
end

end