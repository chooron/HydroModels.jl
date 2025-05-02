module ihm19
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables SI [description = "current interception storage", unit = "-"]
@variables P [description = "precipitation", unit = "-"]
@variables Ep [description = "potential evaporation", unit = "-"]
@variables EI [description = "interception evaporation", unit = "-"]
@variables PEX [description = "throughfal", unit = "-"]
@variables PEXMP [description = "pacropore excess precipitation", unit = "-"]
@variables SMP [description = "macropore storage", unit = "-"]
@variables PEXS1 [description = "soil storage", unit = "-"]

@variables FMP [description = "macropore infiltration", unit = "-"]
@variables QMP [description = "macropore storage", unit = "-"]
@variables QEXMP [description = "macropore excess flow", unit = "-"]
@variables PQEXS1 [description = "first soil layer", unit = "-"]

@variables SS1 [description = "upper soil storage", unit = "-"]
@variables FS1 [description = "soil infiltration", unit = "-"]
@variables QS1 [description = "soil runoff", unit = "-"]
@variables ETAS1 [description = "evapotranspiration", unit = "-"]
@variables FS1 [description = "soil infiltration rate", unit = "-"]
@variables FF [description = "forest fraction", unit = "-"]
@variables Q0 [description = "surface runoff", unit = "-"]
@variables QMPS1 [description = "first soil storage", unit = "-"]
@variables SS2 [description = "ower soil store", unit = "-"]

@variables PC [description = "percolation", unit = "-"]
@variables QH [description = "produced between the soil layers", unit = "-"]
@variables QS2 [description = "second soil layer", unit = "-"]
@variables Qt [description = "total runoff", unit = "-"]
@variables Q0R [description = "surface runoff", unit = "-"]

# Model parameters
@parameters SIMAX [description = "maximum interception storage", bounds = (2, 5), unit = "mm"]
@parameters A [description = "splitting coeffcient for excess precipitation", bounds = (0.9, 1), unit = "-"]
@parameters FF [description = "forest fraction", bounds = (0.4, 0.95), unit = "-"]
@parameters SMPMAX [description = "maximum storage macropores", bounds = (0.05, 5), unit = "mm"]
@parameters CQMP [description = "runoff time parameter (fast/slow runnoff) first soil layer", bounds = (0, 1), unit = "1/d"]
@parameters XQMP [description = "runoff scale parameter first soil layer", bounds = (1, 5), unit = "-"]
@parameters SS1MAX [description = "maximum soil moisture storage first soil layer", bounds = (400, 600), unit = "mm"]
@parameters FCCS1 [description = "field capacity as fraction of maximum storage first soil layer", bounds = (0.3, 0.7), unit = "-"]
@parameters CFS1 [description = "maximum infiltration rate first soil layer", bounds = (0, 1000), unit = "mm/d"]
@parameters XFS1 [description = "infiltration loss exponent first soil layer", bounds = (0, 15), unit = "-"]
@parameters CQS1 [description = "runoff time parameter for (fast/slow runnoff) first soil layer", bounds = (0, 1), unit = "1/d"]
@parameters XQS1 [description = "runoff scale parameter first soil layer", bounds = (1, 5), unit = "-"]
@parameters SS2MAX [description = "maximum soil moisture storage second soil layer", bounds = (300, 500), unit = "mm"]
@parameters CQS2 [description = "runoff time parameter for (fast/slow runnoff) second soil layer", bounds = (0, 1), unit = "1/d"]
@parameters XQS2 [description = "runoff scale parameter second soil layer", bounds = (1, 5), unit = "-"]
@parameters D0 [description = "Flow delay before surface runoff", bounds = (2, 5), unit = "d"]

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux EI ~ step_func(SI) * Ep
        @hydroflux PEX ~ step_func(SI - SIMAX) * P
        @hydroflux PEXMP ~ (1 - A) * PEX
        @hydroflux PEXS1 ~ A * PEX
    end
    dfluxes = begin
        @stateflux SI ~ P - EI - PEX
    end
end

bucket2 = @hydrobucket :bucket2 begin
    fluxes = begin
        @hydroflux FMP ~ step_func(SMPMAX - SMP) * PEXMP
        @hydroflux QMP ~ CQMP * max(0.0, SMP)^XQMP
        @hydroflux QEXMP ~ PEXMP - FMP
        @hydroflux PQEXS1 ~ PEXS1 + QEXMP
    end
    dfluxes = begin
        @stateflux SMP ~ FMP - QMP
    end
end

bucket3 = @hydrobucket :bucket3 begin
    fluxes = begin
        @hydroflux FS1 ~ min(PQEXS1, step_func(SS1MAX - SS1) * CFS1 * exp(-XFS1 * SS1 / SS1MAX))
        @hydroflux ETAS1 ~ (1 - FF) * SS1 / SS1MAX * Ep + FF * min(1.0, SS1 / (FCCS1 * SS1MAX)) * Ep
        @hydroflux QS1 ~ CQS1 * max(0.0, SS1 - FCCS1 * SS1MAX)^XQS1
        @hydroflux Q0 ~ PQEXS1 - FS1
        @hydroflux QMPS1 ~ QMP + QS1
    end
    dfluxes = begin
        @stateflux SS1 ~ FS1 - ETAS1 - QS1
    end
end

@variables Q0_routed t

uh_1 = @unithydro begin
    uh_func = begin
        D0 => (t / D0)^2.5
    end
    uh_vars = [Q0]
    configs = (solvetype=:SPARSE, suffix=:_routed)
end


bucket4 = @hydrobucket :bucket4 begin
    fluxes = begin
        @hydroflux PC ~ step_func(SS2MAX - SS2) * QMPS1
        @hydroflux QS2 ~ CQS2 * max(0.0, SS2)^XQS2
        @hydroflux QH ~ QMPS1 - PC
        @hydroflux Qt ~ Q0_routed + QH + QS2
    end
    dfluxes = begin
        @stateflux SS2 ~ PC - QS2
    end
end

model = @hydromodel :ihm19 begin
    bucket1
    bucket2
    bucket3
    uh_1
    bucket4
end

end