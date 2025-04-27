module lascam
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables P [description = "precipitation", unit = "mm/d"]
@variables Ep [description = "potential rate", unit = "mm/d"]

@variables F [description = "the current storage in the unsaturated infiltration store", unit = "mm"]
@variables B [description = "the relative storage in groundwater", unit = "mm"]
@variables A [description = "the current storage in the more permeable upper zone (above less permeable lower zone F)", unit = "mm/d"]

@variables fa [description = "actual infiltration", unit = "mm/d"]
@variables Aprop [description = "(A - amin) / (amax - amin)", unit = "-"]
@variables rf [description = "recharge", unit = "mm/d"]
@variables Ei [description = "leaf evaporation", unit = "mm/d"]
@variables Eb [description = "evaporation", unit = "mm/d"]
@variables φss [description = "the fraction saturated catchment area", unit = "-"]
@variables φc [description = "overland flow", unit = "-"]
@variables fss [description = "catchment-scale infiltration capacity", unit = "mm/d"]
@variables F [description = "the relative infiltration volume in the catchment", unit = "mm"]
@variables Pc [description = "the lesser of throughfall rate Pg minus saturation excess qse", unit = "mm/d"]
@variables Pg [description = "throughfall rate", unit = "mm/d"]
@variables qse [description = "saturation excess", unit = "mm/d"]
@variables fs [description = "catchment infiltration capacity", unit = "mm/d"]
@variables Ef [description = "Evaporation", unit = "mm/d"]
@variables Ep [description = "potential rate", unit = "mm/d"]
@variables qsse [description = "sub-surface saturation excess", unit = "mm/d"]
@variables qsie [description = "sub-surface infiltration excess", unit = "mm/d"]
@variables qb [description = "discharge from groundwater", unit = "mm/d"]
@variables Ea [description = "evaporation", unit = "mm/d"]
@variables qa [description = "subsurface stormflow", unit = "mm/d"]
@variables ra [description = "recharge", unit = "mm/d"]
@variables qie [description = "infiltration excess on the surface.", unit = "mm/d"]
@variables Qt [description = "total runoff", unit = "mm/d"]

# Model parameters
@parameters αf [description = "Catchment-scale infiltration parameter", bounds = (0, 200), unit = "mm/d"]
@parameters βf [description = "Catchment-scale infiltration non-linearity parameter", bounds = (0, 5), unit = "-"]
@parameters stot [description = "Total catchment storage", bounds = (1, 2000), unit = "mm"]
@parameters xa [description = "Fraction of Stot that is Amax", bounds = (0.01, 0.99), unit = "-"]
@parameters xf [description = "Fraction of Stot-Amax that is depth Fmax", bounds = (0.01, 0.99), unit = "-"]
@parameters na [description = "Fraction of Amax that is Amin", bounds = (0.01, 0.99), unit = "-"]
@parameters αc [description = "Variable contributing area scaling", bounds = (0, 5), unit = "-"]
@parameters βc [description = "Variable contributing area non-linearity", bounds = (0, 10), unit = "-"]
@parameters αss [description = "Subsurface saturation area scaling", bounds = (0, 5), unit = "-"]
@parameters βss [description = "Subsurface saturation area non-linearity", bounds = (0, 10), unit = "-"]
@parameters c [description = "Maximum infiltration rate", bounds = (0, 200), unit = "mm/d"]
@parameters αg [description = "Interception base parameter", bounds = (0, 5), unit = "mm/d"]
@parameters βg [description = "Interception fraction parameter", bounds = (0, 1), unit = "-"]
@parameters γf [description = "F-store evaporation scaling", bounds = (0, 1), unit = "-"]
@parameters δf [description = "F-store evaporation non-linearity", bounds = (0, 10), unit = "-"]
@parameters rd [description = "Recharge time parameter", bounds = (0, 1), unit = "d-1"]
@parameters αb [description = "Groundwater flow scaling", bounds = (0, 1), unit = "-"]
@parameters βb [description = "Groundwater flow base rate", bounds = (0.01, 200), unit = "mm/d"]
@parameters γa [description = "A-store evaporation scaling", bounds = (0, 1), unit = "-"]
@parameters δa [description = "A-store evaporation non-linearity", bounds = (0, 10), unit = "-"]
@parameters αa [description = "Subsurface storm flow rate", bounds = (0.01, 200), unit = "mm/d"]
@parameters βa [description = "Subsurface storm flow non-linearity", bounds = (1, 5), unit = "-"]
@parameters γb [description = "B-store evaporation scaling", bounds = (0, 1), unit = "-"]
@parameters δb [description = "B-store evaporation non-linearity", bounds = (0, 10), unit = "-"]

amax = xa * stot
fmax = xf * (stot - amax)
bmax = (1 - xf) * (stot - amax)
amin = na * amax

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Aprop ~ max(0.0, A - amin) / (amax - amin)
        @hydroflux φc ~ min(1.0, αc * Aprop^βc)
        @hydroflux φss ~ min(1.0, αss * Aprop^βss)
        @hydroflux fss ~ clamp(αf * (1 - (B / bmax)) * max(0.0, F / fmax)^(-βf), 0, 1e6)
        @hydroflux Pg ~ max(βg * P - αg, 0)
        @hydroflux Ei ~ P - Pg
        @hydroflux qse ~ φc * Pg
        @hydroflux Pc ~ min(Pg - qse, c)
        @hydroflux fa ~ min(Pc * max(1, (1 - φss) / (1 - φc)), fss)
        @hydroflux Ef ~ min(F, γf * Ep * (max(0.0, F / fmax)^δf))
        @hydroflux rf ~ rd * F

        @hydroflux qsse ~ (φss - φc) / (1 - φc) * Pc
        @hydroflux qsie ~ max(Pc * (1 - φss) / (1 - φc) - fss, 0)
        @hydroflux qb ~ min(B, βb * (exp(αb * B / bmax) - 1))
        @hydroflux Ea ~ min(A, φc * Ep + γa * Ep * (max(0.0, A / amax)^δa))
        @hydroflux qa ~ min(A, αa * Aprop^βa)
        @hydroflux ra ~ min(A, φss * fss)

        @hydroflux Eb ~ min(B, γb * Ep * (max(0.0, B / bmax)^δb))
        @hydroflux qie ~ Pg - qse - Pc
        @hydroflux Qt ~ qse + qie + qa
    end
    dfluxes = begin
        @stateflux F ~ fa - Ef - rf
        @stateflux A ~ qsse + qsie + qb - Ea - qa - ra
        @stateflux B ~ rf + ra - Eb - qb
    end
end

model = @hydromodel :lascam begin
    bucket1
end

end