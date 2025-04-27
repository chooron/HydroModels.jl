module ihacres
using ..HydroModels
using ..HydroModels: step_func

using Integrals
# Model variables
@variables CMD [description = "current moisture deficit", unit = "mm"]
@variables P [description = "incoming precipitation that reduces the deficit", unit = "mm/d"]
@variables Ea [description = "evaporation that increases the deficit", unit = "mm/d"]
@variables U [description = "effective precipitation that occurs when the deficit is below a threshold d", unit = "mm/d"]
@variables Ep [description = "Evapotranspiration rate", unit = "mm"]
@variables Q [description = "Xs+Xq", unit = "mm"]
@variables t
@variables Uq
@variables Us
@variables Xs
@variables Xq

# Model parameters
@parameters lp [description = "Wilting point", bounds = (1, 2000), unit = "mm"]
@parameters d [description = "Deficit threshold for flow from rain", bounds = (1, 2000), unit = "mm"]
@parameters p [description = "Deficit non-linearity", bounds = (0, 10), unit = "-"]
@parameters alpha [description = "Fraction flow to quick routing", bounds = (0, 1), unit = "-"]
@parameters Tq [description = " Unit Hydrograph time base", bounds = (1, 700), unit = "d"]
@parameters Ts [description = "Unit Hydrograph time base", bounds = (1, 700), unit = "d"]
@parameters Td [description = "Unit Hydrograph delay", bounds = (0, 119), unit = "d"]

@parameters tau_q
@parameters tau_s
@parameters tau_d

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Ea ~ Ep * min(1, exp(2 * (1 - CMD / lp)))
        @hydroflux U ~ P * (1 - min(1, min((CMD / d)^p)))
        @hydroflux Uq ~ alpha * U
        @hydroflux Us ~ (1 - alpha) * U
    end
    dfluxes = begin
        @stateflux CMD ~ -P + Ea + U
    end
end

exp_uh_func(t, ps) = begin
    prob = IntegralProblem((x, p) -> exp(-x), (t, t + 1), ())
    sol = solve(prob, QuadGKJL())
    return sol(t, ps)
end

uh_1 = UnitHydrograph(
    [Uq], [tau_q],
    uh_func=exp_uh_func,
    max_lag_func=(ps) -> ceil(ps.tau_q),
    configs=(solvetype=:SPARSE, outputs=[Xq])
)

uh_2 = UnitHydrograph(
    [Us], [tau_s],
    uh_func=exp_uh_func,
    max_lag_func=(ps) -> ceil(ps.tau_s),
    configs=(solvetype=:SPARSE, outputs=[Xs])
)


model = @hydromodel :ihacres begin
    bucket1
    uh_1
    uh_2
    @variables Qt ~ Xq + Xs
end

end