module mcrm
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables S [description = "current interception storage", unit = "mm"]
@variables D [description = "current storage in soil moisture", unit = "mm"]
@variables Sg [description = "current groundwater storage", unit = "mm"]
@variables Sic [description = "current in-channel storage", unit = "mm"]
@variables Soc [description = "current out-of-channel storage", unit = "mm"]

@variables P [description = "precipitation", unit = "mm/d"]
@variables Ec [description = "evaporation", unit = "mm/d"]
@variables qf [description = "throughfall", unit = "mm/d"]
@variables Ep [description = "evapotranspiration", unit = "mm/d"]

@variables qn [description = "net infiltration", unit = "mm/d"]
@variables Er [description = "evaporation", unit = "mm/d"]
@variables qd [description = "direct runoff", unit = "mm/d"]
@variables qp [description = "percolation", unit = "mm/d"]
@variables qr [description = "rapid runoff", unit = "mm/d"]
@variables qb [description = "groundwater flow", unit = "mm/d"]
@variables Qt [description = "total flow", unit = "mm/d"]

@variables uib [description = "delayed flow", unit = "mm/d"]
@variables qic [description = "in-channel flow", unit = "mm/d"]
@variables uob [description = "out-of-bank flow", unit = "mm/d"]
@variables qoc [description = "out-of-channel flow", unit = "mm/d"]

# Model parameters
@parameters Smax [description = "Maximum interception storage", bounds = (0, 5), unit = "mm"]
@parameters cmax [description = "Maximum fraction of area contributing to rapid runoff", bounds = (0.01, 0.99), unit = "-"]
@parameters ct [description = "Fraction of cmax that is the minimum contributing area", bounds = (0.01, 0.99), unit = "-"]
@parameters c1 [description = "Shape parameter for rapid flow distribution", bounds = (0, 2), unit = "1/mm"]
@parameters ce [description = "Shape parameter for evaporation", bounds = (0, 1), unit = "1/mm"]
@parameters dsurp [description = "Threshold for direct runoff", bounds = (1, 2000), unit = "mm"]
@parameters kd [description = "Direct runoff time parameter", bounds = (0, 1), unit = "1/d"]
@parameters gamd [description = "Direct runoff flow non-linearity", bounds = (1, 5), unit = "-"]
@parameters qpmax [description = "Maximum percolation rate", bounds = (0, 20), unit = "mm/d"]
@parameters kg [description = "Groundwater time parameter", bounds = (0, 1), unit = "1/d"]
@parameters tau [description = "Routing delay  ", bounds = (1, 120), unit = "d"]
@parameters sbf [description = "Maximum routing store depth", bounds = (1, 300), unit = "mm"]
@parameters kcr [description = "Channel flow time parameter", bounds = (0, 1), unit = "1/d"]
@parameters gamcr [description = "Channel flow non-linearity", bounds = (1, 5), unit = "-"]
@parameters kor [description = "Out-of-bank flow time parameter", bounds = (0, 1), unit = "1/d"]
@parameters gamor [description = "Out-of-bank flow non-linearity", bounds = (1, 5), unit = "-"]


bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Ec ~ step_func(S) * Ep
        @hydroflux qf ~ step_func(S - Smax) * P
    end
    dfluxes = begin
        @stateflux S ~ P - Ec - qf
    end
end

bucket2 = @hydrobucket :bucket2 begin
    fluxes = begin
        @hydroflux qr ~ min(cmax, ct + ct * exp(c1 * D)) * qf
        @hydroflux qn ~ qf - qr
        @hydroflux Er ~ 1 / (1 + exp(-c1 * D)) * (Ep - Ec)
        @hydroflux qd ~ kd * max(0.0, D - dsurp)^gamd
        @hydroflux qp ~ clamp(D / dsurp, 0, 1) * qpmax
    end
    dfluxes = begin
        @stateflux D ~ qn - Er - qd - qp
    end
end

bucket3 = @hydrobucket :bucket3 begin
    fluxes = begin
        @hydroflux qb ~ kg * max(0.0, Sg)^1.5
    end
    dfluxes = begin
        @stateflux Sg ~ qp - qb
    end
end

bucket4 = @hydrobucket :bucket4 begin
    fluxes = begin
        @hydroflux uib ~ qb
        @hydroflux uob ~ step_func(Sic - sbf) * uib
        @hydroflux qic ~ min(3 / 4 * Sic, kcr * max(0.0, Sic)^gamcr)
    end
    dfluxes = begin
        @stateflux Sic ~ uib - uob - qic
    end
end


bucket5 = @hydrobucket :bucket5 begin
    fluxes = begin
        @hydroflux qoc ~ min(3 / 4 * Soc, kor * max(0.0, Soc)^gamor)
    end
    dfluxes = begin
        @stateflux Soc ~ uob - qoc
    end
end

flux1 = @hydroflux Qt ~ qoc + qic

model = @hydromodel :mcrm begin
    bucket1
    bucket2
    bucket3
    bucket4
    bucket5
    flux1
end

end