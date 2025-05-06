module penman
using ..HydroModels

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

# Model variables
@variables Srz [description = "current storage in the root zone", unit = "mm"]
@variables P [description = "precipitation", unit = "mm/d"]
@variables Ea [description = "evaporation", unit = "mm/d"]
@variables qex [description = "moisture excess", unit = "mm/d"]
@variables Ep [description = "the potential rate", unit = "mm/d"]
@variables Sdef [description = "the current moisture deficit ", unit = "mm"]
@variables Et [description = " evaporation", unit = "mm/d"]
@variables u2 [description = "Saturation ", unit = "mm/d"]
@variables q12 [description = "inflow", unit = "mm/d"]
@variables Cres [description = "the current storage in the routing reservoir", unit = "mm"]
@variables u1 [description = "the fraction phi of qex", unit = "mm/d"]
@variables Qt [description = "runoff", unit = "mm/d"]

# Model parameters
@parameters Smax [description = "Maximumsoil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters phi [description = "Fraction of saturation excess that is direct runof", bounds = (0.0, 1.0), unit = "-"]
@parameters gam [description = " Evaporation reduction factor", bounds = (0.0, 1.0), unit = "-"]
@parameters k1 [description = "Runoffcoefficient", bounds = (0.0, 1.0), unit = "1/d"]

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Ea ~ step_func(Srz) * Ep
        @hydroflux qex ~ step_func(Srz - Smax) * P
        @hydroflux q12 ~ (1 - phi) * qex

        @hydroflux Et ~ (1 - step_func(Srz)) * gam * Ep
        @hydroflux u2 ~ min(Sdef, step_func(Sdef) * q12)
    end
    dfluxes = begin
        @stateflux Srz ~ P - Ea - qex
        @stateflux Sdef ~ q12 - (Et + u2)
    end
end

bucket2 = @hydrobucket :bucket2 begin
    fluxes = begin
        @hydroflux u1 ~ phi * qex
        @hydroflux Qt ~ k1 * Cres
    end
    dfluxes = begin
        @stateflux Cres ~ u1 + u2 - Qt
    end
end

model = @hydromodel :penman begin
    bucket1
    bucket2
end

end