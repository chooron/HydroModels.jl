module tcm
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables Pin [description = "infiltrated precipitation", unit = "mm/d"]
@variables Ea [description = "evaporation", unit = "mm/d"]
@variables qex1 [description = "storage excess flow", unit = "mm/d"]
@variables Pn [description = "net precipitation", unit = "mm/d"]
@variables P [description = "difference between precipitation", unit = "mm/d"]
@variables Ep [description="potential evapotranspiration", unit="mm/d"]
@variables Sdef [description = "the current storage in the soil moisture deficit store", unit = "mm"]
@variables Et [description = "evaporation", unit = "mm/d"]
@variables qex2 [description = "percolation", unit = "mm/d"]
@variables Suz [description = " current storage in the unsaturated zone", unit = "mm"]
@variables Pby [description = "preferential recharge ", unit = "mm/d"]
@variables quz [description = "drained by groundwater flow", unit = "mm/d"]
@variables Qt [description = "has a quadratic relation with storage", unit = "1/(mm*d)"]

@variables Ssz [description = "the current storage in the saturated zone", unit = "mm"]
@variables Suz [description = "the current storage in the unsaturated zone", unit = "mm"]
@variables Srz [description = "the current storage in the root zone", unit = "mm"]


# Model parameters
@parameters ca [description = "Abstraction rate", unit = "mm/d", bounds = (0, 100)]
@parameters phi [description = "Fraction of net precipitation that is preferntial flow", bounds = (0, 1), unit = "-"]
@parameters rc [description = " Maximum root zone storage", bounds = (1, 2000), unit = "mm"]
@parameters gam [description = " Transpiration reduction factor", bounds = (0, 1), unit = "-"]
@parameters k1 [description = " Runoff coefficient", bounds = (0, 1), unit = "1/d"]
@parameters fa [description = "Abstraction rate", bounds = (0, 1), unit = "mm/d"]
@parameters k2 [description = " Runoff coefficient", bounds = (0, 1), unit = "1/d"]

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Pn ~ max(P - Ep, 0)
        @hydroflux Pin ~ (1 - phi) * Pn
        @hydroflux Ea ~ step_func(Srz) * Ep
        @hydroflux qex1 ~ step_func(Srz - rc) * Pin
    end
    dfluxes = begin
        @stateflux Srz ~ Pin - Ea - qex1
    end
end

bucket2 = @hydrobucket :bucket2 begin
    fluxes = begin
        @hydroflux Et ~ (1 - step_func(Srz)) * gam * Ep
        @hydroflux qex2 ~ (1 - step_func(Sdef)) * qex1
    end
    dfluxes = begin
        @stateflux Sdef ~ Et + qex1 - qex2
    end
end

bucket3 = @hydrobucket :bucket3 begin
    fluxes = begin
        @hydroflux Pby ~ phi * Pn
        @hydroflux quz ~ k1 * Suz
    end
    dfluxes = begin
        @stateflux Suz ~ Pby + qex2 - quz
    end
end

bucket4 = @hydrobucket :bucket4 begin
    fluxes = begin
        @hydroflux Qt ~ min(step_func(Ssz) * k2 * (Ssz^2), Ssz)
    end
    dfluxes = begin
        @stateflux Ssz ~ quz - min(ca, Ssz) - Qt
    end
end

model = @hydromodel :tcm begin
    bucket1
    bucket2
    bucket3
    bucket4
end

end