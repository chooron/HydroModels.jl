module susannah1
using ..HydroModels

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

# Model variables
@variables P [description = "precipitation input", unit = "mm/d"]
@variables Ep [description = "The potential rate", unit = "mm/d"]

@variables Suz [description = "current storage in the upper zone", unit = "mm"]
@variables Ebs [description = "bare soil evaporation", unit = "mm/d"]
@variables Eveg [description = "transpiration from vegetation", unit = "-"]
@variables Qse [description = "saturation excess flow ", unit = "mm/d"]
@variables Qss [description = " non-linear subsurface flow", unit = "-"]
@variables S
@variables Sgw [description = "groundwater storage ", unit = "mm"]
@variables Qr
@variables Qb [description = " baseflow flux ", unit = "mm/d"]
@variables Qt

# Model parameters
@parameters Sb [description = "Maximumsoil moisture storage", bounds = (1, 2000), unit = "mm"]
@parameters Sfc [description = "Field capacity", bounds = (0.05, 0.95), unit = "mm"]
@parameters M [description = "forest fraction", bounds = (0.05, 0.95), unit = "-"]
@parameters a [description = "Runoff time coefficient", bounds = (1, 50), unit = "d"]
@parameters b [description = "Runoff nonlinearity", bounds = (0.2, 1), unit = "-"]
@parameters r [description = " Fraction subsurface flow to groundwater", bounds = (0, 1), unit = "-"]

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Ebs ~ Suz / Sb * (1 - M) * Ep
        @hydroflux Eveg ~ min(1.0, Suz / Sfc) * M * Ep
        @hydroflux Qse ~ step_func(Suz - Sb) * P
        @hydroflux Qss ~ (max(0.0, Suz - Sfc) / a)^(1 / b)
        @hydroflux Ebs ~ Suz / Sb * (1 - M) * Ep
        @hydroflux Qr ~ r * Qss
        @hydroflux Qb ~ ((1 / a) * max(0.0, Sgw))^(1 / b)
    end
    dfluxes = begin
        @stateflux Suz ~ P - Ebs - Eveg - Qse - Qss
        @stateflux Sgw ~ Qr - Qb
    end
end

model = @hydromodel :susannah1 begin
    bucket1
    @hydroflux Qt ~ Qse + (Qss - Qr) + Qb
end

end