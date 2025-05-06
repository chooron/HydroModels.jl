module nam
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables S [description = "Snow storage", unit = "mm"]
@variables Ps [description = "Precipitation falling as snow", unit = "mm/d"]
@variables M [description = "Snowmelt rate", unit = "mm/d"]
@variables T [description = "Temperature", unit = "°C"]
@variables P [description = "Precipitation", unit = "mm/d"]
@variables U [description = "Upper zone storage", unit = "mm"]
@variables Pr [description = "Precipitation as rain", unit = "mm/d"]
@variables Eu [description = "Evaporation from upper zone", unit = "mm/d"]
@variables If [description = "Interflow from upper zone", unit = "mm/d"]
@variables Pn [description = "Net precipitation excess", unit = "mm/d"]
@variables Ep [description = "Potential evapotranspiration", unit = "mm/d"]
@variables L [description = "Lower zone storage", unit = "mm"]
@variables Dl [description = "Infiltration to lower zone", unit = "mm/d"]
@variables Et [description = "Evapotranspiration from lower zone", unit = "mm/d"]
@variables Of [description = "Overland flow", unit = "mm/d"]
@variables O [description = "Overland flow routing storage", unit = "mm"]
@variables Qo [description = "Routed overland flow", unit = "mm/d"]
@variables I [description = "Interflow routing storage", unit = "mm"]
@variables Qi [description = "Routed interflow", unit = "mm/d"]
@variables G [description = "Groundwater storage", unit = "mm"]
@variables Gw [description = "Groundwater recharge", unit = "mm/d"]
@variables Qb [description = "Baseflow", unit = "mm/d"]
@variables Qt [description = "Total streamflow", unit = "mm/d"]

# Model parameters
@parameters Cs [description = "Degree-day factor", bounds = (0, 20), unit = "mm/oC/d"]
@parameters Cif [description = "Runoff coefficient", bounds = (0, 1), unit = "1/d"]
@parameters stot [description = "Maximum lower zone storage", bounds = (1, 2000), unit = "mm"]
@parameters fl [description = "Fraction of total soil depth that makes up lstar", bounds = (0.01, 0.99), unit = "mm"]
@parameters CL1 [description = "Fractional threshold for interflow generation", bounds = (0, 0.99), unit = "-"]
@parameters Cof [description = "Runoff coefficient", bounds = (0, 1), unit = "1/d"]
@parameters CL2 [description = "Fractional threshold for overland flow generation", bounds = (0, 0.99), unit = "-"]
@parameters K0 [description = "Runoff coefficient", bounds = (0, 1), unit = "1/d"]
@parameters K1 [description = "Runoff coefficient", bounds = (0, 1), unit = "1/d"]
@parameters Kb [description = "Runoff coefficient", bounds = (0, 1), unit = "1/d"]

Lmax = fl * stot     # Maximum lower zone storage [mm]
Umax = (1 - fl) * stot  # Upper zone maximum storage [mm]

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Ps ~ step_func(-T) * P
        @hydroflux M ~ min(Cs * max(0.0, T), S)
    end
    dfluxes = begin
        @stateflux S ~ Ps - M
    end
end

bucket2 = @hydrobucket :bucket2 begin
    fluxes = begin
        @hydroflux Pr ~ step_func(T) * P
        @hydroflux Eu ~ step_func(U) * Ep
        @hydroflux If ~ Cif * max(0.0, (L / Lmax - CL1) / (1 - CL1)) * U
        @hydroflux Pn ~ step_func(U - Umax) * (Pr + M)

        @hydroflux Of ~ step_func(L / Lmax - CL2) * (Cof * (L / Lmax - CL2) / (1 - CL2) * Pn)
        @hydroflux Dl ~ (Pn - Of) * (1 - L / Lmax)
        @hydroflux Et ~ step_func(U) * ((L / Lmax) * Ep)
    end
    dfluxes = begin
        @stateflux U ~ Pr + M - Eu - If - Pn
        @stateflux L ~ Dl - Et
    end
end

bucket3 = @hydrobucket :bucket3 begin
    fluxes = begin
        @hydroflux Qo ~ K0 * O
        @hydroflux Qi ~ K1 * I
        @hydroflux Gw ~ (Pn - Of) * (L / Lmax)
        @hydroflux Qb ~ Kb * O
        @hydroflux Qt ~ Qo + Qi + Qb
    end
    dfluxes = begin
        @stateflux O ~ Of - Qo
        @stateflux I ~ If - Qi
        @stateflux G ~ Gw - Qb
    end
end

model = @hydromodel :nam begin
    bucket1
    bucket2
    bucket3
end

end