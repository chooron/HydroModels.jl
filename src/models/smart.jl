module smart
#! cite: https://thibhlln.github.io/smartpy/model_description.html

using ..HydroModels
using ..HydroModels: step_func
# Model variables
@variables P [description = "Precipitation as rainfall", unit = "mm/d"]
@variables Ep [description = "Potential evapotranspiration", unit = "mm/d"]
@variables Qt [description = "Total flow", unit = "mm/d"]
@variables Ea [description = "Actual evapotranspiration", unit = "mm/d"]

@variables Pc [description = "corrected precipitation", unit = "mm/d"]
@variables Pe [description = "effective precipitation", unit = "mm/d"]
@variables Sprop [description = "proportion of soil moisture", unit = "-"]

@variables Sof [description = "surface flow", unit = "mm"]
@variables Sdf [description = "drainage flow", unit = "mm"]
@variables Sif [description = "interflow", unit = "mm"]
@variables Ssgw [description = "shallow groundwater flow", unit = "mm"]
@variables Sdgw [description = "deep groundwater flow", unit = "mm"]
@variables Sr [description = "flow", unit = "mm"]

Qs = @variables Q0 Q1 Q2 Q3 Q4 Q5 Q6
Soils = @variables S1 S2 S3 S4 S5 S6
Sifs = @variables Sif1 Sif2 Sif3 Sif4 Sif5 Sif6
Ssgws = @variables Ssgw1 Ssgw2 Ssgw3 Ssgw4 Ssgw5 Ssgw6
Sdgws = @variables Sdgw1 Sdgw2 Sdgw3 Sdgw4 Sdgw5 Sdgw6
Bs = @variables B0 B1 B2 B3 B4 B5 B6
Es = @variables E0 E1 E2 E3 E4 E5 E6
Ds = @variables D0 D1 D2 D3 D4 D5 D6

@variables cond [description = "condition", unit = "-"]
@variables Rdf [description = "quick runoff as drainflow", unit = "mm/d"]
@variables Rif [description = "Interflow", unit = "mm/d"]
@variables Rsgw [description = "shallow groundwater flow", unit = "mm/d"]
@variables Rdgw [description = "deep groundwater flow", unit = "mm/d"]

@variables Qof [description = "surface flow", unit = "mm/d"]
@variables Qdf [description = "drainage flow", unit = "mm/d"]
@variables Qif [description = "interflow", unit = "mm/d"]
@variables Qsgw [description = "shallow groundwater flow", unit = "mm/d"]
@variables Qdgw [description = "deep groundwater flow", unit = "mm/d"]

# Model parameters
@parameters racf [description = "Rainfall aerial correction factor", bounds = (0.9, 1.1), unit = "-"]
@parameters edc [description = "Evaporation decay coefficient", bounds = (0.0, 1.0), unit = "-"]
@parameters qrr [description = "Quick runoff ratio", bounds = (0, 0.3), unit = "-"]
@parameters dfr [description = "Drain flow ratio", bounds = (0, 1.0), unit = "-"]
@parameters soc [description = "Soil outflow coefficient", bounds = (0.0, 1.0), unit = "-"]
@parameters efd [description = "Effective soil depth", bounds = (15.0, 1500.0), unit = "mm"]
@parameters SK [description = "Surface reservoir residence time", bounds = (0.0, 1.0), unit = "-"]
@parameters FK [description = "Interflow reservoir residence time", bounds = (0.0, 1.0), unit = "-"]
@parameters GK [description = "Groundwater reservoir residence time", bounds = (0.0, 1.0), unit = "-"]
@parameters RK [description = "Channel reservoir residence time", bounds = (0.0, 1.0), unit = "-"]

# Soil water component
bucket_1 = @hydrobucket :bucket_1 begin
    fluxes = begin
        @hydroflux Pc ~ racf * P
        @hydroflux cond ~ step_func(Pc - Ep)
        @hydroflux Pe ~ max(Pc - Ep, 0)
        @hydroflux Ea ~ min(Ep, Pc)
        @hydroflux Sprop ~ clamp((S1 + S2 + S3 + S4 + S5 + S6) / efd, 0.0, 1.0)
        @hydroflux Q0 ~ (1 - qrr * Sprop) * Pe
        @hydroflux for i in 2:7
            Qs[i] ~ clamp(Qs[i-1] - (efd / 6 - Soils[i-1]), 0, Soils[i-1])
        end
        @hydroflux for i in 1:6
            Sifs[i] ~ cond * Soils[i] * (soc * Sprop)^i
            Ssgws[i] ~ cond * Soils[i] * (soc * Sprop) / i
            Sdgws[i] ~ cond * Soils[i] * (soc * Sprop)^(7 - i)
        end
        @hydroflux Rdf ~ Q6 * dfr
        @hydroflux Rif ~ cond * (Sif1 + Sif2 + Sif3 + Sif4 + Sif5 + Sif6)
        @hydroflux Rsgw ~ cond * (Ssgw1 + Ssgw2 + Ssgw3 + Ssgw4 + Ssgw5 + Ssgw6)
        @hydroflux Rdgw ~ cond * (Sdgw1 + Sdgw2 + Sdgw3 + Sdgw4 + Sdgw5 + Sdgw6)
        @hydroflux D0 ~ (1 - cond) * (Ep - Ea)
        @hydroflux B0 ~ (1 - cond)
        @hydroflux for i in 2:7
            Es[i] ~ Bs[i-1] * min(Soils[i-1], Ds[i-1])
            Ds[i] ~ max(0.0, Ds[i-1] - Soils[i-1]) * edc
            Bs[i] ~ step_func(Ds[i-1] - Soils[i-1]) * Bs[i-1]
        end
    end

    dfluxes = begin
        @stateflux S1 ~ Q0 - Sif1 - Ssgw1 - Sdgw1 - E1
        @stateflux S2 ~ Q1 - Sif2 - Ssgw2 - Sdgw2 - E2
        @stateflux S3 ~ Q2 - Sif3 - Ssgw3 - Sdgw3 - E3
        @stateflux S4 ~ Q3 - Sif4 - Ssgw4 - Sdgw4 - E4
        @stateflux S5 ~ Q4 - Sif5 - Ssgw5 - Sdgw5 - E5
        @stateflux S6 ~ Q5 - Sif6 - Ssgw6 - Sdgw6 - E6
    end
end

bucket_2 = @hydrobucket :bucket_2 begin
    fluxes = begin
        @hydroflux Qof ~ Sof * SK
        @hydroflux Qdf ~ Sdf * SK
        @hydroflux Qif ~ Sif * FK
        @hydroflux Qsgw ~ Ssgw * GK
        @hydroflux Qdgw ~ Sdgw * GK
        @hydroflux Qt ~ Sr * RK
    end

    dfluxes = begin
        @stateflux Sof ~ Pe - Q0 - Qof
        @stateflux Sdf ~ Rdf - Qdf
        @stateflux Sif ~ Q6 - Rdf + Rif - Qif
        @stateflux Ssgw ~ Rsgw - Qsgw
        @stateflux Sdgw ~ Rdgw - Qdgw
        @stateflux Sr ~ Qof + Qdf + Qif + Qsgw + Qdgw - Qt
    end
end

model = @hydromodel :smart begin
    bucket_1
    bucket_2
end

end