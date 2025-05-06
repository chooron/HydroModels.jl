module sacramento
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables UZTW [description = "upper zone tension water", unit = "mm"]
@variables Peff [description = "effective precipitation", unit = "mm/d"]
@variables Ru [description = "redistribution of free water", unit = "mm/d"]
@variables Euztw [description = "drained by evaporation", unit = "mm/d"]
@variables Twexu [description = "tension water excess", unit = "mm/d"]
@variables P [description = "precipitation", unit = "mm/d"]
@variables Qdir [description = "direct runoff from impervious area", unit = "mm/d"]
@variables UZFW [description = "upper zone free water", unit = "mm"]
@variables Ep [description = "potential evapotranspiration", unit = "mm/d"]
@variables Euzfw [description = "evaporation from upper zone free water", unit = "mm/d"]
@variables Qsur [description = "surface runoff", unit = "mm/d"]
@variables Qint [description = "interflow", unit = "mm/d"]
@variables Pc [description = "percolation to deeper groundwater", unit = "mm/d"]
@variables Pcdemand [description = "percolation demand", unit = "mm/d"]
@variables LZdeficiency [description = "lower zone deficiency", unit = "mm"]
@variables LZcapacity [description = "lower zone capacity", unit = "mm"]
@variables Rprop1 Rprop2 Rprop3 [description = "redistribution proportions", unit = "-"]
@variables LZTW [description = "lower zone tension water", unit = "mm"]
@variables LZFWP [description = "lower zone primary free water storage", unit = "mm"]
@variables LZFWS [description = "lower zone supplemental free water storage", unit = "mm"]
@variables Rlp [description = "redistribution from lower zone primary free water", unit = "-"]
@variables Rls [description = "redistribution from lower zone supplemental free water", unit = "-"]
@variables Elztw [description = "evaporation from lower zone tension water", unit = "mm/d"]
@variables Twexl [description = "lower zone tension water excess", unit = "mm/d"]
@variables Qbfp [description = "primary baseflow", unit = "mm/d"]
@variables Qbfs [description = "supplemental baseflow", unit = "mm/d"]
@variables Qt [description = "total flow", unit = "mm/d"]
# Model parameters
@parameters PCTIM [description = "Fraction impervious area", bounds = (0, 1), unit = "-"]
@parameters UZTWM [description = "Maximumupperzonetension water storage", bounds = (0, 1), unit = "mm"]
@parameters kuz [description = "Runoff coefficient", bounds = (0, 1), unit = "1/d"]
@parameters REXP [description = "Base percolation rate non-linearity factor", bounds = (0, 7), unit = "-"]
@parameters PFREE [description = " Fraction of percolation to free storage", bounds = (0, 1), unit = "-"]
@parameters klzp [description = "Primary baseflow runoff coefficient", bounds = (0, 1), unit = "1/d"]
@parameters klzs [description = "Supplemental baseflow runoff coefficient", bounds = (0, 1), unit = "1/d"]
@parameters Smax [description = "Maximum total storage depth", bounds = (1, 2000), unit = "mm"]
@parameters f1 [description = "fraction of smax that is Maximum upper zone tension water storage", bounds = (0.005, 0.995), unit = "mm"]
@parameters f2 [description = "fraction of smax-S1max that is Maximum upper zone free water storage", bounds = (0.005, 0.995), unit = "mm"]
@parameters f3 [description = "fraction of smax-S1max-S2max that is  Maximum lower zone tension water storage", bounds = (0.005, 0.995), unit = "mm"]
@parameters f4 [description = "fraction of smax-S1max-S2max-S3max that is  Maximum lower zone primary free water storage", bounds = (0.005, 0.995), unit = "mm"]

UZFWM = max(0.005 / 4, f2 * (Smax - UZTWM)) # Maximum upper zone free water storage [mm]
LZTWM = max(0.005 / 4, f3 * (Smax - UZTWM - UZFWM)) # Maximum lower zone tension water storage [mm]
LZFWPM = max(0.005 / 4, f4 * (Smax - UZTWM - UZFWM - LZTWM))# Maximum lower zone primary free water storage [mm]
LZFWSM = max(0.005 / 4, (1 - f4) * (Smax - UZTWM - UZFWM - LZTWM))# Maximum lower zone supplemental free water storage [mm]
PBASE = LZFWPM * klzp + LZFWSM * klzs # Base percolation rate [mm/d]  a base percolation rate
ZPERC = min(100000, (LZTWM + LZFWSM * (1 - klzs)) / (LZFWSM * klzs + LZFWPM * klzp) + (LZFWPM * (1 - klzp)) / (LZFWSM * klzs + LZFWPM * klzp)) # Base percolation rate multiplication factor [-]

split_flux = @hydroflux begin
    Peff ~ (1 - PCTIM) * P
    Qdir ~ PCTIM * P
end

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux Euztw ~ UZTW / UZTWM * Ep
        @hydroflux Euzfw ~ max(0.0, Ep - Euztw) * step_func(UZFW)
        @hydroflux Elztw ~ step_func(LZTW) * max(0.0, Ep - Euztw - Euzfw) / (UZTWM + UZFWM)

        @hydroflux Pcdemand ~ PBASE * (1 + ZPERC * (max(0.0, ((LZTWM - LZTW) + (LZFWPM - LZFWP) + (LZFWSM - LZFWS)) / (LZTWM + LZFWPM + LZFWSM))^(1 + REXP)))
        @hydroflux Pc ~ Pcdemand * UZFW / UZFWM

        @hydroflux Twexu ~ step_func(UZTW - UZTWM) * Peff
        @hydroflux Twexl ~ step_func(LZTW - LZTWM) * (1 - PFREE) * Pc

        @hydroflux Rprop1 ~ clamp((-LZTW * (LZFWPM + LZFWSM) + LZTWM * (LZFWP + LZFWS)) / ((LZFWPM + LZFWSM) * (LZTWM + LZFWPM + LZFWSM)), 0.0, 1.0)
        @hydroflux Rprop2 ~ clamp((LZFWPM - LZFWP) / (LZFWPM * ((LZFWPM - LZFWP) / LZFWPM + (LZFWSM - LZFWS) / LZFWSM)), 0.0, 1.0)
        @hydroflux Rprop3 ~ clamp((LZFWSM - LZFWS) / (LZFWSM * ((LZFWPM - LZFWP) / LZFWPM + (LZFWSM - LZFWS) / LZFWSM)), 0.0, 1.0)

        @hydroflux Ru ~ max(0.0, (UZTWM * UZFW - UZFWM * UZTW) / (UZTWM + UZFWM))
        @hydroflux Rlp ~ LZFWPM * Rprop1
        @hydroflux Rls ~ LZFWSM * Rprop1

        @hydroflux Qsur ~ step_func(UZFW - UZFWM) * Twexu
        @hydroflux Qint ~ kuz * UZFW
        @hydroflux Qbfp ~ klzp * LZFWP
        @hydroflux Qbfs ~ klzs * LZFWS
    end
    dfluxes = begin
        @stateflux UZTW ~ Peff + Ru - Euztw - Twexu
        @stateflux UZFW ~ Twexu - Euzfw - Qsur - Qint - Pc - Ru
        @stateflux LZTW ~ (1 - PFREE) * Pc + Rlp + Rls - Elztw - Twexl
        @stateflux LZFWP ~ Rprop2 * (PFREE * Pc) + Rprop2 * Twexl - Qbfp
        @stateflux LZFWS ~ Rprop3 * (PFREE * Pc) + Rprop3 * Twexl - Qbfs
    end
end

flux1 = @hydroflux Qt ~ Qdir + Qsur + Qint + Qbfp + Qbfs

model = @hydromodel :sacramento begin
    split_flux
    bucket1
    flux1
end

end