@reexport module Flexi
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4
@variables prcp pet temp
@variables peff ei ru eur ps rf rs rfl rsl qf qs

@parameters smax [description = "Maximum soil moisture storage [mm]"]
@parameters beta [description = "Unsaturated zone shape parameter [-]"]
@parameters d [description = "Fast/slow runoff distribution parameter [-]"]
@parameters percmax [description = "Maximum percolation rate [mm/d]"]
@parameters lp [description = "Wilting point as fraction of s1max [-]"]
@parameters nlagf [description = "Flow delay before fast runoff [d]"]
@parameters nlags [description = "Flow delay before slow runoff [d]"]
@parameters kf [description = "Fast runoff coefficient [d-1]"]
@parameters ks [description = "Slow runoff coefficient [d-1]"]
@parameters imax [description = "Maximum interception storage [mm]"]

# Bucket 1: Interception store (S1)
funcs_1 = [
    InterceptionFlux((S=S1, in=prcp), (Smax=imax,), Val(1), interception=peff),
    EvaporationFlux((S=S1, in=pet), NamedTuple(), Val(1), evaporation=ei)
]

dfuncs_1 = [
    StateFlux([prcp] => [peff, ei], S1)
]

# Bucket 2: Unsaturated zone (S2)
funcs_2 = [
    SaturationFlux((S=S2, in=peff), (Smax=smax, p1=beta), Val(3), saturation=ru),
    EvaporationFlux((S=S2, in=pet), (p1=lp, Smax=smax), Val(3), evaporation=eur),
    PercolationFlux((S=S2,), (p1=percmax, Smax=smax), Val(2), percolation=ps)
]

dfuncs_2 = [
    StateFlux([ru] => [eur, ps], S2)
]

# Bucket 3: Fast response reservoir (S3)
funcs_3 = [
    HydroFlux([peff, ru] => [rf], [d], exprs=[(1-d) * (peff - ru)]),
    BaseflowFlux((S=S3,), (p1=kf,), Val(1), baseflow=qf)
]

dfuncs_3 = [
    StateFlux([rfl] => [qf], S3)
]

# Bucket 4: Slow response reservoir (S4)
funcs_4 = [
    HydroFlux([peff, ru] => [rs], [d], exprs=[d * (peff - ru)]),
    BaseflowFlux((S=S4,), (p1=ks,), Val(1), baseflow=qs)
]

dfuncs_4 = [
    StateFlux([rsl] => [qs], S4)
]

# Unit Hydrograph routing
uh_fast = UnitHydrograph([rf] => [rfl], [nlagf], uhfunc=UHFunction(:UH_3_HALF), solvetype=:SPARSE)
uh_slow = UnitHydrograph([ps, rs] => [rsl], [nlags], uhfunc=UHFunction(:UH_3_HALF), solvetype=:SPARSE)

# Total flow calculation
q_flux = HydroFlux([qf, qs] => [q], exprs=[qf + qs])

# Create buckets and model
flexi_bucket_1 = HydroBucket(name=:flexi_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
flexi_bucket_2 = HydroBucket(name=:flexi_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
flexi_bucket_3 = HydroBucket(name=:flexi_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
flexi_bucket_4 = HydroBucket(name=:flexi_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)
flexi_model = HydroModel(name=:flexi, components=[flexi_bucket_1, flexi_bucket_2, flexi_bucket_3, flexi_bucket_4, uh_fast, uh_slow, q_flux])

export flexi_bucket_1, flexi_bucket_2, flexi_bucket_3, flexi_bucket_4, flexi_model
end
