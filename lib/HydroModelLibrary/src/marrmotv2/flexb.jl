@reexport module Flexb
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3
@variables prcp pet temp
@variables ru eur ps rf rs rfl rsl qf qs

@parameters s1max [description="Maximum soil moisture storage [mm]"]
@parameters beta [description="Unsaturated zone shape parameter [-]"]
@parameters d [description="Fast/slow runoff distribution parameter [-]"]
@parameters percmax [description="Maximum percolation rate [mm/d]"]
@parameters lp [description="Wilting point as fraction of s1max [-]"]
@parameters kf [description="Fast runoff coefficient [d-1]"]
@parameters ks [description="Slow runoff coefficient [d-1]"]
@parameters nlagf [description="Flow delay before fast runoff [d]"]
@parameters nlags [description="Flow delay before slow runoff [d]"]

# Bucket 1: Unsaturated zone (S1)
funcs_1 = [
    SaturationFlux((S=S1, in=prcp), (Smax=s1max, p1=beta), Val(3), saturation=ru),
    EvaporationFlux((S=S1, in=pet), (p1=lp, Smax=s1max), Val(3), evaporation=eur),
    PercolationFlux((S=S1,), (p1=percmax, Smax=s1max), Val(2), percolation=ps)
]

dfuncs_1 = [
    StateFlux([ru] => [eur, ps], S1)
]

# Bucket 2: Fast response reservoir (S2)
funcs_2 = [
    HydroFlux([prcp, ru] => [rf], [d], exprs=[(1-d) * (prcp - ru)]),
    BaseflowFlux((S=S2,), (p1=kf,), Val(1), baseflow=qf)
]

dfuncs_2 = [
    StateFlux([rfl] => [qf], S2)
]

# Bucket 3: Slow response reservoir (S3)
funcs_3 = [
    HydroFlux([prcp, ru] => [rs], [d], exprs=[d * (prcp - ru)]),
    BaseflowFlux((S=S3,), (p1=ks,), Val(1), baseflow=qs)
]

dfuncs_3 = [
    StateFlux([rsl] => [qs], S3)
]

# Unit Hydrograph routing
uh_fast = UnitHydrograph([rf] => [rfl], [nlagf], uhfunc=UHFunction(:UH_1_HALF), solvetype=:SPARSE)
uh_slow = UnitHydrograph([ps, rs] => [rsl], [nlags], uhfunc=UHFunction(:UH_2_FULL), solvetype=:SPARSE)

# Total flow calculation
q_flux = HydroFlux([qf, qs] => [q], exprs=[qf + qs])

# Create buckets and model
flexb_bucket_1 = HydroBucket(name=:flexb_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
flexb_bucket_2 = HydroBucket(name=:flexb_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
flexb_bucket_3 = HydroBucket(name=:flexb_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
flexb_model = HydroModel(name=:flexb, components=[flexb_bucket_1, flexb_bucket_2, flexb_bucket_3, uh_fast, uh_slow, q_flux])

export flexb_bucket_1, flexb_bucket_2, flexb_bucket_3, flexb_model
end
