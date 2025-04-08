@reexport module Flexis
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4 S5
@variables prcp pet temp
@variables ps pi m peff ei ru eur rp rf rs rfl rsl qf qs

# Basic parameters
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
@parameters tt [description = "Threshold temperature for snowfall/snowmelt [oC]"]
@parameters ddf [description = "Degree-day factor for snowmelt [mm/d/oC]"]

# Bucket 1: Snow storage
funcs_1 = [
    SnowfallFlux((in=prcp, T=temp), (tt=tt,), Val(1), snowfall=ps),
    SnowmeltFlux((S=S1, T=temp), (ddf=ddf, tt=tt), Val(1), snowmelt=m)
]

dfuncs_1 = [
    StateFlux([ps] => [m], S1)
]

# Bucket 2: Interception storage
funcs_2 = [
    RainfallFlux((in=prcp, T=temp), (tt=tt,), Val(1), rainfall=pi),
    InterceptionFlux((S=S2, in=m+pi), (Smax=imax,), Val(1), interception=peff),
    EvaporationFlux((S=S2, in=pet), NamedTuple(), Val(1), evaporation=ei)
]

dfuncs_2 = [
    StateFlux([m, pi] => [peff, ei], S2)
]

# Bucket 3: Soil moisture storage
funcs_3 = [
    SaturationFlux((S=S3, in=peff), (Smax=smax, p1=beta), Val(3), saturation=ru),
    EvaporationFlux((S=S3, in=pet), (Smax=smax, p1=lp), Val(3), evaporation=eur),
    PercolationFlux((S=S3, in=smax), (p1=percmax,), Val(2), percolation=rp),
    SplitFlux((in=peff-ru), (p1=1-d,), Val(1), split=rf),
    SplitFlux((in=peff-ru), (p1=d,), Val(1), split=rs)
]

dfuncs_3 = [
    StateFlux([ru] => [eur, rp], S3)
]

# Bucket 4: Fast response reservoir
funcs_4 = [
    RouteFlux((in=rf,), (nlag=nlagf,), Val(1), route=rfl),
    BaseflowFlux((S=S4,), (p1=kf,), Val(1), baseflow=qf)
]

dfuncs_4 = [
    StateFlux([rfl] => [qf], S4)
]

# Bucket 5: Slow response reservoir
funcs_5 = [
    RouteFlux((in=rs+rp,), (nlag=nlags,), Val(1), route=rsl),
    BaseflowFlux((S=S5,), (p1=ks,), Val(1), baseflow=qs)
]

dfuncs_5 = [
    StateFlux([rsl] => [qs], S5)
]

# Total flow calculation
q_flux = HydroFlux([qf, qs] => [q], exprs=[qf + qs])

# Create buckets and model
flexis_bucket_1 = HydroBucket(name=:flexis_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
flexis_bucket_2 = HydroBucket(name=:flexis_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
flexis_bucket_3 = HydroBucket(name=:flexis_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
flexis_bucket_4 = HydroBucket(name=:flexis_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)
flexis_bucket_5 = HydroBucket(name=:flexis_bucket_5, funcs=funcs_5, dfuncs=dfuncs_5)

flexis_model = HydroModel(name=:flexis, components=[flexis_bucket_1, flexis_bucket_2, flexis_bucket_3, 
                                                   flexis_bucket_4, flexis_bucket_5, q_flux])

export flexis_bucket_1, flexis_bucket_2, flexis_bucket_3, flexis_bucket_4, flexis_bucket_5, flexis_model
end
