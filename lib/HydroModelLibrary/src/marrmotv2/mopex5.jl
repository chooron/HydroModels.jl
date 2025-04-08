@reexport module Mopex5
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4 S5
@variables prcp pet temp t
@variables epc ps pr qn et1 i q1f qw et2 q2f q2u qf qs

# Basic parameters
@parameters tcrit [description = "Snowfall & snowmelt temperature [oC]"]
@parameters ddf [description = "Degree-day factor for snowmelt [mm/oC/d]"]
@parameters s2max [description = "Maximum soil moisture storage [mm]"]
@parameters tw [description = "Groundwater leakage time [d-1]"]
@parameters i_alpha [description = "Intercepted fraction of Pr [-]"]
@parameters i_s [description = "Maximum Leaf Area Index timing [d]"]
@parameters tmin [description = "Growing Season Index minimum temperature"]
@parameters trange [description = "Growing Season Index temperature range"]
@parameters tu [description = "Slow flow routing response time [d-1]"]
@parameters se [description = "Root zone storage capacity as fraction of s3max [-]"]
@parameters s3max [description = "Maximum groundwater storage [mm]"]
@parameters tc [description = "Mean residence time [d-1]"]

# Auxiliary parameters
tmax = 365.25  # Duration of seasonal cycle [d]

# Bucket 1: Snow storage
funcs_1 = [
    PhenologyFlux((T=temp), (tmin=tmin, tmax=tmin+trange, in=pet), Val(1), phenology=epc),
    SnowfallFlux((in=prcp, T=temp), (tt=tcrit,), Val(1), snowfall=ps),
    SnowmeltFlux((S=S1, T=temp), (ddf=ddf, tt=tcrit), Val(1), snowmelt=qn)
]

dfuncs_1 = [
    StateFlux([ps] => [qn], S1)
]

# Bucket 2: Upper soil moisture storage
funcs_2 = [
    RainfallFlux((in=prcp, T=temp), (tt=tcrit,), Val(1), rainfall=pr),
    EvaporationFlux((S=S2, in=epc), (Smax=s2max,), Val(7), evaporation=et1),
    InterceptionFlux((in=pr, t=t), (p1=i_alpha, p2=i_s, p3=tmax), Val(4), interception=i),
    SaturationFlux((S=S2, in=pr+qn), (Smax=s2max,), Val(1), saturation=q1f),
    RechargeFlux((S=S2,), (p1=tw,), Val(3), recharge=qw)
]

dfuncs_2 = [
    StateFlux([pr, qn] => [et1, i, q1f, qw], S2)
]

# Bucket 3: Lower soil moisture storage
funcs_3 = [
    EvaporationFlux((S=S3, in=epc), (Smax=se*s3max,), Val(7), evaporation=et2),
    SaturationFlux((S=S3, in=qw), (Smax=s3max,), Val(1), saturation=q2f),
    BaseflowFlux((S=S3,), (p1=tu,), Val(1), baseflow=q2u)
]

dfuncs_3 = [
    StateFlux([qw] => [et2, q2f, q2u], S3)
]

# Bucket 4: Fast flow routing
funcs_4 = [
    BaseflowFlux((S=S4,), (p1=tc,), Val(1), baseflow=qf)
]

dfuncs_4 = [
    StateFlux([q1f, q2f] => [qf], S4)
]

# Bucket 5: Slow flow routing
funcs_5 = [
    BaseflowFlux((S=S5,), (p1=tc,), Val(1), baseflow=qs)
]

dfuncs_5 = [
    StateFlux([q2u] => [qs], S5)
]

# Total flow calculation
q_flux = HydroFlux([qf, qs] => [q], exprs=[qf + qs])

# Create buckets and model
mopex5_bucket_1 = HydroBucket(name=:mopex5_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
mopex5_bucket_2 = HydroBucket(name=:mopex5_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
mopex5_bucket_3 = HydroBucket(name=:mopex5_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
mopex5_bucket_4 = HydroBucket(name=:mopex5_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)
mopex5_bucket_5 = HydroBucket(name=:mopex5_bucket_5, funcs=funcs_5, dfuncs=dfuncs_5)

mopex5_model = HydroModel(name=:mopex5, components=[mopex5_bucket_1, mopex5_bucket_2, mopex5_bucket_3, 
                                                   mopex5_bucket_4, mopex5_bucket_5, q_flux])

export mopex5_bucket_1, mopex5_bucket_2, mopex5_bucket_3, mopex5_bucket_4, mopex5_bucket_5, mopex5_model
end
