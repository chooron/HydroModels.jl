@reexport module Mopex3
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4 S5
@variables prcp pet temp
@variables ps pr qn et1 q1f qw et2 q2f q2u qf qs

@parameters tcrit [description = "Snowfall & snowmelt temperature [oC]"]
@parameters ddf [description = "Degree-day factor for snowmelt [mm/oC/d]"]
@parameters s2max [description = "Maximum soil moisture storage [mm]"]
@parameters tw [description = "Groundwater leakage time [d-1]"]
@parameters tu [description = "Slow flow routing response time [d-1]"]
@parameters se [description = "Root zone storage capacity as fraction of s3max [-]"]
@parameters s3max [description = "Maximum groundwater storage [mm]"]
@parameters tc [description = "Mean residence time [d-1]"]

# Bucket 1: Snow storage
funcs_1 = [
    SnowfallFlux((in=prcp, T=temp), (p1=tcrit,), Val(1), snowfall=ps),
    RainfallFlux((in=prcp, T=temp), (p1=tcrit,), Val(1), rainfall=pr),
    MeltFlux((S=S1, T=temp), (p1=ddf, p2=tcrit), Val(1), melt=qn)
]

dfuncs_1 = [
    StateFlux([ps] => [qn], S1)
]

# Bucket 2: Soil moisture storage
funcs_2 = [
    EvaporationFlux((S=S2, in=pet), (Smax=s2max,), Val(7), evaporation=et1),
    SaturationFlux((S=S2, in=:($(pr)+$(qn))), (Smax=s2max,), Val(1), saturation=q1f),
    RechargeFlux((S=S2,), (p1=tw,), Val(3), recharge=qw)
]

dfuncs_2 = [
    StateFlux([pr, qn] => [et1, q1f, qw], S2)
]

# Bucket 3: Groundwater storage
funcs_3 = [
    EvaporationFlux((S=S3, in=pet), (Smax=:($(se)*$(s3max)),), Val(7), evaporation=et2),
    SaturationFlux((S=S3, in=qw), (Smax=s3max,), Val(1), saturation=q2f),
    BaseflowFlux((S=S3,), (p1=tu,), Val(1), baseflow=q2u)
]

dfuncs_3 = [
    StateFlux([qw] => [et2, q2f, q2u], S3)
]

# Linear reservoirs combined
funcs_4 = [
    BaseflowFlux((S=S4,), (p1=tc,), Val(1), baseflow=qf),
    BaseflowFlux((S=S5,), (p1=tc,), Val(1), baseflow=qs)
]

dfuncs_4 = [
    StateFlux([q1f, q2f] => [qf], S4),
    StateFlux([q2u] => [qs], S5)
]

# Total flow calculation
q_flux = HydroFlux([qf, qs] => [q], exprs=[qf + qs])

# Create buckets and model
mopex3_bucket_1 = HydroBucket(name=:mopex3_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
mopex3_bucket_2 = HydroBucket(name=:mopex3_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
mopex3_bucket_3 = HydroBucket(name=:mopex3_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
mopex3_bucket_4 = HydroBucket(name=:mopex3_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)

mopex3_model = HydroModel(name=:mopex3, components=[mopex3_bucket_1, mopex3_bucket_2, mopex3_bucket_3, mopex3_bucket_4, q_flux])

export mopex3_bucket_1, mopex3_bucket_2, mopex3_bucket_3, mopex3_bucket_4, mopex3_model
end
