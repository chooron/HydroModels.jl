@reexport module Mopex1
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4
@variables prcp pet temp
@variables et1 q1f qw et2 q2u qf qs

@parameters s1max [description = "Maximum soil moisture storage [mm]"]
@parameters tw [description = "Groundwater leakage time [d-1]"]
@parameters tu [description = "Slow flow routing response time [d-1]"]
@parameters se [description = "Root zone storage capacity [mm]"]
@parameters tc [description = "Mean residence time [d-1]"]

# Bucket 1: Surface store (S1)
funcs_1 = [
    EvaporationFlux((S=S1, in=pet), (Smax=s1max,), Val(7), evaporation=et1),
    SaturationFlux((S=S1, in=prcp), (Smax=s1max,), Val(1), saturation=q1f),
    RechargeFlux((S=S1,), (p1=tw,), Val(3), recharge=qw)
]

dfuncs_1 = [
    StateFlux([prcp] => [et1, q1f, qw], S1)
]

# Bucket 2: Root zone store (S2)
funcs_2 = [
    EvaporationFlux((S=S2, in=pet), (Smax=se,), Val(7), evaporation=et2),
    BaseflowFlux((S=S2,), (p1=tu,), Val(1), baseflow=q2u)
]

dfuncs_2 = [
    StateFlux([qw] => [et2, q2u], S2)
]

# Bucket 3: Fast response store (S3)
funcs_3 = [
    BaseflowFlux((S=S3,), (p1=tc,), Val(1), baseflow=qf)
]

dfuncs_3 = [
    StateFlux([q1f] => [qf], S3)
]

# Bucket 4: Slow response store (S4)
funcs_4 = [
    BaseflowFlux((S=S4,), (p1=tc,), Val(1), baseflow=qs)
]

dfuncs_4 = [
    StateFlux([q2u] => [qs], S4)
]

# Create buckets and model
mopex1_bucket_1 = HydroBucket(name=:mopex1_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
mopex1_bucket_2 = HydroBucket(name=:mopex1_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
mopex1_bucket_3 = HydroBucket(name=:mopex1_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
mopex1_bucket_4 = HydroBucket(name=:mopex1_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)
mopex1_model = HydroModel(name=:mopex1, components=[mopex1_bucket_1, mopex1_bucket_2, mopex1_bucket_3, mopex1_bucket_4])

export mopex1_bucket_1, mopex1_bucket_2, mopex1_bucket_3, mopex1_bucket_4, mopex1_model
end
