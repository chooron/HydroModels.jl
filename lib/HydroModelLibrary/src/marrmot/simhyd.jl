@reexport module Simhyd
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3
@variables prcp pet
@variables ei exc inf int rec et gwf bas srun qt smf
@parameters insc coeff sq smsc sub crak k

# Bucket 1: Interception store (S1)
funcs_1 = [
    EvaporationFlux((S=S1, in=pet), NamedTuple(), Val(1), evaporation=ei),
    InterceptionFlux((S=S1, in=prcp), (Smax=insc,), Val(1), interception=exc),
]

dfuncs_1 = [
    StateFlux([prcp] => [ei, exc], S1)
]

# Bucket 2: Soil moisture store (S2)
funcs_2 = [
    InfiltrationFlux((S=S2, in=exc), (p1=coeff, p2=sq, Smax=smsc), Val(1), infiltration=inf),
    InterflowFlux((S=S2, in=inf), (p1=sub, Smax=smsc), Val(1), interflow=int),
    RechargeFlux((S=S2, in=inf - int), (p1=crak, Smax=smsc), Val(1), recharge=rec),
    EvaporationFlux((S=S2, in=pet), (p1=10, Smax=smsc), Val(2), evaporation=et),
    SaturationFlux((S=S2, in=inf - int - rec), (Smax=smsc,), Val(1), saturation=gwf),
    HydroFlux([inf, int, rec] => [smf], exprs=[inf - int - rec])
]

dfuncs_2 = [
    StateFlux([smf] => [et, gwf], S2)
]

# Bucket 3: Groundwater store (S3)
funcs_3 = [
    BaseflowFlux((S=S3,), (p1=k,), Val(1), baseflow=bas)
]

dfuncs_3 = [
    StateFlux([rec, gwf] => [bas], S3)
]

# Surface runoff and total flow calculations
funcs_4 = [
    HydroFlux([exc, inf] => [srun], exprs=[exc - inf]),
    HydroFlux([srun, int, bas] => [qt], exprs=[srun + int + bas])
]

# Create buckets and model
simhyd_bucket_1 = HydroBucket(name=:simhyd_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
simhyd_bucket_2 = HydroBucket(name=:simhyd_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
simhyd_bucket_3 = HydroBucket(name=:simhyd_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
simhyd_model = HydroModel(name=:simhyd, components=[simhyd_bucket_1, simhyd_bucket_2, simhyd_bucket_3, funcs_4])

export simhyd_bucket_1, simhyd_bucket_2, simhyd_bucket_3, simhyd_model
end
