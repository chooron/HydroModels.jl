@reexport module GSFB
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3
@variables prcp pet temp
@variables ea qs f qb dp qdr
@parameters c ndc smax emax frate b dpf sdrmax

# Bucket 1: Surface store (S1)
funcs_1 = [
    EvaporationFlux((S=S1, in=pet), (p1=emax, p2=ndc, Smax=smax), Val(20), evaporation=ea),
    SaturationFlux((S=S1, in=prcp), (Smax=smax,), Val(1), saturation=qs),
    InterflowFlux((S=S1,), (p1=frate, p2=ndc * smax), Val(11), interflow=f),
    RechargeFlux((S=S3, S2=S1), (p1=c, p2=ndc * smax), Val(5), recharge=qdr)
]

dfuncs_1 = [
    StateFlux([prcp, qdr] => [ea, qs, f], S1)
]

# Bucket 2: Subsurface store (S2)
funcs_2 = [
    BaseflowFlux((S=S2,), (p1=b * dpf, p2=sdrmax), Val(9), baseflow=qb),
    BaseflowFlux((S=S2,), (p1=(1 - b) * dpf,), Val(1), baseflow=dp)
]

dfuncs_2 = [
    StateFlux([f] => [qb, dp], S2)
]

# Bucket 3: Deep store (S3)
funcs_3 = [
    RechargeFlux((S=S3, S2=S1), (p1=c, p2=ndc * smax), Val(5), recharge=qdr)
]

dfuncs_3 = [
    StateFlux([dp] => [qdr], S3)
]

# Create buckets and model
gsfb_bucket_1 = HydroBucket(name=:gsfb_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
gsfb_bucket_2 = HydroBucket(name=:gsfb_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
gsfb_bucket_3 = HydroBucket(name=:gsfb_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
gsfb_model = HydroModel(name=:gsfb, components=[gsfb_bucket_1, gsfb_bucket_2, gsfb_bucket_3])

export gsfb_bucket_1, gsfb_bucket_2, gsfb_bucket_3, gsfb_model
end
