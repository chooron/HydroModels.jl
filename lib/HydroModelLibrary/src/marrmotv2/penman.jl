@reexport module Penman
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3
@variables prcp pet temp
@variables ea qex u1 q12 et u2 q
@parameters smax phi gam k1

funcs_1 = [
    EvaporationFlux((S=S1, in=pet), NamedTuple(), Val(1), evaporation=ea),
    SaturationFlux((S=S1, in=prcp), (Smax=smax,), Val(1), saturation=qex),
]

dfuncs_1 = [
    StateFlux([prcp] => [ea, qex], S1),
]

funcs_2 = [
    HydroFlux([qex] => [u1, q12], [phi], exprs=[phi * qex, (1 - phi) * qex]),
    EvaporationFlux((S1=Inf, S2=S1, in=pet), (p1=gam, Smin=0.01), Val(16), evaporation=et),
    SaturationFlux((S=S2, in=q12), NamedTuple(), Val(9), saturation=u2),
    BaseflowFlux((S=S3,), (p1=k1,), Val(1), baseflow=q),
]

dfuncs_2 = [
    StateFlux([et, u2] => [q12], S2),
]

funcs_3 = [
    SaturationFlux((S=S2, in=q12), NamedTuple(), Val(9), saturation=u2),
    BaseflowFlux((S=S3,), (p1=k1,), Val(1), baseflow=q),
]

dfuncs_3 = [
    StateFlux([u1, u2] => [q], S3)
]

penman_bucket_1 = HydroBucket(name=:penman_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
penman_bucket_2 = HydroBucket(name=:penman_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
penman_bucket_3 = HydroBucket(name=:penman_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
penman_model = HydroModel(name=:penman, components=[penman_bucket_1, penman_bucket_2, penman_bucket_3])

export penman_bucket_1, penman_bucket_2, penman_bucket_3, penman_model
end
