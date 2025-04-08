@reexport module Collie2
using ..HydroModels
using ..HydroModelLibrary: EvaporationFlux, SaturationFlux, InterflowFlux

@variables S pet prcp petb petv eb ev qse qss
@parameters Smax M d sfc a

funcs = [
    HydroFlux([pet] => [petb, petv], [M], exprs=[(1 - M) * pet, M * pet]),
    EvaporationFlux((S=S, in=petb), (Smax=Smax, p1=d,), Val(7), evaporation=eb),
    EvaporationFlux((S=S, in=petv), (Smax=Smax, p1=sfc), Val(3), evaporation=ev),
    SaturationFlux((S=S, in=prcp), (Smax=Smax,), Val(1), saturation=qse),
    # TODO HydroFlux need to be able to support formula input
    InterflowFlux((S=S,), (Smax=Smax * sfc, p1=a,), Val(8), saturation=qss),
]

dfuncs = [StateFlux([prcp] => [eb, ev, qse, qss], S)]

collie2_bucket = HydroBucket(name=:collie2_bucket, funcs=funcs, dfuncs=dfuncs)
collie2_model = HydroModel(name=:collie2, components=[collie2_bucket])

export collie2_bucket, collie2_model
end