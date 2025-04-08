@reexport module Collie1
using ..HydroModels
using ..HydroModelLibrary: EvaporationFlux, SaturationFlux

@variables S pet prcp evaporation q
@parameters Smax

funcs = [
    EvaporationFlux((S=S, pet=pet), (Smax=Smax,), Val(7)),
    SaturationFlux((S=S, in=prcp), (Smax=Smax,), Val(1), saturation=q)
]

dfuncs = [StateFlux([prcp] => [evaporation, q], S)]

collie_bucket = HydroBucket(name=:collie1_bucket, funcs=funcs, dfuncs=dfuncs)

collie_model = HydroModel(name=:collie1, components=[collie_bucket])
export collie_bucket, collie_model
end