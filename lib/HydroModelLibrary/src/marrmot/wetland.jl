@reexport module Wetland
using ..HydroModels
using ..HydroModelLibrary: EvaporationFlux, SaturationFlux, InterceptionFlux, BaseflowFlux

@variables S pet prcp evap pe qs qg
@parameters Smax d beta k

funcs = [
    InterceptionFlux((S=S,), (p1=d,), Val(2), interception=pe),
    EvaporationFlux((S=S, in=pet), NamedTuple(), Val(1), evaporation=evap),
    SaturationFlux((S=S, in=pe), (Smax=Smax, p1=beta), Val(2), saturation=qs),
    BaseflowFlux((S=S,), (p1=k,), Val(1), saturation=qg),
]

dfuncs = [StateFlux([pe] => [evap, qs, qg], S)]

wetland_bucket = HydroBucket(name=:wetland_bucket, funcs=funcs, dfuncs=dfuncs)
wetland_model = HydroModel(name=:wetland, components=[wetland_bucket])

export wetland_bucket, wetland_model
end