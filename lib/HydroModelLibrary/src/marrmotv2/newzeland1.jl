@reexport module Newzeland1
using ..HydroModels
using ..HydroModelLibrary: EvaporationFlux, SaturationFlux, BaseflowFlux, InterflowFlux

@variables S pet prcp veg ebs qse qss qbf
@parameters Smax sfc m a b tcbf

funcs = [
    EvaporationFlux((S=S, in=pet), (Smax=Smax, p1=m, p2=sfc), Val(6), evaporation=veg),
    EvaporationFlux((S=S, in=pet), (Smax=Smax, p1=m), Val(5), evaporation=ebs),
    SaturationFlux((S=S, in=prcp), (Smax=Smax,), Val(1), saturation=qse),
    InterflowFlux((S=S,), (p1=a, p2=sfc * Smax, p3=b), Val(9), saturation=qss),
    BaseflowFlux((S=S,), (p1=tcbf,), Val(1), saturation=qbf),
]

dfuncs = [StateFlux([prcp] => [veg, ebs, qse, qss, qbf], S)]

newzeland1_bucket = HydroBucket(name=:newzeland1_bucket, funcs=funcs, dfuncs=dfuncs)
newzeland1 = HydroModel(name=:newzeland1, components=[newzeland1_bucket])

export newzeland1_bucket, newzeland1
end