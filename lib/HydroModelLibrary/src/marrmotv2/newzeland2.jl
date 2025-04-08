@reexport module Newzeland2
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2
@variables prcp pet temp
@variables eint qtf veg ebs qse qss qbf qt q
@parameters s1max s2max sfc m a b tcbf

funcs_1 = [
    EvaporationFlux((S=S1, in=pet), NamedTuple(), Val(1), evaporation=eint),
    InterceptionFlux((S=S1, in=prcp), (Smax=s1max,), Val(1), interception=qtf),
]
dfuncs_1 = [StateFlux([prcp, cap] => [eint, qtf], S1)]

funcs_2 = [
    EvaporationFlux((S=S2, in=pet), (Smax=s2max, p1=m, p2=sfc), Val(6), evaporation=veg),
    EvaporationFlux((S=S2, in=pet), (Smax=s2max, p1=m), Val(5), evaporation=ebs),
    SaturationFlux((S=S2, in=qtf), (Smax=s2max,), Val(1), saturation=qse),
    InterflowFlux((S=S2,), (p1=a, p2=sfc * s2max, p3=b), Val(9), interflow=qss),
    BaseflowFlux((S=S2,), (p1=tcbf,), Val(1), baseflow=qbf),
    HydroFlux([qse, qss, qbf] => [qt], exprs=[qse + qss + qbf])
]
dfuncs_2 = [StateFlux([qtf] => [veg, ebs, qse, qss, qbf], S2)]

newzeland2_bucket_1 = HydroBucket(name=:newzeland2_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
newzeland2_bucket_2 = HydroBucket(name=:newzeland2_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)

uh = UnitHydrograph([qt] => [q], [d], uhfunc=UHFunction(:UH_4_FULL), solvetype=:SPARSE)
newzeland2_model = HydroModel(name=:newzeland2, components=[newzeland2_bucket_1, newzeland2_bucket_2, uh])

export newzeland2_bucket, newzeland2_model
end