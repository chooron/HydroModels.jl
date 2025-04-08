@reexport module Hillslope
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2
@variables prcp pet
@variables pe ei ea qse qses qseg cap qhgw qhsrf q
@parameters d beta Smax a c kh

funcs_1 = [
    InterceptionFlux((P=prcp,), (p1=d,), Val(2), interception=pe),
    HydroFlux([prcp, pe] => [ei], exprs=[prcp - pe]),
    EvaporationFlux((S=S1, in=pet), NamedTuple(), Val(1), evaporation=ea),
    SaturationFlux((S=S1, in=pe), (Smax=Smax, p1=beta), Val(2), saturation=qse),
    HydroFlux([qse] => [qses, qseg], [a], exprs=[a * qse, (1 - a) * qse]),
    CapillaryFlux((S=S2,), (p1=c,), Val(2), capillary=cap),
    BaseflowFlux((S=S2,), (p1=kh,), Val(1), baseflow=qhgw),
]

dfuncs_1 = [StateFlux([pe, cap] => [ea, qse], S1), StateFlux([qseg] => [cap, qhgw], S2)]

hillslope_bucket = HydroBucket(name=:hillslope_bucket, funcs=funcs_1, dfuncs=dfuncs_1)

uh = UnitHydrograph([qses] => [qhsrf], [th], uhfunc=UHFunction(:UH_3_HALF), solvetype=:SPARSE)
q_flux = HydroFlux([qhsrf, qhgw] => [q], exprs=[qhsrf + qhgw])
hillslope_model = HydroModel(name=:hillslope, components=[hillslope_bucket, uh, q_flux])

export hillslope_bucket, hillslope_model
end