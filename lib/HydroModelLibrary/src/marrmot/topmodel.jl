@reexport module TopModel
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2
@variables prcp pet temp
@variables qof peff ea qex qv qb
@parameters suzmax st kd q0 f chi phi

lambda = 3 + chi * phi

funcs_1 = [
    SaturationFlux((S=S2, in=prcp), (p1=chi, p2=phi, p3=3, p4=lambda, p5=f), Val(7), saturation=qof),
    HydroFlux([qof, prcp] => [peff], exprs=[prcp - qof]),
    EvaporationFlux((S=S1, in=pet), (Smax=suzmax, p1=st), Val(3), evaporation=ea),
    SaturationFlux((S=S1, in=peff), (Smax=suzmax,), Val(1), saturation=qex),
    InterflowFlux((S=S1,), (p1=kd, p2=st * suzmax, p3=suzmax - st * suzmax), Val(10), interflow=qv),
    BaseflowFlux((S=S2,), (p1=q0, p2=f), Val(4), baseflow=qb),
]

dfuncs_1 = [StateFlux([peff] => [ea, qex, qv], S1), StateFlux([qb] => [qv], S2)]

topmodel_bucket = HydroBucket(name=:topmodel_bucket, funcs=funcs_1, dfuncs=dfuncs_1)

q_flux = HydroFlux([qof, qex, qb] => [q], exprs=[qof + qex + qb])
topmodel_model = HydroModel(name=:topmodel, components=[topmodel_bucket, q_flux])

export topmodel_bucket, topmodel_model
end