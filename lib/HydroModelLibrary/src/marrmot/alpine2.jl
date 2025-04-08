@reexport module Alpine2
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 pet prcp temp ps pr qn ea qseqin qbf q
@parameters tt ddf Smax cfc tcin tcbf

funcs_1 = [
    SnowfallFlux((P=prcp, T=temp), (p1=tt,), Val(1), interception=ps),
    RainfallFlux((P=prcp, T=temp), (p1=tt,), Val(1), interception=pr),
    MeltFlux((S=S1, T=temp), (p1=ddf, p2=tt), Val(1), interception=qn)
]

dfuncs_1 = [StateFlux([ps] => [qn], S1)]

funcs_2 = [
    EvaporationFlux((S=S2, in=pet), NamedTuple(), Val(1), evaporation=ea),
    SaturationFlux((S=S2, in=qn + prcp), (Smax=Smax,), Val(1), saturation=qse),
    InterflowFlux((S=S2,), (p1=tcin, p2=cfc * Smax), Val(8), saturation=qin),
    BaseflowFlux((S=S2,), (p1=tcbf,), Val(1), saturation=qbf),
]
dfuncs_2 = [StateFlux([pr, qn] => [ea, qse, qin, qbf], S2)]

alpine2_bucket_1 = HydroBucket(name=:alpine2_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
alpine2_bucket_2 = HydroBucket(name=:alpine2_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)

q_flux = HydroFlux([qse, qin, qbf] => [q], exprs=[qse + qin - qbf])
alpine2_model = HydroModel(name=:alpine2, components=[alpine2_bucket_1, alpine2_bucket_2])

export alpine2_bucket_1, alpine2_bucket_2, alpine2_model
end