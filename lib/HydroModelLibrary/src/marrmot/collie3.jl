@reexport module Collie3
using ..HydroModels
using ..HydroModelLibrary: EvaporationFlux, SaturationFlux, InterflowFlux, BaseflowFlux

@variables S1 S2
@variables pet prcp temp eb ev qse qss qsss qsg qt
@parameters S1max sfc a M b lambda

funcs_1 = [
    EvaporationFlux((S=S, in=(1 - M) * pet), (Smax=S1max, p1=d,), Val(7), evaporation=eb),
    EvaporationFlux((S=S, in=M * pet), (Smax=S1max, p1=sfc), Val(3), evaporation=ev),
    SaturationFlux((S=S, in=prcp), (Smax=S1max,), Val(1), saturation=qse),
    InterflowFlux((S=S1,), (Smax=S1max * sfc, p1=a, p2=b), Val(9), saturation=qss),
    HydroFlux([qss] => [qsss], [lambda], exprs=[qss * lambda]),
    BaseflowFlux((S=S2,), (p1=1 / a, p2=1 / b), Val(2), saturation=qsg),
    HydroFlux([qse, qss, qsg] => [qt], exprs=[qse + (1 - lambda) * flux_qss + flux_qsg])
]

dfuncs_1 = [StateFlux([prcp] => [eb, ev, qse, qss], S1)]

funcs_2 = [
    HydroFlux([qss] => [qsss], [lambda], exprs=[qss * lambda]),
    BaseflowFlux((S=S2,), (p1=1 / a, p2=1 / b), Val(2), saturation=qsg),
]

dfuncs_2 = [StateFlux([qsss] => [qsg], S2)]

collie3_bucket_1 = HydroBucket(name=:collie3_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
collie3_bucket_2 = HydroBucket(name=:collie3_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
collie3_model = HydroModel(name=:collie3, components=[collie3_bucket])

export collie3_bucket, collie3_model
end