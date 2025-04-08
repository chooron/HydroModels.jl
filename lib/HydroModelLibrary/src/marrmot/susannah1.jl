@reexport module Suannah1
using ..HydroModels
using ..HydroModelLibrary: EvaporationFlux, SaturationFlux, InterceptionFlux, BaseflowFlux, InterflowFlux

@variables S1 S2 pet prcp temp ebs evg qse qss qr qb qt
@parameters sb sfc m a b r

funcs_1 = [
    EvaporationFlux((S1=S1, in=pet), (Smax=sb, p1=m), Val(5), evaporation=ebs),
    EvaporationFlux((S1=S1, in=pet), (Smax=sb, p1=m, p2=sfc), Val(6), evaporation=eveg),
    SaturationFlux((S=S1, in=prcp), (Smax=sb,), Val(1), saturation=qse),
    InterflowFlux((S=S1,), (Smax=sb, p1=sfc, p2=a, p3=b), Val(7), saturation=qss),
]

dfuncs_1 = [StateFlux([prcp] => [ebs, eveg, qse, qss], S1)]

funcs_2 = [
    BaseflowFlux((in=qss,), (p1=r,), Val(1), saturation=qr),
    BaseflowFlux((S=S2,), (p1=a, p2=b), Val(2), saturation=qb),
    HydroFlux([qse, qss, qr, qb] => [qt], exprs=[qse + qss - qr + qb])
]
dfuncs_2 = [StateFlux([qr] => [qb], S2)]

suannah1_soil_bucket = HydroBucket(name=:suannah1_soil_bucket, funcs=funcs_1, dfuncs=dfuncs_1)
suannah1_zone_bucket = HydroBucket(name=:suannah1_zone_bucket, funcs=funcs_2, dfuncs=dfuncs_2)
suannah1_model = HydroModel(name=:suannah1, components=[suannah1_soil_bucket, suannah1_zone_bucket])

export suannah1_soil_bucket, suannah1_zone_bucket, suannah1_model
end