@reexport module Suannah2
using ..HydroModels
using ..HydroModelLibrary: EvaporationFlux, SaturationFlux, ExcessFlux, InterceptionFlux, BaseflowFlux, InterflowFlux

@variables S1 S2 pet prcp temp eus rg se esat qse qss qr qt
@parameters sb phi fc r c d

funcs_1 = [
    EvaporationFlux((S1=S1, in=pet), (Smax=sb,), Val(7), evaporation=eus),
    SaturationFlux((S=S1, in=prcp), (Smax=(sb - S2) * fc / phi,), Val(1), saturation=rg),
    ExcessFlux((S=S1,), (Smax=(sb - S2) * fc / phi,), Val(1), saturation=se),
]

dfuncs_1 = [StateFlux([prcp] => [eus, rg, se], S1)]

funcs_2 = [
    EvaporationFlux((S1=S2, in=pet), (Smax=sb,), Val(7), evaporation=esat),
    SaturationFlux((S=S2, in=rg + se), (Smax=sb,), Val(1), saturation=qse),
    InterflowFlux((S=S2,), (p1=(1 - r) * c, p2=d), Val(3), saturation=qss),
    InterflowFlux((S=S2,), (p1=r * c, p2=d), Val(3), saturation=qr),
    HydroFlux([qse, qss] => [qt], exprs=[qse + qss])
]
dfuncs_2 = [StateFlux([rg, se] => [esat, qse, qss, qr], S2)]

suannah2_bucket_1 = HydroBucket(name=:suannah2_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
suannah2_bucket_2 = HydroBucket(name=:suannah2_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
suannah2_model = HydroModel(name=:suannah2, components=[suannah2_bucket_1, suannah2_bucket_2])

export suannah2_bucket_1, suannah2_bucket_2, suannah2_model
end