@reexport module US1
using ..HydroModels
using ..HydroModelLibrary: EvaporationFlux, SaturationFlux, InterceptionFlux, BaseflowFlux

@variables S1 S2 pet prcp temp eusei eusveg eusbs esatveg esatbs rg se qse qss
@parameters alpha_ei m Smax fc alpha_ss

funcs_1 = [
    InterceptionFlux((in = prcp), (p1=alpha_ei,), Val(3), saturation=eusei),
    EvaporationFlux((S1=S1, S2=S2, in=pet), (Smax=Smax, p1=m), Val(8), evaporation=eusveg),
    EvaporationFlux((S1=S1, S2=S2, in=pet), (Smax=Smax, p1=m), Val(9), evaporation=eusbs),
    EvaporationFlux((S1=S2, in=pet), (Smax=S1 + S2, p1=m), Val(10), evaporation=esatveg),
    EvaporationFlux((S1=S2, in=pet), (Smax=S1 + S2, p1=m), Val(5), evaporation=esatbs),
    SaturationFlux((S=S1, in=prcp), (Smax=fc * (Smax - S2),), Val(1), saturation=rg),
    HydroFlux([S1] => [se], [Smax], exprs=[max(0, S1 - Smax)]),
]

dfuncs_1 = [StateFlux([prcp] => [eusei, eusveg, eusbs, rg, se], S1)]

funcs_2 = [
    SaturationFlux((S=S2, in=rg + se), (Smax=Smax,), Val(1), saturation=rg),
    BaseflowFlux((S=S2,), (p1=alpha_ss,), Val(1), saturation=rg),
]
dfuncs_2 = [StateFlux([rg, rs] => [esatveg, esatbs, qse, qss], S2)]

us1_soil_bucket = HydroBucket(name=:us1_soil_bucket, funcs=funcs_1, dfuncs=dfuncs_1)
us1_zone_bucket = HydroBucket(name=:us1_zone_bucket, funcs=funcs_2, dfuncs=dfuncs_2)
us1_model = HydroModel(name=:us1, components=[us1_soil_bucket, us1_zone_bucket])

export us1_soil_bucket, us1_zone_bucket, us1_model
end