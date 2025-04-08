@reexport module Alpine1
using ..HydroModels
using ..HydroModelLibrary: SnowfallFlux, RainfallFlux, MeltFlux, EvaporationFlux, SaturationFlux, InterceptionFlux, BaseflowFlux

@variables S1 S2 pet prcp temp ps pr qn ea qse qss
@parameters tt ddf Smax tc

funcs_1 = [
    ifelse(temp > tt, 0, prcp),
    SnowfallFlux((P=prcp, T=temp), (p1=tt,), Val(1), output=ps),
    RainfallFlux((P=prcp, T=temp), (p1=tt,), Val(1), output=pr),
    MeltFlux((S=S1, T=temp), (p1=ddf, p2=tt), Val(1), output=qn)
]

dfuncs_1 = [StateFlux([ps] => [qn], S1)]

funcs_2 = [
    EvaporationFlux((S=S2, in=pet), NamedTuple(), Val(1), evaporation=ea),
    SaturationFlux((S=S2, in=qn + prcp), (Smax=Smax,), Val(1), saturation=qse),
    BaseflowFlux((S=S2,), (p1=tc,), Val(1), saturation=qss),
]
dfuncs_2 = [StateFlux([pr, qn] => [ea, qse, qss], S2)]

alpine1_snow_bucket = HydroBucket(name=:alpine1_snow_bucket, funcs=funcs_1, dfuncs=dfuncs_1)
alpine1_soil_bucket = HydroBucket(name=:alpine1_soil_bucket, funcs=funcs_2, dfuncs=dfuncs_2)
alpine1_model = HydroModel(name=:alpine1, components=[alpine1_snow_bucket, alpine1_soil_bucket])

export alpine1_snow_bucket, alpine1_soil_bucket, alpine1_model
end