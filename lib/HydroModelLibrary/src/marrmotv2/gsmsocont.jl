@reexport module Gsmsocont
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4 S5 S6
@variables prcp pet temp
@variables pice pices picer mis pirs piri mii qis qii pni pnis pnir mnis peq peff pinf et qsl qqu

# Basic parameters
@parameters fice [description = "Fraction of catchment covered by glacier [-]"]
@parameters t0 [description = "Threshold temperature for snowfall [oC]"]
@parameters asnow [description = "Degree-day factor for snow melt [mm/oC/d]"]
@parameters tm [description = "Threshold temperature for snow melt [oC]"]
@parameters ks [description = "Runoff coeficient for snow melt on glacier [d-1]"]
@parameters aice [description = "Threshold temperature for ice melt [oC]"]
@parameters ki [description = "Runoff coeficient for ice melt on glacier [d-1]"]
@parameters a [description = "Maximum soil moisture storage [mm]"]
@parameters x [description = "Evaporation non-linearity [-]"]
@parameters y [description = "Infiltration non-linearity [-]"]
@parameters ksl [description = "Runoff coefficient for baseflow [d-1]"]
@parameters beta [description = "Runoff coefficient for quick flow [mm^(4/3)/d]"]

# Bucket 1: Snow on ice storage
funcs_1 = [
    SplitFlux((in=prcp,), (p1=fice,), Val(1), split=pice),
    SnowfallFlux((in=pice, T=temp), (tt=t0,), Val(1), snowfall=pices),
    MeltFlux((S=S1, T=temp), (ddf=asnow, tt=tm), Val(1), melt=mis)
]

dfuncs_1 = [
    StateFlux([pices] => [mis], S1)
]

# Bucket 2: Snow melt on ice storage
funcs_2 = [
    RainfallFlux((in=pice, T=temp), (tt=t0,), Val(1), rainfall=picer),
    SaturationFlux((in=picer, S=S1), (Smax=0.01,), Val(9), saturation=pirs),
    BaseflowFlux((S=S2,), (p1=ks,), Val(1), baseflow=qis)
]

dfuncs_2 = [
    StateFlux([pirs, mis] => [qis], S2)
]

# Bucket 3: Ice melt storage
funcs_3 = [
    EffectiveFlux((in=picer,), (p1=pirs,), Val(1), effective=piri),
    MeltFlux((S=S1, T=temp), (ddf=aice, tt=tm, Smax=0.01), Val(3), melt=mii),
    BaseflowFlux((S=S3,), (p1=ki,), Val(1), baseflow=qii)
]

dfuncs_3 = [
    StateFlux([piri, mii] => [qii], S3)
]

# Bucket 4: Snow on soil storage
funcs_4 = [
    SplitFlux((in=prcp,), (p1=1-fice,), Val(1), split=pni),
    SnowfallFlux((in=pni, T=temp), (tt=t0,), Val(1), snowfall=pnis),
    RainfallFlux((in=pni, T=temp), (tt=t0,), Val(1), rainfall=pnir),
    MeltFlux((S=S4, T=temp), (ddf=asnow, tt=tm), Val(1), melt=mnis)
]

dfuncs_4 = [
    StateFlux([pnis] => [mnis], S4)
]

# Bucket 5: Soil moisture storage
funcs_5 = [
    EffectiveFlux((in=pnir+mnis,), NamedTuple(), Val(1), effective=peq),
    InfiltrationFlux((S=S5, in=peq), (p1=1, p2=y, Smax=a), Val(6), infiltration=peff),
    EffectiveFlux((in=peq,), (p1=peff,), Val(1), effective=pinf),
    EvaporationFlux((S=S5, in=pet), (p1=1, p2=x, Smax=a), Val(19), evaporation=et),
    BaseflowFlux((S=S5,), (p1=ksl,), Val(1), baseflow=qsl)
]

dfuncs_5 = [
    StateFlux([pinf] => [et, qsl], S5)
]

# Bucket 6: Quick flow storage
funcs_6 = [
    InterflowFlux((S=S6,), (p1=beta, p2=5/3), Val(3), interflow=qqu)
]

dfuncs_6 = [
    StateFlux([peff] => [qqu], S6)
]

# Total flow calculation
q_flux = HydroFlux([qis, qii, qsl, qqu] => [q], exprs=[qis + qii + qsl + qqu])

# Create buckets and model
gsmsocont_bucket_1 = HydroBucket(name=:gsmsocont_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
gsmsocont_bucket_2 = HydroBucket(name=:gsmsocont_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
gsmsocont_bucket_3 = HydroBucket(name=:gsmsocont_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
gsmsocont_bucket_4 = HydroBucket(name=:gsmsocont_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)
gsmsocont_bucket_5 = HydroBucket(name=:gsmsocont_bucket_5, funcs=funcs_5, dfuncs=dfuncs_5)
gsmsocont_bucket_6 = HydroBucket(name=:gsmsocont_bucket_6, funcs=funcs_6, dfuncs=dfuncs_6)

gsmsocont_model = HydroModel(name=:gsmsocont, components=[gsmsocont_bucket_1, gsmsocont_bucket_2,
                                                         gsmsocont_bucket_3, gsmsocont_bucket_4,
                                                         gsmsocont_bucket_5, gsmsocont_bucket_6, q_flux])

export gsmsocont_bucket_1, gsmsocont_bucket_2, gsmsocont_bucket_3
export gsmsocont_bucket_4, gsmsocont_bucket_5, gsmsocont_bucket_6
export gsmsocont_model
end
