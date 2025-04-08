@reexport module Prms
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4 S5 S6 S7
@variables prcp pet temp
@variables ps pr pim psm pby pin ptf m mim msm sas sro inf pc excs sep qres gad ras bas snk ein eim ea et

# Basic parameters
@parameters tt [description = "Temperature threshold for snowfall and melt [oC]"]
@parameters ddf [description = "Degree-day factor for snowmelt [mm/oC/d]"]
@parameters alpha [description = "Fraction of rainfall on soil moisture going to interception [-]"]
@parameters beta [description = "Fraction of catchment where rain goes to soil moisture [-]"]
@parameters stor [description = "Maximum interception capcity [mm]"]
@parameters retip [description = "Maximum impervious area storage [mm]"]
@parameters fscn [description = "Fraction of SCX where SCN is located [-]"]
@parameters scx [description = "Maximum contributing fraction area to saturation excess flow [-]"]
@parameters flz [description = "Fraction of total soil moisture that is the lower zone [-]"]
@parameters stot [description = "Total soil moisture storage [mm]: REMX+SMAX"]
@parameters cgw [description = "Constant drainage to deep groundwater [mm/d]"]
@parameters resmax [description = "Maximum flow routing reservoir storage [mm]"]
@parameters k1 [description = "Groundwater drainage coefficient [d-1]"]
@parameters k2 [description = "Groundwater drainage non-linearity [-]"]
@parameters k3 [description = "Interflow coefficient 1 [d-1]"]
@parameters k4 [description = "Interflow coefficient 2 [mm-1 d-1]"]
@parameters k5 [description = "Baseflow coefficient [d-1]"]
@parameters k6 [description = "Groundwater sink coefficient [d-1]"]

# Auxiliary parameters
scn = fscn * scx      # Minimum contributing fraction area to saturation excess flow [-]
remx = (1-flz) * stot # Maximum upper soil moisture storage [mm]
smax = flz * stot     # Maximum lower soil moisture storage [mm]

# Bucket 1: Snow storage
funcs_1 = [
    SnowfallFlux((in=prcp, T=temp), (tt=tt,), Val(1), snowfall=ps),
    MeltFlux((S=S1, T=temp), (ddf=ddf, tt=tt), Val(1), melt=m)
]

dfuncs_1 = [
    StateFlux([ps] => [m], S1)
]

# Bucket 2: Interception storage
funcs_2 = [
    RainfallFlux((in=prcp, T=temp), (tt=tt,), Val(1), rainfall=pr),
    SplitFlux((in=pr,), (p1=beta,), Val(1), split=psm),
    SplitFlux((in=psm,), (p1=alpha,), Val(1), split=pin),
    InterceptionFlux((in=pin, S=S2), (Smax=stor,), Val(1), interception=ptf),
    EvaporationFlux((S=S2, in=beta*pet), NamedTuple(), Val(1), evaporation=ein)
]

dfuncs_2 = [
    StateFlux([pin] => [ein, ptf], S2)
]

# Bucket 3: Impervious area storage
funcs_3 = [
    SplitFlux((in=pr,), (p1=1-beta,), Val(1), split=pim),
    SplitFlux((in=m,), (p1=1-beta,), Val(1), split=mim),
    SaturationFlux((in=pim+mim, S=S3), (Smax=retip,), Val(1), saturation=sas),
    EvaporationFlux((S=S3, in=(1-beta)*pet), NamedTuple(), Val(1), evaporation=eim)
]

dfuncs_3 = [
    StateFlux([pim, mim] => [eim, sas], S3)
]

# Bucket 4: Upper soil zone storage
funcs_4 = [
    SplitFlux((in=psm,), (p1=1-alpha,), Val(1), split=pby),
    SplitFlux((in=m,), (p1=beta,), Val(1), split=msm),
    SaturationFlux((S=S4, in=msm+ptf+pby), (p1=scn, p2=scx, Smax=remx), Val(8), saturation=sro),
    EffectiveFlux((in=msm+ptf+pby,), (p1=sro,), Val(1), effective=inf),
    SaturationFlux((in=inf, S=S4), (Smax=remx,), Val(1), saturation=pc),
    EvaporationFlux((S=S4, in=pet-ein-eim), (Smax=remx,), Val(7), evaporation=ea)
]

dfuncs_4 = [
    StateFlux([inf] => [ea, pc], S4)
]

# Bucket 5: Lower soil zone storage
funcs_5 = [
    SaturationFlux((in=pc, S=S5), (Smax=smax,), Val(1), saturation=excs),
    RechargeFlux((in=excs,), (p1=cgw,), Val(7), recharge=sep),
    EffectiveFlux((in=excs,), (p1=sep,), Val(1), effective=qres),
    EvaporationFlux((S=S5, S2=S4, in=pet-ein-eim-ea), (Smax=smax, p1=pet-ein-eim), Val(15), evaporation=et)
]

dfuncs_5 = [
    StateFlux([pc] => [et, excs], S5)
]

# Bucket 6: Subsurface reservoir
funcs_6 = [
    RechargeFlux((S=S6,), (p1=k2, p2=resmax, p3=k1), Val(2), recharge=gad),
    InterflowFlux((S=S6,), (p1=k3, p2=k4), Val(4), interflow=ras)
]

dfuncs_6 = [
    StateFlux([qres] => [gad, ras], S6)
]

# Bucket 7: Groundwater reservoir
funcs_7 = [
    BaseflowFlux((S=S7,), (p1=k5,), Val(1), baseflow=bas),
    BaseflowFlux((S=S7,), (p1=k6,), Val(1), baseflow=snk)
]

dfuncs_7 = [
    StateFlux([sep, gad] => [bas, snk], S7)
]

# Total flow calculation
q_flux = HydroFlux([sas, sro, ras, bas] => [q], exprs=[sas + sro + ras + bas])

# Create buckets and model
prms_bucket_1 = HydroBucket(name=:prms_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
prms_bucket_2 = HydroBucket(name=:prms_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
prms_bucket_3 = HydroBucket(name=:prms_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
prms_bucket_4 = HydroBucket(name=:prms_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)
prms_bucket_5 = HydroBucket(name=:prms_bucket_5, funcs=funcs_5, dfuncs=dfuncs_5)
prms_bucket_6 = HydroBucket(name=:prms_bucket_6, funcs=funcs_6, dfuncs=dfuncs_6)
prms_bucket_7 = HydroBucket(name=:prms_bucket_7, funcs=funcs_7, dfuncs=dfuncs_7)

prms_model = HydroModel(name=:prms, components=[prms_bucket_1, prms_bucket_2, prms_bucket_3,
                                               prms_bucket_4, prms_bucket_5, prms_bucket_6,
                                               prms_bucket_7, q_flux])

export prms_bucket_1, prms_bucket_2, prms_bucket_3, prms_bucket_4
export prms_bucket_5, prms_bucket_6, prms_bucket_7, prms_model
end
