@reexport module Hycymodel
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4 S5 S6
@variables prcp pet temp
@variables rc rg eic qie qis rt eis rs rn esu re qin esb qb qh qc qt ec

# Basic parameters
@parameters c [description = "Fraction area that is channel [-]"]
@parameters imax [description = "Maximum total interception storage [mm]"]
@parameters a [description = "Fraction stem/trunk interception [-]"]
@parameters fi2 [description = "Fraction of total interception that is trunk/stem interception [mm]"]
@parameters kin [description = "Infiltration runoff coefficient [d-1]"]
@parameters d50 [description = "Soil depth where 50% of area contributes to effective flow [mm]"]
@parameters fd16 [description = "Fraction of D50 that is D16 [-]"]
@parameters sbc [description = "Soil depth where evaporation rate starts to decline [mm]"]
@parameters kb [description = "Baseflow runoff coefficient [d-1]"]
@parameters pb [description = "Baseflow non-linearity [-]"]
@parameters kh [description = "Hillslope runoff coefficient [d-1]"]
@parameters kc [description = "Channel runoff coefficient [d-1]"]

# Auxiliary parameters
i1max = (1-fi2)*imax  # Maximum canopy interception [mm]
i2max = fi2*imax      # Maximum trunk/stem interception [mm]
d16 = fd16*d50       # Soil depth where 16% of area contributes to effective flow [mm]

# Bucket 1: Canopy interception storage
funcs_1 = [
    SplitFlux((in=prcp,), (p1=1-c,), Val(1), split=rg),
    EvaporationFlux((S=S1, in=(1-c)*pet), NamedTuple(), Val(1), evaporation=eic),
    InterceptionFlux((in=rg, S=S1), (Smax=i1max,), Val(1), interception=qie)
]

dfuncs_1 = [
    StateFlux([rg] => [eic, qie], S1)
]

# Bucket 2: Trunk/stem interception storage
funcs_2 = [
    SplitFlux((in=qie,), (p1=a,), Val(1), split=qis),
    SplitFlux((in=qie,), (p1=1-a,), Val(1), split=rt),
    EvaporationFlux((S=S2, in=(1-c)*pet), NamedTuple(), Val(1), evaporation=eis),
    InterceptionFlux((in=qis, S=S2), (Smax=i2max,), Val(1), interception=rs)
]

dfuncs_2 = [
    StateFlux([qis] => [eis, rs], S2)
]

# Bucket 3: Upper soil storage
funcs_3 = [
    EffectiveFlux((in=rt+rs,), NamedTuple(), Val(1), effective=rn),
    EvaporationFlux((S=S3, in=(1-c)*pet), NamedTuple(), Val(1), evaporation=esu),
    SaturationFlux((S=S3, in=rn), (d50=d50, d16=d16), Val(13), saturation=re),
    RechargeFlux((S=S3,), (p1=kin,), Val(3), recharge=qin)
]

dfuncs_3 = [
    StateFlux([rn] => [re, esu, qin], S3)
]

# Bucket 4: Lower soil storage
funcs_4 = [
    EvaporationFlux((S=S4, in=max(0,(1-c)*pet-esu)), (Smax=sbc,), Val(3), evaporation=esb),
    BaseflowFlux((S=S4,), (p1=kb, p2=pb), Val(7), baseflow=qb)
]

dfuncs_4 = [
    StateFlux([qin] => [esb, qb], S4)
]

# Bucket 5: Hillslope storage
funcs_5 = [
    InterflowFlux((S=S5,), (p1=kh, p2=5/3), Val(3), interflow=qh)
]

dfuncs_5 = [
    StateFlux([re] => [qh], S5)
]

# Bucket 6: Channel storage
funcs_6 = [
    SplitFlux((in=prcp,), (p1=c,), Val(1), split=rc),
    InterflowFlux((S=S6,), (p1=kc, p2=5/3), Val(3), interflow=qc)
]

dfuncs_6 = [
    StateFlux([rc] => [qc], S6)
]

# Total flow calculation
q_flux = HydroFlux([qb, qh, qc, pet] => [qt, ec], exprs=[
    qb + qh + qc - (qb + qh + qc - qt),  # qt
    (qb + qh + qc) - qt                   # ec
])

# Create buckets and model
hycymodel_bucket_1 = HydroBucket(name=:hycymodel_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
hycymodel_bucket_2 = HydroBucket(name=:hycymodel_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
hycymodel_bucket_3 = HydroBucket(name=:hycymodel_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
hycymodel_bucket_4 = HydroBucket(name=:hycymodel_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)
hycymodel_bucket_5 = HydroBucket(name=:hycymodel_bucket_5, funcs=funcs_5, dfuncs=dfuncs_5)
hycymodel_bucket_6 = HydroBucket(name=:hycymodel_bucket_6, funcs=funcs_6, dfuncs=dfuncs_6)

hycymodel_model = HydroModel(name=:hycymodel, components=[hycymodel_bucket_1, hycymodel_bucket_2, 
                                                         hycymodel_bucket_3, hycymodel_bucket_4, 
                                                         hycymodel_bucket_5, hycymodel_bucket_6, q_flux])

export hycymodel_bucket_1, hycymodel_bucket_2, hycymodel_bucket_3
export hycymodel_bucket_4, hycymodel_bucket_5, hycymodel_bucket_6
export hycymodel_model
end
