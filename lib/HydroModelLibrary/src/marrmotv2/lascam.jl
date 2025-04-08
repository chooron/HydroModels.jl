@reexport module Lascam
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3
@variables prcp pet temp
@variables phiss phic fss pg ei qse pc qie qsse fa qsie ef rf ea1 ea2 qa ra qb eb

# Basic parameters
@parameters af [description = "Catchment-scale infiltration parameter [mm/d]"]
@parameters bf [description = "Catchment-scale infiltration non-linearity parameter [-]"]
@parameters stot [description = "Total catchment storage [mm]"]
@parameters xa [description = "Fraction of Stot that is Amax [-]"]
@parameters xf [description = "Fraction of Stot-Amx that is depth Fmax [-]"]
@parameters na [description = "Fraction of Amax that is Amin [-]"]
@parameters ac [description = "Variable contributing area scaling [-]"]
@parameters bc [description = "Variable contributing area non-linearity [-]"]
@parameters ass [description = "Subsurface saturation area scaling [-]"]
@parameters bss [description = "Subsurface saturation area non-linearity [-]"]
@parameters c [description = "Maximum infiltration rate [mm/d]"]
@parameters ag [description = "Interception base parameter [mm/d]"]
@parameters bg [description = "Interception fraction parameter [-]"]
@parameters gf [description = "F-store evaporation scaling [-]"]
@parameters df [description = "F-store evaporation non-linearity [-]"]
@parameters td [description = "Recharge time parameter [d-1]"]
@parameters ab [description = "Groundwater flow scaling [-]"]
@parameters bb [description = "Groundwater flow base rate [mm/d]"]
@parameters ga [description = "A-store evaporation scaling [-]"]
@parameters da [description = "A-store evaporation non-linearity [-]"]
@parameters aa [description = "Subsurface storm flow rate [mm/d]"]
@parameters ba [description = "Subsurface storm flow non-linearity [-]"]
@parameters gb [description = "B-store evaporation scaling [-]"]
@parameters db [description = "B-store evaporation non-linearity [-]"]

# Derived parameters
amax = xa * stot            # Maximum contributing area depth [mm]
fmax = xf * (stot - amax)   # Infiltration depth scaling [mm]
bmax = (1 - xf) * (stot - amax) # Groundwater depth scaling [mm]
amin = na * amax            # Minimum contributing area depth [mm]

# Bucket 1: F-store (Infiltration store)
funcs_1 = [
    AreaFlux((S=S2,), (p1=ass, p2=bss, p3=amin, p4=amax), Val(1), area=phiss),
    AreaFlux((S=S2,), (p1=ac, p2=bc, p3=amin, p4=amax), Val(1), area=phic),
    InfiltrationFlux((S=S3, S2=S1), (p1=af, p2=bf, Smax=bmax, Smax2=fmax), Val(5), infiltration=fss),
    InterceptionFlux((in=prcp,), (p1=bg, p2=ag), Val(5), interception=pg),
    HydroFlux([prcp, pg] => [ei], exprs=[prcp - pg]),
    SaturationFlux((S=S2, in=pg), (p1=ac, p2=bc, p3=amin, p4=amax), Val(11), saturation=qse),
    InfiltrationFlux((in=pg - qse,), (p1=c,), Val(4), infiltration=pc),
    HydroFlux([pg, qse, pc] => [qie], exprs=[pg - qse - pc]),
    EvaporationFlux((S=S1, in=pet), (p1=gf, p2=df, Smax=fmax), Val(19), evaporation=ef),
    RechargeFlux((S=S1,), (p1=td,), Val(3), recharge=rf)
]

dfuncs_1 = [
    StateFlux([fa] => [ef, rf], S1)
]

# Bucket 2: A-store (Variable contributing area store)
funcs_2 = [
    SaturationFlux((S=S2,), (p1=phiss, p2=phic, in=pc), Val(12), saturation=qsse),
    InfiltrationFlux((in=pc,), (p1=fss, p2=phiss, p3=phic), Val(4), infiltration=fa),
    HydroFlux([pc, fa, qsse] => [qsie], exprs=[pc - fa - qsse]),
    EvaporationFlux((S=S2, in=phic * pet), NamedTuple(), Val(1), evaporation=ea1),
    EvaporationFlux((S=S2, in=pet), (p1=ga, p2=da, Smax=amax), Val(19), evaporation=ea2),
    SaturationFlux((S=S2,), (p1=aa, p2=ba, p3=amin, p4=amax), Val(11), saturation=qa),
    RechargeFlux((S=S2,), (p1=phic, p2=fss), Val(4), recharge=ra)
]

dfuncs_2 = [
    StateFlux([qsse, qsie, qb] => [ea1, ea2, ra, qa], S2)
]

# Bucket 3: B-store (Groundwater store)
funcs_3 = [
    BaseflowFlux((S=S3,), (p1=bb, p2=ab, Smax=bmax), Val(8), baseflow=qb),
    EvaporationFlux((S=S3, in=pet), (p1=gb, p2=db, Smax=bmax), Val(19), evaporation=eb)
]

dfuncs_3 = [
    StateFlux([rf, ra] => [eb, qb], S3)
]

# Create buckets and model
lascam_bucket_1 = HydroBucket(name=:lascam_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
lascam_bucket_2 = HydroBucket(name=:lascam_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
lascam_bucket_3 = HydroBucket(name=:lascam_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
lascam_model = HydroModel(name=:lascam, components=[lascam_bucket_1, lascam_bucket_2, lascam_bucket_3])

export lascam_bucket_1, lascam_bucket_2, lascam_bucket_3, lascam_model
end
