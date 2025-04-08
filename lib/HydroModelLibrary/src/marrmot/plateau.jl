@reexport module Plateau
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2
@variables prcp pet temp
@variables pe ei pie pin et r cap qpgw qpieo q
@parameters fmax dp sumax lp p tp c kp

funcs_1 = [
    InterceptionFlux((P=prcp,), (p1=dp,), Val(2), interception=pe),
    HydroFlux([pe, prcp] => [ei], exprs=[prcp - pe]),
    InfiltrationFlux((in=pe,), (p1=fmax,), Val(4), infiltration=pin),
    HydroFlux([pin, pe] => [pie], exprs=[pe - pin]),
    EvaporationFlux((S=S1, in=pet), (Smax=suzmax, p1=p, p2=lp), Val(4), evaporation=et),
    CapillaryFlux((S=S2,), (p1=c,), Val(2), capillary=cap),
    SaturationFlux((S=S1, in=pin + cap), (Smax=suzmax,), Val(1), saturation=r),
]
dfuncs_1 = [StateFlux([pin, cap] => [et, r], S1)]

funcs_2 = [
    BaseflowFlux((S=S2,), (p1=kp,), Val(1), baseflow=qpgw),
]
dfuncs_2 = [StateFlux([r] => [qpgw, cap], S2)]

plateau_bucket_1 = HydroBucket(name=:plateau_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
plateau_bucket_2 = HydroBucket(name=:plateau_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)

uh = UnitHydrograph([pie] => [qpieo], [tp], uhfunc=UHFunction(:UH_3_HALF), solvetype=:SPARSE)
q_flux = HydroFlux([qpgw, qpieo] => [q], exprs=[qpgw + qpieo])
plateau_model = HydroModel(name=:plateau, components=[plateau_bucket_1, plateau_bucket_2, uh, q_flux])

export plateau_bucket, plateau_model
end