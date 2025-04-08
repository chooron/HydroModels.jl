@reexport module Smar
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4 S5 S6
@variables prcp pet temp
@variables pstar estar evap r1 i r2 e1 e2 e3 e4 e5 q1 q2 q3 q4 r3 rg r3star qr qg

# Basic parameters
@parameters h [description = "Maximum fraction of direct runoff [-]"]
@parameters y [description = "Infiltration rate [mm/d]"]
@parameters smax [description = "Maximum soil moisture storage [mm]"]
@parameters c [description = "Evaporation reduction coefficient [-]"]
@parameters g [description = "Groundwater recharge coefficient [-]"]
@parameters kg [description = "Groundwater time parameter [d-1]"]
@parameters n [description = "Number of Nash cascade reservoirs [-]"]
@parameters nk [description = "Routing delay [d]"]

# Bucket 1: First soil moisture layer
funcs_1 = [
    EffectiveFlux((in=prcp,), (p1=pet,), Val(1), effective=pstar),
    EffectiveFlux((in=pet,), (p1=prcp,), Val(1), effective=estar),
    EffectiveFlux((in=pet,), (p1=prcp,), Val(1), effective=evap, exprs=[min(pet, prcp)]),
    SaturationFlux((in=pstar,), (p1=h, S=S1+S2+S3+S4+S5, Smax=smax), Val(6), saturation=r1),
    InfiltrationFlux((in=pstar-r1,), (p1=y,), Val(4), infiltration=i),
    EffectiveFlux((in=pstar-r1,), (p1=i,), Val(1), effective=r2),
    EvaporationFlux((S=S1, in=estar), (p1=c, p2=0), Val(13), evaporation=e1),
    SaturationFlux((S=S1, in=i), (Smax=smax/5,), Val(1), saturation=q1)
]

dfuncs_1 = [
    StateFlux([i] => [e1, q1], S1)
]

# Bucket 2: Second soil moisture layer
funcs_2 = [
    EvaporationFlux((S=S2, S2=S1, in=estar), (p1=c, p2=1, p3=0.1), Val(14), evaporation=e2),
    SaturationFlux((S=S2, in=q1), (Smax=smax/5,), Val(1), saturation=q2)
]

dfuncs_2 = [
    StateFlux([q1] => [e2, q2], S2)
]

# Bucket 3: Third soil moisture layer
funcs_3 = [
    EvaporationFlux((S=S3, S2=S2, in=estar), (p1=c, p2=2, p3=0.1), Val(14), evaporation=e3),
    SaturationFlux((S=S3, in=q2), (Smax=smax/5,), Val(1), saturation=q3)
]

dfuncs_3 = [
    StateFlux([q2] => [e3, q3], S3)
]

# Bucket 4: Fourth soil moisture layer
funcs_4 = [
    EvaporationFlux((S=S4, S2=S3, in=estar), (p1=c, p2=3, p3=0.1), Val(14), evaporation=e4),
    SaturationFlux((S=S4, in=q3), (Smax=smax/5,), Val(1), saturation=q4)
]

dfuncs_4 = [
    StateFlux([q3] => [e4, q4], S4)
]

# Bucket 5: Fifth soil moisture layer
funcs_5 = [
    EvaporationFlux((S=S5, S2=S4, in=estar), (p1=c, p2=4, p3=0.1), Val(14), evaporation=e5),
    SaturationFlux((S=S5, in=q4), (Smax=smax/5,), Val(1), saturation=r3),
    SplitFlux((in=r3,), (p1=g,), Val(1), split=rg),
    SplitFlux((in=r3,), (p1=1-g,), Val(1), split=r3star)
]

dfuncs_5 = [
    StateFlux([q4] => [e5, r3], S5)
]

# Bucket 6: Groundwater storage
funcs_6 = [
    BaseflowFlux((S=S6,), (p1=kg,), Val(1), baseflow=qg)
]

dfuncs_6 = [
    StateFlux([rg] => [qg], S6)
]

# Unit hydrograph routing
uh_flux = UnitHydrograph((in=r1 + r2 + r3star,), (n=n, k=nk/n), Val(6), route=qr)

# Total flow calculation
q_flux = HydroFlux([qr, qg] => [q], exprs=[qr + qg])

# Create buckets and model
smar_bucket_1 = HydroBucket(name=:smar_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
smar_bucket_2 = HydroBucket(name=:smar_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
smar_bucket_3 = HydroBucket(name=:smar_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
smar_bucket_4 = HydroBucket(name=:smar_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)
smar_bucket_5 = HydroBucket(name=:smar_bucket_5, funcs=funcs_5, dfuncs=dfuncs_5)
smar_bucket_6 = HydroBucket(name=:smar_bucket_6, funcs=funcs_6, dfuncs=dfuncs_6)

smar_model = HydroModel(name=:smar, components=[smar_bucket_1, smar_bucket_2, smar_bucket_3, 
                                                     smar_bucket_4, smar_bucket_5, smar_bucket_6, 
                                                     uh_flux, q_flux])

export smar_bucket_1, smar_bucket_2, smar_bucket_3, smar_bucket_4, smar_bucket_5, 
       smar_bucket_6, smar_model
end
