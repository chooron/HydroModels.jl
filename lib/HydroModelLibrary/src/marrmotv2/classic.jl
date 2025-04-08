@reexport module Classic
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4 S5 S6 S7 S8
@variables prcp pet temp
@variables pp ps pi epx ppx epy ppe q psd psi esx psx esy pse psq pss xq xs pie ei u

# Basic parameters
@parameters fap [description = "Fraction of catchment area that has permeable soils [-]"]
@parameters fdp [description = "Fraction of depth of permeable soil that is store Px [-]"]
@parameters dp [description = "Depth of permeable soil [mm]"]
@parameters cq [description = "Runoff coefficient for permeable soil [d-1]"]
@parameters d1 [description = "Fraction of Ps that infiltrates into semi-permeable soil [-]"]
@parameters tf [description = "Fraction of (1-fap) that is fas [-]"]
@parameters fds [description = "Fraction of depth of semi-permeable soil that is store Sx [-]"]
@parameters ds [description = "Depth of semi-permeable soil [mm]"]
@parameters d2 [description = "Fraction effective precipitation in semi-permeable soils that goes to quick flow [-]"]
@parameters cxq [description = "Quick runoff coefficient for semi-permeable soil [d-1]"]
@parameters cxs [description = "Slow runoff coefficient for semi-permeable soil [d-1]"]
@parameters cu [description = "Runoff coefficient for impermeable soil [d-1]"]

# Auxiliary parameters
fas = (1-fap)*tf      # Fraction of catchment area that has semi-permeable soils [-]
fai = 1-fap-fas       # Fraction of catchment area that has impermeable soils [-]
pxm = fdp*dp         # Depth of store Px [mm]
pym = (1-fdp)*dp     # Depth of store Py [mm]
sxm = fds*ds         # Depth of store Sx [mm]
sym = (1-fds)*ds     # Depth of store Sy [mm]

# Bucket 1: Upper permeable soil storage (Px)
funcs_1 = [
    SplitFlux((in=prcp,), (p1=fap,), Val(1), split=pp),
    EvaporationFlux((S=S1, in=fap*pet), NamedTuple(), Val(1), evaporation=epx),
    SaturationFlux((in=pp, S=S1), (Smax=pxm,), Val(1), saturation=ppx)
]

dfuncs_1 = [
    StateFlux([pp] => [epx, ppx], S1)
]

# Bucket 2: Lower permeable soil storage (Py)
funcs_2 = [
    EvaporationFlux((S=S2+pxm, in=fap*pet-epx), (p1=1.9, p2=0.6523), Val(18), evaporation=epy),
    SaturationFlux((in=ppx, S=S2), (Smax=0.01,), Val(9), saturation=ppe)
]

dfuncs_2 = [
    StateFlux([ppx] => [epy, ppe], S2)
]

# Bucket 3: Permeable soil runoff storage
funcs_3 = [
    BaseflowFlux((S=S3,), (p1=cq,), Val(1), baseflow=q)
]

dfuncs_3 = [
    StateFlux([ppe] => [q], S3)
]

# Bucket 4: Upper semi-permeable soil storage (Sx)
funcs_4 = [
    SplitFlux((in=prcp,), (p1=fas,), Val(1), split=ps),
    SplitFlux((in=ps,), (p1=d1,), Val(1), split=psi),
    EvaporationFlux((S=S4, in=fas*pet), NamedTuple(), Val(1), evaporation=esx),
    SaturationFlux((in=psi, S=S4), (Smax=sxm,), Val(1), saturation=psx)
]

dfuncs_4 = [
    StateFlux([psi] => [esx, psx], S4)
]

# Bucket 5: Lower semi-permeable soil storage (Sy)
funcs_5 = [
    EvaporationFlux((S=S4+sxm, in=fas*pet-esx), (p1=1.9, p2=0.6523), Val(18), evaporation=esy),
    SaturationFlux((in=psx, S=S5), (Smax=0.01,), Val(9), saturation=pse)
]

dfuncs_5 = [
    StateFlux([psx] => [esy, pse], S5)
]

# Bucket 6: Quick flow storage from semi-permeable soil
funcs_6 = [
    SplitFlux((in=ps,), (p1=1-d1,), Val(1), split=psd),
    SplitFlux((in=pse+psd,), (p1=d2,), Val(1), split=psq),
    BaseflowFlux((S=S6,), (p1=cxq,), Val(1), baseflow=xq)
]

dfuncs_6 = [
    StateFlux([psq] => [xq], S6)
]

# Bucket 7: Slow flow storage from semi-permeable soil
funcs_7 = [
    SplitFlux((in=pse+psd,), (p1=1-d2,), Val(1), split=pss),
    BaseflowFlux((S=S7,), (p1=cxs,), Val(1), baseflow=xs)
]

dfuncs_7 = [
    StateFlux([pss] => [xs], S7)
]

# Bucket 8: Impermeable soil storage
funcs_8 = [
    SplitFlux((in=prcp,), (p1=fai,), Val(1), split=pi),
    EffectiveFlux((in=fai*pi,), (p1=0.5,), Val(1), effective=pie),
    EffectiveFlux((in=pi,), (p1=pie,), Val(1), effective=ei),
    BaseflowFlux((S=S8,), (p1=cu,), Val(1), baseflow=u)
]

dfuncs_8 = [
    StateFlux([pie] => [u], S8)
]

# Total flow calculation
q_flux = HydroFlux([q, xq, xs, u] => [qt], exprs=[q + xq + xs + u])

# Create buckets and model
classic_bucket_1 = HydroBucket(name=:classic_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
classic_bucket_2 = HydroBucket(name=:classic_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
classic_bucket_3 = HydroBucket(name=:classic_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
classic_bucket_4 = HydroBucket(name=:classic_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)
classic_bucket_5 = HydroBucket(name=:classic_bucket_5, funcs=funcs_5, dfuncs=dfuncs_5)
classic_bucket_6 = HydroBucket(name=:classic_bucket_6, funcs=funcs_6, dfuncs=dfuncs_6)
classic_bucket_7 = HydroBucket(name=:classic_bucket_7, funcs=funcs_7, dfuncs=dfuncs_7)
classic_bucket_8 = HydroBucket(name=:classic_bucket_8, funcs=funcs_8, dfuncs=dfuncs_8)

classic_model = HydroModel(name=:classic, components=[classic_bucket_1, classic_bucket_2,
                                                     classic_bucket_3, classic_bucket_4,
                                                     classic_bucket_5, classic_bucket_6,
                                                     classic_bucket_7, classic_bucket_8, q_flux])

export classic_bucket_1, classic_bucket_2, classic_bucket_3, classic_bucket_4
export classic_bucket_5, classic_bucket_6, classic_bucket_7, classic_bucket_8
export classic_model
end
