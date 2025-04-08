@reexport module Nam
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4 S5 S6
@variables prcp pet temp
@variables ps pr m eu pn of inf iff dl gw el qo qi qb

# Basic parameters
@parameters cs [description = "Degree-day factor for snowmelt [mm/oC/d]"]
@parameters cif [description = "Runoff coefficient for interflow [d-1]"]
@parameters stot [description = "Maximum total soil moisture depth [mm]"]
@parameters cl1 [description = "Lower zone filling threshold for interflow generation [-]"]
@parameters fl [description = "Fraction of total soil depth that makes up lstar"]
@parameters cof [description = "Runoff coefficient for overland flow [d-1]"]
@parameters cl2 [description = "Lower zone filling threshold for overland flow generation [-]"]
@parameters k0 [description = "Overland flow routing delay [d-1]"]
@parameters k1 [description = "Interflow routing delay [d-1]"]
@parameters kb [description = "Baseflow routing delay [d-1]"]

# Auxiliary parameters
lstar = fl * stot        # Maximum lower zone storage [mm]
ustar = (1 - fl) * stot  # Upper zone maximum storage [mm]

# Bucket 1: Snow storage
funcs_1 = [
    SnowfallFlux((in=prcp, T=temp), (tt=0,), Val(1), snowfall=ps),
    SnowmeltFlux((S=S1, T=temp), (ddf=cs, tt=0), Val(1), snowmelt=m)
]

dfuncs_1 = [
    StateFlux([ps] => [m], S1)
]

# Bucket 2: Upper zone storage
funcs_2 = [
    RainfallFlux((in=prcp, T=temp), (tt=0,), Val(1), rainfall=pr),
    EvaporationFlux((S=S2, in=pet), NamedTuple(), Val(1), evaporation=eu),
    SaturationFlux((S=S2, in=pr+m), (Smax=ustar,), Val(1), saturation=pn),
    InterflowFlux((S=S2, S2=S3), (p1=cif, p2=cl1, Smax=lstar), Val(6), interflow=iff)
]

dfuncs_2 = [
    StateFlux([pr, m] => [eu, iff, pn], S2)
]

# Bucket 3: Lower zone storage
funcs_3 = [
    InterflowFlux((in=pn, S=S3), (p1=cof, p2=cl2, Smax=lstar), Val(6), interflow=of),
    EffectiveFlux((in=pn,), (p1=of,), Val(1), effective=inf),
    SplitFlux((in=inf,), (p1=1-S3/lstar,), Val(1), split=dl),
    SplitFlux((in=inf,), (p1=S3/lstar,), Val(1), split=gw),
    EvaporationFlux((S=S3, S2=S2, in=pet), (Smax=lstar, p1=0.01), Val(15), evaporation=el)
]

dfuncs_3 = [
    StateFlux([dl] => [el], S3)
]

# Bucket 4: Overland flow storage
funcs_4 = [
    InterflowFlux((S=S4,), (p1=k0,), Val(5), interflow=qo)
]

dfuncs_4 = [
    StateFlux([of] => [qo], S4)
]

# Bucket 5: Interflow storage
funcs_5 = [
    InterflowFlux((S=S5,), (p1=k1,), Val(5), interflow=qi)
]

dfuncs_5 = [
    StateFlux([iff] => [qi], S5)
]

# Bucket 6: Groundwater storage
funcs_6 = [
    BaseflowFlux((S=S6,), (p1=kb,), Val(1), baseflow=qb)
]

dfuncs_6 = [
    StateFlux([gw] => [qb], S6)
]

# Total flow calculation
q_flux = HydroFlux([qo, qi, qb] => [q], exprs=[qo + qi + qb])

# Create buckets and model
nam_bucket_1 = HydroBucket(name=:nam_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
nam_bucket_2 = HydroBucket(name=:nam_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
nam_bucket_3 = HydroBucket(name=:nam_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
nam_bucket_4 = HydroBucket(name=:nam_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)
nam_bucket_5 = HydroBucket(name=:nam_bucket_5, funcs=funcs_5, dfuncs=dfuncs_5)
nam_bucket_6 = HydroBucket(name=:nam_bucket_6, funcs=funcs_6, dfuncs=dfuncs_6)

nam_model = HydroModel(name=:nam, components=[nam_bucket_1, nam_bucket_2, nam_bucket_3, 
                                             nam_bucket_4, nam_bucket_5, nam_bucket_6, q_flux])

export nam_bucket_1, nam_bucket_2, nam_bucket_3, nam_bucket_4, nam_bucket_5, nam_bucket_6, nam_model
end
