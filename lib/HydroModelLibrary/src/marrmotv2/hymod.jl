@reexport module Hymod
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4 S5
@variables prcp pet temp
@variables ea pe pf ps qf1 qf2 qf3 qs

@parameters smax [description = "Maximum soil moisture storage [mm]"]
@parameters b [description = "Soil depth distribution parameter [-]"]
@parameters a [description = "Runoff distribution fraction [-]"]
@parameters kf [description = "Fast runoff coefficient [d-1]"]
@parameters ks [description = "Slow runoff coefficient [d-1]"]

# Bucket 1: Soil moisture storage
funcs_1 = [
    EvaporationFlux((S=S1, in=pet), (Smax=smax,), Val(7), evaporation=ea),
    SaturationFlux((S=S1, in=prcp), (p1=b, Smax=smax), Val(2), saturation=pe),
    SplitFlux((in=pe), (p1=a,), Val(1), split=pf),
    SplitFlux((in=pe), (p1=1-a,), Val(1), split=ps)
]

dfuncs_1 = [
    StateFlux([prcp] => [ea, pe], S1)
]

# Bucket 2: Fast reservoir 1
funcs_2 = [
    BaseflowFlux((S=S2,), (p1=kf,), Val(1), baseflow=qf1)
]

dfuncs_2 = [
    StateFlux([pf] => [qf1], S2)
]

# Bucket 3: Fast reservoir 2
funcs_3 = [
    BaseflowFlux((S=S3,), (p1=kf,), Val(1), baseflow=qf2)
    BaseflowFlux((S=S4,), (p1=kf,), Val(1), baseflow=qf3)
    BaseflowFlux((S=S5,), (p1=ks,), Val(1), baseflow=qs)
]

dfuncs_3 = [
    StateFlux([qf1] => [qf2], S3)
    StateFlux([qf2] => [qf3], S4)
    StateFlux([ps] => [qs], S5)
]

# Create buckets and model
hymod_bucket_1 = HydroBucket(name=:hymod_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
hymod_bucket_2 = HydroBucket(name=:hymod_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
hymod_bucket_3 = HydroBucket(name=:hymod_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)

# Total flow calculation
q_flux = HydroFlux([qf3, qs] => [q], exprs=[qf3 + qs])

hymod_model = HydroModel(name=:hymod, components=[hymod_bucket_1, hymod_bucket_2, hymod_bucket_3, q_flux])

export hymod_bucket_1, hymod_bucket_2, hymod_bucket_3, hymod_model
end
