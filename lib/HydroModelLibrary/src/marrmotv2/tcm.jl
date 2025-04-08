@reexport module TCM
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4
@variables prcp pet temp
@variables pn en pby pin ea et qex1 qex2 quz a q

@parameters phi [description = "Fraction preferential recharge [-]"]
@parameters rc [description = "Maximum soil moisture depth [mm]"]
@parameters gam [description = "Fraction of Ep reduction with depth [-]"]
@parameters k1 [description = "Runoff coefficient [d-1]"]
@parameters fa [description = "Fraction of average P abstracted [-]"]
@parameters k2 [description = "Runoff coefficient [mm-1 d-1]"]

# Derived parameter (calculated in init)
# ca = fa * mean(P)    # Abstraction rate [mm/day]
# We'll need to implement a way to calculate this dynamically

# Bucket 1: Upper soil store (S1)
funcs_1 = [
    HydroFlux([prcp, pet] => [pn, en], exprs=[max(0.0, prcp - pet), prcp - max(0.0, prcp - pet)]),
    HydroFlux([pn] => [pby, pin], [phi], exprs=[phi * pn, (1 - phi) * pn]),
    EvaporationFlux((S=S1, in=pet), NamedTuple(), Val(1), evaporation=ea),
    SaturationFlux((S=S1, in=pin), (Smax=rc,), Val(1), saturation=qex1)
]

dfuncs_1 = [
    StateFlux([pin] => [ea, qex1], S1)
]

# Bucket 2: Lower soil store (S2)
funcs_2 = [
    EvaporationFlux((S1=Inf, S2=S1, in=pet), (p1=gam, Smin=0.01), Val(16), evaporation=et),
    SaturationFlux((S=S2, in=qex1), (Smax=0.01,), Val(9), saturation=qex2)
]

dfuncs_2 = [
    StateFlux([et, qex2] => [qex1], S2)
]

# Bucket 3: Upper groundwater store (S3)
funcs_3 = [
    BaseflowFlux((S=S3,), (p1=k1,), Val(1), baseflow=quz)
]

dfuncs_3 = [
    StateFlux([qex2, pby] => [quz], S3)
]

# Bucket 4: Lower groundwater store (S4)
funcs_4 = [
    AbstractionFlux((in=1.0,), NamedTuple(), Val(1), abstraction=a),  # ca will be set dynamically
    BaseflowFlux((S=S4,), (p1=k2, p2=0.0), Val(6), baseflow=q)
]

dfuncs_4 = [
    StateFlux([quz] => [a, q], S4)
]

# Create buckets and model
tcm_bucket_1 = HydroBucket(name=:tcm_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
tcm_bucket_2 = HydroBucket(name=:tcm_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
tcm_bucket_3 = HydroBucket(name=:tcm_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
tcm_bucket_4 = HydroBucket(name=:tcm_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)
tcm_model = HydroModel(name=:tcm, components=[tcm_bucket_1, tcm_bucket_2, tcm_bucket_3, tcm_bucket_4])

export tcm_bucket_1, tcm_bucket_2, tcm_bucket_3, tcm_bucket_4, tcm_model
end
