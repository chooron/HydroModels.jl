@reexport module Tank1
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4
@variables prcp pet temp
@variables y1 y2 y3 y4 y5 e1 e2 e3 e4 f12 f23 f34

# Time parameters for drainage and runoff
@parameters a0 [description = "Time parameter for drainage 1>2 [d-1]"]
@parameters b0 [description = "Time parameter for drainage 2>3 [d-1]"]
@parameters c0 [description = "Time parameter for drainage 3>4 [d-1]"]
@parameters a1 [description = "Time parameter for surface runoff 1 [d-1]"]
@parameters fa [description = "Fraction of a1 that is a2 [-]"]
@parameters fb [description = "Fraction of a2 that is b1 [-]"]
@parameters fc [description = "Fraction of b1 that is c1 [-]"]
@parameters fd [description = "Fraction of c1 that is d1 [-]"]

# Storage and threshold parameters
@parameters st [description = "Maximum soil depth (sum of runoff thresholds) [mm]"]
@parameters f2 [description = "Fraction of st that constitutes threshold t2 [-]"]
@parameters f1 [description = "Fraction of st-t2 that is added to t2 to find threshold 1 [-]"]
@parameters f3 [description = "Fraction of st-t1-t2 that constitutes threshold 3 [-]"]

# Threshold calculations
t2 = f2*st                     # Threshold before surface runoff 2 starts [mm]
t1 = t2 + f1*(st-t2)          # Threshold before surface runoff 1 starts [mm]
t3 = f3*(st-t1)               # Threshold before intermediate runoff starts [mm]
t4 = st - t1 - t3             # Threshold before sub-base runoff starts [mm]

# Time parameter calculations
a2 = fa*a1                     # Time parameter for surface runoff 2 [d-1]
b1 = fb*a2                     # Time parameter for intermediate runoff 1 [d-1]
c1 = fc*b1                     # Time parameter for sub-base runoff 1 [d-1]
d1 = fd*c1                     # Time parameter for base runoff 1 [d-1]

# Bucket 1: Surface layer
funcs_1 = [
    InterflowFlux((S=S1, in=prcp), (p1=a1, Smax=t1), Val(8), interflow=y1),
    InterflowFlux((S=S1, in=prcp), (p1=a2, Smax=t2), Val(8), interflow=y2),
    RechargeFlux((S=S1,), (p1=a0,), Val(3), recharge=f12),
    EvaporationFlux((S=S1, in=pet), NamedTuple(), Val(1), evaporation=e1)
]

dfuncs_1 = [
    StateFlux([prcp] => [y1, y2, f12, e1], S1)
]

# Bucket 2: Intermediate layer 1
funcs_2 = [
    InterflowFlux((S=S2, in=f12), (p1=b1, Smax=t3), Val(8), interflow=y3),
    RechargeFlux((S=S2,), (p1=b0,), Val(3), recharge=f23),
    EvaporationFlux((S=S2, in=pet), NamedTuple(), Val(1), evaporation=e2)
]

dfuncs_2 = [
    StateFlux([f12] => [y3, f23, e2], S2)
]

# Bucket 3: Intermediate layer 2
funcs_3 = [
    InterflowFlux((S=S3, in=f23), (p1=c1, Smax=t4), Val(8), interflow=y4),
    RechargeFlux((S=S3,), (p1=c0,), Val(3), recharge=f34),
    EvaporationFlux((S=S3, in=pet), NamedTuple(), Val(1), evaporation=e3)
]

dfuncs_3 = [
    StateFlux([f23] => [y4, f34, e3], S3)
]

# Bucket 4: Base layer
funcs_4 = [
    BaseflowFlux((S=S4,), (p1=d1,), Val(1), baseflow=y5),
    EvaporationFlux((S=S4, in=pet), NamedTuple(), Val(1), evaporation=e4)
]

dfuncs_4 = [
    StateFlux([f34] => [y5, e4], S4)
]

# Total flow calculation
q_flux = HydroFlux([y1, y2, y3, y4, y5] => [q], exprs=[y1 + y2 + y3 + y4 + y5])

# Create buckets and model
tank1_bucket_1 = HydroBucket(name=:tank1_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
tank1_bucket_2 = HydroBucket(name=:tank1_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
tank1_bucket_3 = HydroBucket(name=:tank1_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
tank1_bucket_4 = HydroBucket(name=:tank1_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)
tank1_model = HydroModel(name=:tank1, components=[tank1_bucket_1, tank1_bucket_2, tank1_bucket_3, tank1_bucket_4, q_flux])

export tank1_bucket_1, tank1_bucket_2, tank1_bucket_3, tank1_bucket_4, tank1_model
end
