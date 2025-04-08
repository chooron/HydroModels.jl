@reexport module Tank2
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4 S5
@variables prcp pet temp
@variables t1 t2 y1 y2 y3 y4 y5 e1 e2 e3 e4 f12 f23 f34

# Basic parameters
@parameters a0 [description = "Time parameter for drainage 1>2 [d-1]"]
@parameters b0 [description = "Time parameter for drainage 2>3 [d-1]"]
@parameters c0 [description = "Time parameter for drainage 3>4 [d-1]"]
@parameters a1 [description = "Time parameter for surface runoff 1 [d-1]"]
@parameters fa [description = "Fraction of a1 that is a2 [-]"]
@parameters fb [description = "Fraction of a2 that is b1 [-]"]
@parameters fc [description = "Fraction of b1 that is c1 [-]"]
@parameters fd [description = "Fraction of c1 that is d1 [-]"]
@parameters st [description = "Maximum soil depth (sum of runoff thresholds) [mm]"]
@parameters f2 [description = "Fraction of st-sm1 that is added to sm1 to find threshold t2 [-]"]
@parameters f1 [description = "Fraction of st-t2 that is added to t2 to find threshold 1 [-]"]
@parameters f3 [description = "Fraction of st-t1-sm2 that consitutes threshold 3 [-]"]
@parameters k1 [description = "Base rate of capillary rise [mm/d]"]
@parameters k2 [description = "Base rate of soil moisture exchange [mm/d]"]
@parameters z1 [description = "Fraction st that is sm1 [-]"]
@parameters z2 [description = "Fraction of st-t1 that is sm2 [-]"]

# Auxiliary parameters (calculated once during initialization)
sm1 = z1*st                # Size of primary soil moisture store, threshold before F12 starts [mm]
t2_val = sm1+f2*(st-sm1)  # Threshold before surface runoff Y2 starts [mm]
t1_val = t2_val+f1*(st-t2_val)  # Thresold before surface runoff Y1 starts [mm]
sm2 = z2*(st-t1_val)      # Size of secondary soil moisture store S5 [mm]
t3_val = f3*(st-t1_val-sm2)  # Threshold before intermediate runoff starts [mm]
t4_val = st-t1_val-sm2-t3_val  # Threshold before sub-base runoff starts [mm]
a2_val = fa*a1            # Time parameter for surface runoff 2 [d-1]
b1_val = fb*a2_val        # Time parameter for intermediate runoff 1 [d-1]
c1_val = fc*b1_val        # Time parameter for sub-base runoff 1 [d-1]
d1_val = fd*c1_val        # Time parameter for base runoff 1 [d-1]

# Bucket 1: Primary soil moisture storage
funcs_1 = [
    CapillaryFlux((S=S1, S2=S2), (p1=k1, Smax=sm1), Val(3), capillary=t1),
    ExchangeFlux((S=S1, S2=S5), (p1=k2, Smax=sm1, S2max=sm2), Val(2), exchange=t2),
    InterflowFlux((S=S1,), (p1=a1, Smax=t1_val), Val(8), interflow=y1),
    InterflowFlux((S=S1,), (p1=a2_val, Smax=t2_val), Val(8), interflow=y2),
    EvaporationFlux((S=S1, in=pet), NamedTuple(), Val(1), evaporation=e1),
    InterflowFlux((S=S1,), (p1=a0, Smax=sm1), Val(8), interflow=f12)
]

dfuncs_1 = [
    StateFlux([prcp, t1] => [t2, e1, f12, y1, y2], S1)
]

# Bucket 2: Secondary soil moisture storage
funcs_2 = [
    InterflowFlux((S=S2,), (p1=b1_val, Smax=t3_val), Val(8), interflow=y3),
    EvaporationFlux((S=S2, in=max(0,pet-e1)), NamedTuple(), Val(1), evaporation=e2),
    RechargeFlux((S=S2,), (p1=b0,), Val(3), recharge=f23)
]

dfuncs_2 = [
    StateFlux([f12] => [t1, e2, f23, y3], S2)
]

# Bucket 3: Intermediate storage
funcs_3 = [
    InterflowFlux((S=S3,), (p1=c1_val, Smax=t4_val), Val(8), interflow=y4),
    EvaporationFlux((S=S3, in=max(0,pet-e1-e2)), NamedTuple(), Val(1), evaporation=e3),
    RechargeFlux((S=S3,), (p1=c0,), Val(3), recharge=f34)
]

dfuncs_3 = [
    StateFlux([f23] => [e3, f34, y4], S3)
]

# Bucket 4: Base storage
funcs_4 = [
    BaseflowFlux((S=S4,), (p1=d1_val,), Val(1), baseflow=y5),
    EvaporationFlux((S=S4, in=max(0,pet-e1-e2-e3)), NamedTuple(), Val(1), evaporation=e4)
]

dfuncs_4 = [
    StateFlux([f34] => [e4, y5], S4)
]

# Bucket 5: Secondary soil moisture store
dfuncs_5 = [
    StateFlux([t2] => [], S5)
]

# Total flow calculation
q_flux = HydroFlux([y1, y2, y3, y4, y5] => [q], exprs=[y1 + y2 + y3 + y4 + y5])

# Create buckets and model
tank2_bucket_1 = HydroBucket(name=:tank2_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
tank2_bucket_2 = HydroBucket(name=:tank2_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
tank2_bucket_3 = HydroBucket(name=:tank2_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
tank2_bucket_4 = HydroBucket(name=:tank2_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)
tank2_bucket_5 = HydroBucket(name=:tank2_bucket_5, funcs=[], dfuncs=dfuncs_5)

tank2_model = HydroModel(name=:tank2, components=[tank2_bucket_1, tank2_bucket_2, tank2_bucket_3, 
                                                 tank2_bucket_4, tank2_bucket_5, q_flux])

export tank2_bucket_1, tank2_bucket_2, tank2_bucket_3, tank2_bucket_4, tank2_bucket_5, tank2_model
end
