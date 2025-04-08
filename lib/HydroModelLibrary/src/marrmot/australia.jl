@reexport module Australia
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3
@variables prcp pet temp
@variables eus rg se esat qse qss qr qbf
@parameters sb phi fc alpha_ss beta_ss k_deep alpha_bf beta_bf

# Bucket 1: Upper soil store (S1)
funcs_1 = [
    EvaporationFlux((S=S1, in=pet), (Smax=sb,), Val(7), evaporation=eus),
    SaturationFlux((S=S1, in=prcp), (Smax=(sb-S2)*fc/phi,), Val(1), saturation=rg),
    ExcessFlux((S=S1,), (Smax=(sb-S2)*fc/phi,), Val(1), excess=se)
]

dfuncs_1 = [
    StateFlux([prcp] => [eus, rg, se], S1)
]

# Bucket 2: Saturated store (S2)
funcs_2 = [
    EvaporationFlux((S=S2, in=pet), (Smax=sb,), Val(7), evaporation=esat),
    SaturationFlux((S=S2, in=rg + se), (Smax=sb,), Val(1), saturation=qse),
    InterflowFlux((S=S2,), (p1=alpha_ss, p2=beta_ss), Val(3), interflow=qss),
    RechargeFlux((S=S2,), (p1=k_deep,), Val(3), recharge=qr)
]

dfuncs_2 = [
    StateFlux([rg, se] => [esat, qse, qss, qr], S2)
]

# Bucket 3: Groundwater store (S3)
funcs_3 = [
    InterflowFlux((S=S3,), (p1=alpha_bf, p2=beta_bf), Val(3), interflow=qbf)
]

dfuncs_3 = [
    StateFlux([qr] => [qbf], S3)
]

# Create buckets and model
australia_bucket_1 = HydroBucket(name=:australia_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
australia_bucket_2 = HydroBucket(name=:australia_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
australia_bucket_3 = HydroBucket(name=:australia_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
australia_model = HydroModel(name=:australia, components=[australia_bucket_1, australia_bucket_2, australia_bucket_3])

export australia_bucket_1, australia_bucket_2, australia_bucket_3, australia_model
end
