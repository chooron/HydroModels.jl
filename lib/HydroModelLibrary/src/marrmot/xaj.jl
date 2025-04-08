@reexport module Xaj
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4
@variables prcp pet temp
@variables rb pi e r rs ri rg qs qi qg

# Basic parameters
@parameters aim [description = "Fraction impervious area [-]"]
@parameters a [description = "Tension water distribution inflection parameter [-]"]
@parameters b [description = "Tension water distribution shape parameter [-]"]
@parameters stot [description = "Total soil moisture storage (W+S) [mm]"]
@parameters fwmx [description = "Fraction of Stot that is Wmax [-]"]
@parameters flm [description = "Fraction of wmax that is LM [-]"]
@parameters c [description = "Fraction of LM for second evaporation change [-]"]
@parameters ex [description = "Free water distribution shape parameter [-]"]
@parameters ki [description = "Free water interflow parameter [d-1]"]
@parameters kg [description = "Free water baseflow parameter [d-1]"]
@parameters ci [description = "Interflow time coefficient [d-1]"]
@parameters cg [description = "Baseflow time coefficient [d-1]"]

# Calculate auxiliary parameters
wmax = fwmx*stot      # Maximum tension water depth [mm]
smax = (1-fwmx)*stot  # Maximum free water depth [mm]
lm = flm*wmax         # Tension water threshold for evaporation change [mm]

# Bucket 1: Tension water storage (W)
funcs_1 = [
    SplitFlux((in=prcp), (p1=1-aim,), Val(1), split=pi),
    EvaporationFlux((S=S1, in=pet), (p1=lm, p2=c), Val(21), evaporation=e),
    SaturationFlux((S=S1, in=pi), (p1=a, p2=b, Smax=wmax), Val(14), saturation=r)
]

dfuncs_1 = [
    StateFlux([pi] => [e, r], S1)
]

# Bucket 2: Free water storage (S)
funcs_2 = [
    SaturationFlux((S=S2, in=r), (p1=ex, Smax=smax), Val(2), saturation=rs),
    SaturationFlux((S=S2,), (p1=ex, p2=ki, Smax=smax), Val(2), saturation=ri),
    SaturationFlux((S=S2,), (p1=ex, p2=kg, Smax=smax), Val(2), saturation=rg)
]

dfuncs_2 = [
    StateFlux([r] => [rs, ri, rg], S2)
]

# Bucket 3: Interflow storage
funcs_3 = [
    InterflowFlux((S=S3,), (p1=ci,), Val(5), interflow=qi)
]

dfuncs_3 = [
    StateFlux([ri] => [qi], S3)
]

# Bucket 4: Groundwater storage
funcs_4 = [
    BaseflowFlux((S=S4,), (p1=cg,), Val(1), baseflow=qg)
]

dfuncs_4 = [
    StateFlux([rg] => [qg], S4)
]

# Surface runoff calculation
funcs_s = [
    SplitFlux((in=prcp), (p1=aim,), Val(1), split=rb),
    HydroFlux([rb, rs] => [qs], exprs=[rb + rs])
]

# Create buckets and model
xaj_bucket_1 = HydroBucket(name=:xaj_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
xaj_bucket_2 = HydroBucket(name=:xaj_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
xaj_bucket_3 = HydroBucket(name=:xaj_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
xaj_bucket_4 = HydroBucket(name=:xaj_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)
xaj_model = HydroModel(name=:xaj, components=[xaj_bucket_1, xaj_bucket_2, xaj_bucket_3, xaj_bucket_4, funcs_s])

export xaj_bucket_1, xaj_bucket_2, xaj_bucket_3, xaj_bucket_4, xaj_model
end
