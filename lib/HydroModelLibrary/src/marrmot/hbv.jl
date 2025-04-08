@reexport module Hbv
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4 S5 S2old
@variables prcp pet temp t
@variables sf refr melt rf in se cf ea r q0 perc q1 qt

# Basic parameters
@parameters tt [description = "Middle of snow-rain interval [oC]"]
@parameters tti [description = "Interval length of rain-snow spectrum [oC]"]
@parameters ttm [description = "Threshold temperature for snowmelt [oC]"]
@parameters cfr [description = "Coefficient of refreezing of melted snow [-]"]
@parameters cfmax [description = "Degree-day factor of snowmelt and refreezing [mm/oC/d]"]
@parameters whc [description = "Maximum water holding content of snow pack [-]"]
@parameters cflux [description = "Maximum rate of capillary rise [mm/d]"]
@parameters fc [description = "Maximum soil moisture storage [mm]"]
@parameters lp [description = "Wilting point as fraction of FC [-]"]
@parameters beta [description = "Non-linearity coefficient of upper zone recharge [-]"]
@parameters k0 [description = "Runoff coefficient from upper zone [d-1]"]
@parameters alpha [description = "Non-linearity coefficient of runoff from upper zone [-]"]
@parameters perc [description = "Maximum rate of percolation to lower zone [mm/d]"]
@parameters k1 [description = "Runoff coefficient from lower zone [d-1]"]
@parameters maxbas [description = "Flow routing delay [d]"]

# Bucket 1: Snow storage
funcs_1 = [
    SnowfallFlux((in=prcp, T=temp), (tt=tt, tti=tti), Val(2), snowfall=sf),
    RefreezeFlux((S=S2, T=temp), (p1=cfr, p2=cfmax, tt=ttm), Val(1), refreeze=refr),
    SnowmeltFlux((S=S1, T=temp), (ddf=cfmax, tt=ttm), Val(1), snowmelt=melt)
]

dfuncs_1 = [
    StateFlux([sf, refr] => [melt], S1)
]

# Bucket 2: Snow liquid storage
funcs_2 = [
    RainfallFlux((in=prcp, T=temp), (tt=tt, tti=tti), Val(2), rainfall=rf),
    InfiltrationFlux((in=rf+melt, S=S2), (Smax=whc*S1,), Val(3), infiltration=in),
    ExcessFlux((S=S2old, in=whc*S1), NamedTuple(), Val(1), excess=se)
]

dfuncs_2 = [
    StateFlux([rf, melt] => [refr, in, se], S2)
]

# Bucket 3: Soil moisture storage
funcs_3 = [
    CapillaryFlux((S=S3, S2=S4), (Smax=fc, p1=cflux), Val(1), capillary=cf),
    EvaporationFlux((S=S3, in=pet), (Smax=fc, p1=lp), Val(3), evaporation=ea),
    RechargeFlux((S=S3, in=in+se), (Smax=fc, p1=beta), Val(2), recharge=r)
]

dfuncs_3 = [
    StateFlux([in, se, cf] => [ea, r], S3)
]

# Bucket 4: Upper zone storage
funcs_4 = [
    InterflowFlux((S=S4,), (p1=k0, p2=alpha), Val(2), interflow=q0),
    PercolationFlux((S=S4,), (p1=perc,), Val(1), percolation=perc)
]

dfuncs_4 = [
    StateFlux([r] => [cf, q0, perc], S4)
]

# Bucket 5: Lower zone storage
funcs_5 = [
    BaseflowFlux((S=S5,), (p1=k1,), Val(1), baseflow=q1)
]

dfuncs_5 = [
    StateFlux([perc] => [q1], S5)
]

# Total flow calculation with routing
route_flux = RouteFlux((in=q0+q1,), (nlag=maxbas,), Val(4), route=qt)

# Create buckets and model
hbv_bucket_1 = HydroBucket(name=:hbv_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
hbv_bucket_2 = HydroBucket(name=:hbv_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
hbv_bucket_3 = HydroBucket(name=:hbv_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
hbv_bucket_4 = HydroBucket(name=:hbv_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)
hbv_bucket_5 = HydroBucket(name=:hbv_bucket_5, funcs=funcs_5, dfuncs=dfuncs_5)

hbv_model = HydroModel(name=:hbv, components=[hbv_bucket_1, hbv_bucket_2, hbv_bucket_3, 
                                             hbv_bucket_4, hbv_bucket_5, route_flux])

export hbv_bucket_1, hbv_bucket_2, hbv_bucket_3, hbv_bucket_4, hbv_bucket_5, hbv_model
end
