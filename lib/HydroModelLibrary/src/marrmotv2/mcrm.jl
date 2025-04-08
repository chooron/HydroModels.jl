@reexport module Mcrm
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4 S5
@variables prcp pet temp
@variables ec qt qr er qn qd qp qb uib uob qic qoc

# Basic parameters
@parameters smax [description = "Maximum interception storage [mm]"]
@parameters cmax [description = "Maximum fraction of area contributing to rapid runoff [-]"]
@parameters ct [description = "Fraction of cmax that is the minimum contributing area c0 [-]"]
@parameters c1 [description = "Shape parameter for rapid flow distribution [mm-1]"]
@parameters ce [description = "Shape parameter for evaporation [mm-1]"]
@parameters dsurp [description = "Threshold for direct runoff [mm]"]
@parameters kd [description = "Direct runoff time parameter [d-1]"]
@parameters gamd [description = "Direct runoff flow non-linearity [-]"]
@parameters qpmax [description = "Maximum percolation rate [mm/d]"]
@parameters kg [description = "Groundwater time parameter [d-1]"]
@parameters tau [description = "Routing delay [d]"]
@parameters sbf [description = "Maximum routing store depth [mm]"]
@parameters kcr [description = "Channel flow time parameter [d-1]"]
@parameters gamcr [description = "Channel flow non-linearity [-]"]
@parameters kor [description = "Out-of-bank flow time parameter [d-1]"]
@parameters gamor [description = "Out-of-bank flow non-linearity [-]"]

# Auxiliary parameters
c0 = ct * cmax  # Minimum fraction of area contributing to rapid runoff [-]

# Bucket 1: Interception storage
funcs_1 = [
    EvaporationFlux((S=S1, in=pet), NamedTuple(), Val(1), evaporation=ec),
    InterceptionFlux((S=S1, in=prcp), (Smax=smax,), Val(1), interception=qt)
]

dfuncs_1 = [
    StateFlux([prcp] => [ec, qt], S1)
]

# Bucket 2: Soil moisture storage
funcs_2 = [
    SaturationFlux((in=qt,), (p1=cmax, p2=c0, p3=c1, S=S2), Val(10), saturation=qr),
    EvaporationFlux((S=S2, in=pet - ec), (p1=ce,), Val(17), evaporation=er),
    EffectiveFlux((in=qt,), (p1=qr,), Val(1), effective=qn),
    InterflowFlux((S=S2,), (p1=kd, p2=dsurp, p3=gamd), Val(9), interflow=qd),
    PercolationFlux((S=S2,), (p1=qpmax, p2=dsurp), Val(6), percolation=qp)
]

dfuncs_2 = [
    StateFlux([qn] => [er, qd, qp], S2)
]

# Bucket 3: Groundwater storage
funcs_3 = [
    BaseflowFlux((S=S3,), (p1=kg, p2=1.5), Val(7), baseflow=qb)
]

dfuncs_3 = [
    StateFlux([qp] => [qb], S3)
]

# Unit hydrograph routing
uh_flux = UnitHydrograph((in=qr + qd + qb,), (nlag=tau,), Val(7), route=uib)

# Bucket 4: Channel routing storage
funcs_4 = [
    SaturationFlux((S=S4, in=uib), (Smax=sbf,), Val(1), saturation=uob),
    RoutingFlux((S=S4,), (p1=kcr, p2=gamcr, p3=3 / 4), Val(1), routing=qic)
]

dfuncs_4 = [
    StateFlux([uib] => [uob, qic], S4)
]

# Bucket 5: Out-of-bank storage
funcs_5 = [
    RoutingFlux((S=S5,), (p1=kor, p2=gamor, p3=3 / 4), Val(1), routing=qoc)
]

dfuncs_5 = [
    StateFlux([uob] => [qoc], S5)
]


# Total flow calculation
q_flux = HydroFlux([qic, qoc] => [q], exprs=[qic + qoc])

# Create buckets and model
mcrm_bucket_1 = HydroBucket(name=:mcrm_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
mcrm_bucket_2 = HydroBucket(name=:mcrm_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
mcrm_bucket_3 = HydroBucket(name=:mcrm_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
mcrm_bucket_4 = HydroBucket(name=:mcrm_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)
mcrm_bucket_5 = HydroBucket(name=:mcrm_bucket_5, funcs=funcs_5, dfuncs=dfuncs_5)

mcrm_model = HydroModel(name=:mcrm, components=[mcrm_bucket_1, mcrm_bucket_2, mcrm_bucket_3, uh_flux,
    mcrm_bucket_4, mcrm_bucket_5, q_flux])

export mcrm_bucket_1, mcrm_bucket_2, mcrm_bucket_3, mcrm_bucket_4, mcrm_bucket_5, mcrm_model
end
