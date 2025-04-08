@reexport module VIC
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3
@variables prcp pet temp t
@variables imax ei peff iex qie inf et1 qex1 pc et2 qex2 qb
@parameters ibar [description="Mean interception capacity [mm]"]
@parameters idelta [description="Seasonal interception change as fraction of mean [-]"]
@parameters ishift [description="Maximum interception peak timing [-]"]
@parameters b [description="Infiltration excess shape parameter [-]"]
@parameters k1 [description="Percolation time parameter [d-1]"]
@parameters c1 [description="Percolation non-linearity parameter [-]"]
@parameters k2 [description="Baseflow time parameter [d-1]"]
@parameters c2 [description="Baseflow non-linearity parameter [-]"]
@parameters fsm [description="Fraction of total storage allocated to soil moisture [-]"]
@parameters stot [description="Total storage capacity [mm]"]

# Derived parameters
smmax = fsm * stot      # Maximum soil moisture capacity [mm]
gwmax = (1 - fsm) * stot # Maximum groundwater storage [mm]
tmax = 365.25           # Length of one growing cycle [d]

# Bucket 1: Interception store (S1)
funcs_1 = [
    PhenologyFlux((t=t,), (p1=ibar, p2=idelta, p3=ishift, p4=tmax), Val(2), phenology=imax),
    EvaporationFlux((S=S1, in=pet), (Smax=imax,), Val(7), evaporation=ei),
    InterceptionFlux((S=S1, in=prcp), (Smax=imax,), Val(1), interception=peff),
    ExcessFlux((S=S1,), (Smax=imax,), Val(1), excess=iex)
]

dfuncs_1 = [
    StateFlux([prcp] => [ei, peff, iex], S1)
]

# Bucket 2: Soil moisture store (S2)
funcs_2 = [
    SaturationFlux((S=S2, in=peff + iex), (Smax=smmax, p1=b), Val(2), saturation=qie),
    HydroFlux([peff, iex, qie] => [inf], exprs=[peff + iex - qie]),
    EvaporationFlux((S=S2, in=pet - ei), (Smax=smmax,), Val(7), evaporation=et1),
    SaturationFlux((S=S2, in=inf), (Smax=smmax,), Val(1), saturation=qex1),
    PercolationFlux((S=S2,), (p1=k1, p2=c1, Smax=smmax), Val(5), percolation=pc)
]

dfuncs_2 = [
    StateFlux([inf] => [et1, qex1, pc], S2)
]

# Bucket 3: Groundwater store (S3)
funcs_3 = [
    EvaporationFlux((S=S3, in=pet - ei - et1), (Smax=gwmax,), Val(7), evaporation=et2),
    SaturationFlux((S=S3, in=pc), (Smax=gwmax,), Val(1), saturation=qex2),
    BaseflowFlux((S=S3,), (p1=k2, p2=c2, Smax=gwmax), Val(5), baseflow=qb)
]

dfuncs_3 = [
    StateFlux([pc] => [et2, qex2, qb], S3)
]

# Create buckets and model
vic_bucket_1 = HydroBucket(name=:vic_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
vic_bucket_2 = HydroBucket(name=:vic_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
vic_bucket_3 = HydroBucket(name=:vic_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
vic_model = HydroModel(name=:vic, components=[vic_bucket_1, vic_bucket_2, vic_bucket_3])

export vic_bucket_1, vic_bucket_2, vic_bucket_3, vic_model
end
