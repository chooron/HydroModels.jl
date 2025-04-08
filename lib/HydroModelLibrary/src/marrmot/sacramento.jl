@reexport module Sacramento
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4 S5
@variables prcp pet
@variables qdir peff ru euztw twexu qsur qint euzfw pc pctw elztw twexl twexlp twexls pcfwp pcfws rlp rls qbfp qbfs

# Basic parameters
@parameters pctim [description = "Fraction impervious area [-]"]
@parameters smax [description = "Maximum total storage depth [mm]"]
@parameters f1 [description = "Fraction of smax that is Maximum upper zone tension water storage [-]"]
@parameters f2 [description = "Fraction of smax-uztwm that is Maximum upper zone free water storage [-]"]
@parameters kuz [description = "Interflow runoff coefficient [d-1]"]
@parameters rexp [description = "Base percolation rate non-linearity factor [-]"]
@parameters f3 [description = "Fraction of smax-uztwm-uzfwm that is Maximum lower zone tension water storage [-]"]
@parameters f4 [description = "Fraction of smax-uztwm-uzfwm-lztwm that is Maximum lower zone primary free water storage [-]"]
@parameters pfree [description = "Fraction of percolation directed to free water stores [-]"]
@parameters klzp [description = "Primary baseflow runoff coefficient [d-1]"]
@parameters klzs [description = "Supplemental baseflow runoff coefficient [d-1]"]

# Calculate derived parameters
uztwm = f1*smax                                               # Maximum upper zone tension water storage [mm]
uzfwm = max(0.005/4, f2*(smax-uztwm))                        # Maximum upper zone free water storage [mm]
lztwm = max(0.005/4, f3*(smax-uztwm-uzfwm))                 # Maximum lower zone tension water storage [mm]
lzfwpm = max(0.005/4, f4*(smax-uztwm-uzfwm-lztwm))         # Maximum lower zone primary free water storage [mm]
lzfwsm = max(0.005/4, (1-f4)*(smax-uztwm-uzfwm-lztwm))     # Maximum lower zone supplemental free water storage [mm]
pbase = lzfwpm*klzp + lzfwsm*klzs                          # Base percolation rate [mm/d]
zperc = min(100000, (lztwm+lzfwsm*(1-klzs))/(lzfwsm*klzs+lzfwpm*klzp) + (lzfwpm*(1-klzp))/(lzfwsm*klzs+lzfwpm*klzp))  # Base percolation rate multiplication factor [-]

# Bucket 1: Upper zone tension water storage
funcs_1 = [
    SplitFlux((in=prcp), (p1=1-pctim,), Val(1), split=peff),
    SoilMoistureFlux((S=S1, S2=S2), (Smax=uztwm, S2max=uzfwm), Val(1), soilmoisture=ru),
    EvaporationFlux((S=S1, in=pet), (Smax=uztwm,), Val(7), evaporation=euztw),
    SaturationFlux((S=S1, in=peff), (Smax=uztwm,), Val(1), saturation=twexu)
]

dfuncs_1 = [
    StateFlux([peff, ru] => [euztw, twexu], S1)
]

# Bucket 2: Upper zone free water storage
funcs_2 = [
    SaturationFlux((S=S2, in=twexu), (Smax=uzfwm,), Val(1), saturation=qsur),
    InterflowFlux((S=S2,), (p1=kuz,), Val(5), interflow=qint),
    EvaporationFlux((S=S2, in=pet-euztw), NamedTuple(), Val(1), evaporation=euzfw),
    PercolationFlux((S=S2,), (p1=pbase, p2=zperc, p3=rexp, p4=lztwm+lzfwpm+lzfwsm, Smax=uzfwm), Val(4), percolation=pc)
]

dfuncs_2 = [
    StateFlux([twexu] => [qsur, qint, euzfw, ru, pc], S2)
]

# Bucket 3: Lower zone tension water storage
funcs_3 = [
    SplitFlux((in=pc), (p1=1-pfree,), Val(1), split=pctw),
    EvaporationFlux((S=S3, in=pet-euztw-euzfw), (Smax=lztwm,), Val(7), evaporation=elztw),
    SaturationFlux((S=S3, in=pctw), (Smax=lztwm,), Val(1), saturation=twexl),
    SoilMoistureFlux((S=S3, S2=S4), (Smax=lztwm, S2max=lzfwpm, S3=S5, S3max=lzfwsm), Val(2), soilmoisture=rlp),
    SoilMoistureFlux((S=S3, S2=S5), (Smax=lztwm, S2max=lzfwsm, S3=S4, S3max=lzfwpm), Val(2), soilmoisture=rls)
]

dfuncs_3 = [
    StateFlux([pctw, rlp, rls] => [elztw, twexl], S3)
]

# Bucket 4: Lower zone primary free water storage
funcs_4 = [
    BaseflowFlux((S=S4,), (p1=klzp,), Val(1), baseflow=qbfp)
]

dfuncs_4 = [
    StateFlux([twexl*S4/lzfwpm/(S4/lzfwpm+S5/lzfwsm), pc*pfree*S4/lzfwpm/(S4/lzfwpm+S5/lzfwsm)] => [rlp, qbfp], S4)
]

# Bucket 5: Lower zone supplemental free water storage
funcs_5 = [
    BaseflowFlux((S=S5,), (p1=klzs,), Val(1), baseflow=qbfs)
]

dfuncs_5 = [
    StateFlux([twexl*S5/lzfwsm/(S4/lzfwpm+S5/lzfwsm), pc*pfree*S5/lzfwsm/(S4/lzfwpm+S5/lzfwsm)] => [rls, qbfs], S5)
]

# Direct runoff calculation
funcs_d = [
    SplitFlux((in=prcp), (p1=pctim,), Val(1), split=qdir)
]

# Total flow calculation
q_flux = HydroFlux([qdir, qsur, qint, qbfp, qbfs] => [q], exprs=[qdir + qsur + qint + qbfp + qbfs])

# Create buckets and model
sacramento_bucket_1 = HydroBucket(name=:sacramento_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
sacramento_bucket_2 = HydroBucket(name=:sacramento_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
sacramento_bucket_3 = HydroBucket(name=:sacramento_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
sacramento_bucket_4 = HydroBucket(name=:sacramento_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)
sacramento_bucket_5 = HydroBucket(name=:sacramento_bucket_5, funcs=funcs_5, dfuncs=dfuncs_5)

sacramento_model = HydroModel(name=:sacramento, components=[sacramento_bucket_1, sacramento_bucket_2, sacramento_bucket_3, 
                                                          sacramento_bucket_4, sacramento_bucket_5, funcs_d, q_flux])

export sacramento_bucket_1, sacramento_bucket_2, sacramento_bucket_3, sacramento_bucket_4, sacramento_bucket_5, sacramento_model
end
