@reexport module Echo
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4 S5 S6 S3old
@variables prcp pet temp
@variables ei pn ps pr ms fs gs mw ew eq fi rh eps et rd l ls lf rf rs

# Basic parameters
@parameters rho [description = "Maximum interception storage [mm]"]
@parameters ts [description = "Threshold temperature for snowfall [oC]"]
@parameters tm [description = "Threshold temperature for snowmelt [oC]"]
@parameters as [description = "Degree-day factor [mm/oC/d]"]
@parameters af [description = "Refreezing reduction factor [-]"]
@parameters gmax [description = "Maximum melt due to ground-heat flux [mm/d]"]
@parameters the [description = "Water-holding capacity of snow [-]"]
@parameters phi [description = "Maximum infiltration rate [mm/d]"]
@parameters smax [description = "Maximum soil moisture storage [mm]"]
@parameters fsm [description = "Plant stress point as a fraction of Smax [-]"]
@parameters fsw [description = "Wilting point as fraction of sm [-]"]
@parameters ksat [description = "Runoff rate from soil moisture [d-1]"]
@parameters c [description = "Runoff non-linearity from soil moisture [-]"]
@parameters lmax [description = "Groundwater flux [mm/d]"]
@parameters kf [description = "Runoff coefficient [d-1]"]
@parameters ks [description = "Runoff coefficient [d-1]"]

# Auxiliary parameters
sm = fsm * smax  # Plant stress point [mm]
sw = fsw * sm    # Wilting point [mm]

# Bucket 1: Interception storage
funcs_1 = [
    EvaporationFlux((S=S1, in=pet), NamedTuple(), Val(1), evaporation=ei),
    InterceptionFlux((in=prcp, S=S1), (Smax=rho,), Val(1), interception=pn)
]

dfuncs_1 = [
    StateFlux([prcp] => [ei, pn], S1)
]

# Bucket 2: Snow storage
funcs_2 = [
    SnowfallFlux((in=pn, T=temp), (tt=ts,), Val(1), snowfall=ps),
    MeltFlux((S=S2, T=temp), (ddf=as, tt=tm), Val(1), melt=ms),
    RefreezeFlux((S=S3, T=temp), (p1=af, p2=as, tt=tm), Val(1), refreeze=fs),
    MeltFlux((S=S2,), (ddf=gmax,), Val(2), melt=gs)
]

dfuncs_2 = [
    StateFlux([ps, fs] => [ms, gs], S2)
]

# Bucket 3: Snow liquid storage
funcs_3 = [
    RainfallFlux((in=pn, T=temp), (tt=ts,), Val(1), rainfall=pr),
    SaturationFlux((in=pr+ms, S=S3), (Smax=the*S2,), Val(1), saturation=mw),
    ExcessFlux((S=S3old,), (Smax=the*S2,), Val(1), excess=ew),
    EffectiveFlux((in=mw+gs+ew,), NamedTuple(), Val(1), effective=eq),
    InfiltrationFlux((in=eq,), (Smax=phi,), Val(4), infiltration=fi),
    EffectiveFlux((in=eq,), (p1=fi,), Val(1), effective=rh)
]

dfuncs_3 = [
    StateFlux([pr, ms] => [fs, mw, ew], S3)
]

# Bucket 4: Soil moisture storage
funcs_4 = [
    EffectiveFlux((in=pet,), (p1=ei,), Val(1), effective=eps),
    EvaporationFlux((S=S4, in=eps), (p1=sw, p2=sm), Val(22), evaporation=et),
    SaturationFlux((in=fi, S=S4), (Smax=smax,), Val(1), saturation=rd),
    RechargeFlux((S=S4,), (p1=ksat, p2=c), Val(6), recharge=l)
]

dfuncs_4 = [
    StateFlux([fi] => [et, rd, l], S4)
]

# Bucket 5: Fast groundwater storage
funcs_5 = [
    RechargeFlux((in=l,), (p1=lmax,), Val(7), recharge=ls),
    EffectiveFlux((in=l,), (p1=ls,), Val(1), effective=lf),
    BaseflowFlux((S=S5,), (p1=kf,), Val(1), baseflow=rf)
]

dfuncs_5 = [
    StateFlux([lf] => [rf], S5)
]

# Bucket 6: Slow groundwater storage
funcs_6 = [
    BaseflowFlux((S=S6,), (p1=ks,), Val(1), baseflow=rs)
]

dfuncs_6 = [
    StateFlux([ls] => [rs], S6)
]

# Total flow calculation
q_flux = HydroFlux([rh, rf, rs] => [q], exprs=[rh + rf + rs])

# Create buckets and model
echo_bucket_1 = HydroBucket(name=:echo_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
echo_bucket_2 = HydroBucket(name=:echo_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
echo_bucket_3 = HydroBucket(name=:echo_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
echo_bucket_4 = HydroBucket(name=:echo_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)
echo_bucket_5 = HydroBucket(name=:echo_bucket_5, funcs=funcs_5, dfuncs=dfuncs_5)
echo_bucket_6 = HydroBucket(name=:echo_bucket_6, funcs=funcs_6, dfuncs=dfuncs_6)

echo_model = HydroModel(name=:echo, components=[echo_bucket_1, echo_bucket_2, echo_bucket_3,
                                               echo_bucket_4, echo_bucket_5, echo_bucket_6, q_flux])

export echo_bucket_1, echo_bucket_2, echo_bucket_3
export echo_bucket_4, echo_bucket_5, echo_bucket_6
export echo_model
end
