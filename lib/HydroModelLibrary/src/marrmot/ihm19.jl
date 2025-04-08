@reexport module Ihm19
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4
@variables prcp pet temp
@variables ei pex pexmp pexs1 fmp qexmp qmp pqexs1 fs1 etas1 qs1 q0 q0r qmps1 pc qh qs2 qgw

# Basic parameters
@parameters SIMAX [description = "Maximum interception storage [mm]"]
@parameters A [description = "Splitting coeffcient for excess precipitation [-]"]
@parameters FF [description = "Forest fraction [-]"]
@parameters SMPMAX [description = "Maximum storage macropores [mm]"]
@parameters CQMP [description = "Runoff time parameter (fast/slow runnoff) first soil layer [1/d]"]
@parameters XQMP [description = "Runoff scale parameter first soil layer [-]"]
@parameters SS1MAX [description = "Maximum soil moisture storage first soil layer [mm]"]
@parameters FCCS1 [description = "Field capacity coefficient fist soil layer [-]"]
@parameters CFS1 [description = "Maximum infiltration rate first soil layer [-]"]
@parameters XFS1 [description = "Infiltration loss exponent first soil layer [-]"]
@parameters CQS1 [description = "Runoff time parameter for (fast/slow runnoff) first soil layer [1/d]"]
@parameters XQS1 [description = "Runoff scale parameter first soil layer [-]"]
@parameters SS2MAX [description = "Maximum soil moisture storage second soil layer [mm]"]
@parameters CQS2 [description = "Runoff time parameter for (fast/slow runnoff) second soil layer [1/d]"]
@parameters XQS2 [description = "Runoff scale parameter second soil layer [-]"]
@parameters D [description = "Flow delay before surface runoff [d]"]

# Initial conditions
S0 = Dict(
    :S1 => 0.8 * SIMAX,   # Initial interception storage
    :S2 => 0.8 * SMPMAX,  # Initial macropore storage
    :S3 => 0.8 * SS1MAX,  # Initial soil moisture storage 1
    :S4 => 0.8 * SS2MAX   # Initial soil moisture storage 2
)

# Bucket 1: Interception storage
funcs_1 = [
    EvaporationFlux((S=S1, in=pet), NamedTuple(), Val(1), evaporation=ei),
    InterceptionFlux((in=prcp, S=S1), (Smax=SIMAX,), Val(1), interception=pex)
]

dfuncs_1 = [
    StateFlux([prcp] => [ei, pex], S1)
]

# Bucket 2: Macropore storage
funcs_2 = [
    SplitFlux((in=pex,), (p1=A,), Val(1), split=pexmp),
    InfiltrationFlux((in=pexmp, S=S2), (Smax=SMPMAX,), Val(3), infiltration=fmp),
    EffectiveFlux((in=pexmp,), (p1=fmp,), Val(1), effective=qexmp),
    InterflowFlux((S=S2,), (p1=CQMP, p2=XQMP), Val(3), interflow=qmp)
]

dfuncs_2 = [
    StateFlux([fmp] => [qmp], S2)
]

# Bucket 3: First soil layer storage
funcs_3 = [
    SplitFlux((in=pex,), (p1=A,), Val(2), split=pexs1),
    EffectiveFlux((in=pexs1+qexmp,), NamedTuple(), Val(1), effective=pqexs1),
    InfiltrationFlux((S=S3, in=pqexs1), (p1=CFS1, p2=XFS1, Smax=SS1MAX), Val(7), infiltration=fs1),
    EvaporationFlux((S=S3, in=pet), (p1=FF, p2=FCCS1, Smax=SS1MAX), Val(23), evaporation=etas1),
    InterflowFlux((S=S3,), (p1=CQS1, p2=FCCS1*SS1MAX, p3=XQS1), Val(9), interflow=qs1),
    EffectiveFlux((in=pqexs1,), (p1=fs1,), Val(1), effective=q0)
]

dfuncs_3 = [
    StateFlux([fs1] => [etas1, qs1], S3)
]

# Bucket 4: Second soil layer storage
funcs_4 = [
    EffectiveFlux((in=qmp+qs1,), NamedTuple(), Val(1), effective=qmps1),
    InfiltrationFlux((in=qmps1, S=S4), (Smax=SS2MAX,), Val(3), infiltration=pc),
    EffectiveFlux((in=qmps1,), (p1=pc,), Val(1), effective=qh),
    InterflowFlux((S=S4,), (p1=CQS2, p2=XQS2), Val(3), interflow=qs2)
]

dfuncs_4 = [
    StateFlux([pc] => [qs2], S4)
]

# Unit hydrograph routing
routing = UnitHydrograph(name=:routing, delay=D, input=q0, output=q0r)

# Total flow calculation
q_flux = HydroFlux([q0r, qh, qs2] => [q], exprs=[q0r + qh + qs2 + 0.0195])  # Adding constant groundwater flow

# Create buckets and model
ihm19_bucket_1 = HydroBucket(name=:ihm19_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
ihm19_bucket_2 = HydroBucket(name=:ihm19_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
ihm19_bucket_3 = HydroBucket(name=:ihm19_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
ihm19_bucket_4 = HydroBucket(name=:ihm19_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)

ihm19_model = HydroModel(name=:ihm19, components=[ihm19_bucket_1, ihm19_bucket_2,
                                                 ihm19_bucket_3, ihm19_bucket_4,
                                                 routing, q_flux],
                        initial_state=S0)

export ihm19_bucket_1, ihm19_bucket_2, ihm19_bucket_3, ihm19_bucket_4
export ihm19_model
end
