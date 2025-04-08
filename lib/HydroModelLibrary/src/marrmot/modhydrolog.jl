@reexport module Modhydrolog
using ..HydroModels
using ..HydroModelLibrary

@variables S1 S2 S3 S4 S5
@variables prcp pet temp
@variables Ei EXC INF INT REC SMF Et GWF RUN TRAP Ed DINF SEEP SRUN FLOW Q

# Basic parameters
@parameters insc [description = "Maximum interception capacity [mm]"]
@parameters coeff [description = "Maximum infiltration loss parameter [-]"]
@parameters sq [description = "Infiltration loss exponent [-]"]
@parameters smsc [description = "Maximum soil moisture capacity [mm]"]
@parameters sub [description = "Proportionality constant [-]"]
@parameters crak [description = "Proportionality constant [-]"]
@parameters em [description = "Plant-controled maximum evaporation rate [mm/d]"]
@parameters dsc [description = "Maximum depression capacity [mm]"]
@parameters ads [description = "Land fraction functioning as depression storage [-]"]
@parameters md [description = "Depression storage parameter [-]"]
@parameters vcond [description = "Leakage coefficient [mm/d]"]
@parameters dlev [description = "Datum around which groundwater fluctuates relative to river bed [mm]"]
@parameters k1 [description = "Flow exchange parameter [d-1]"]
@parameters k2 [description = "Flow exchange parameter [d-1]"]
@parameters k3 [description = "Flow exchange parameter [d-1]"]

# Bucket 1: Interception storage
funcs_1 = [
    EvaporationFlux((S=S1, in=pet), NamedTuple(), Val(1), evaporation=Ei),
    InterceptionFlux((S=S1, in=prcp), (Smax=insc,), Val(1), interception=EXC)
]

dfuncs_1 = [
    StateFlux([prcp] => [Ei, EXC], S1)
]

# Bucket 2: Soil moisture storage
funcs_2 = [
    InfiltrationFlux((S=S2, in=EXC), (p1=coeff, p2=sq, Smax=smsc), Val(1), infiltration=INF),
    InterflowFlux((S=S2, in=INF), (p1=sub, Smax=smsc), Val(1), interflow=INT),
    RechargeFlux((S=S2, in=INF-INT), (p1=crak, Smax=smsc), Val(1), recharge=REC),
    SoilMoistureFlux((S=S2,), NamedTuple(), Val(1), soilmoisture=SMF, exprs=[INF - INT - REC]),
    EvaporationFlux((S=S2, in=pet), (Smax=smsc, p1=em), Val(2), evaporation=Et),
    SaturationFlux((S=S2, in=SMF), (Smax=smsc,), Val(1), saturation=GWF),
    RunoffFlux((in=EXC,), (p1=INF,), Val(1), runoff=RUN)
]

dfuncs_2 = [
    StateFlux([SMF, DINF] => [Et, GWF], S2)
]

# Bucket 3: Depression storage
funcs_3 = [
    DepressionFlux((S=S3, in=RUN), (p1=ads, p2=md, Smax=dsc), Val(1), depression=TRAP),
    EvaporationFlux((S=S3, in=ads*pet), NamedTuple(), Val(1), evaporation=Ed),
    InfiltrationFlux((S=S2, in=SMF, S2=S3), (p1=coeff, p2=sq, Smax=smsc), Val(2), infiltration=DINF)
]

dfuncs_3 = [
    StateFlux([TRAP] => [Ed, DINF], S3)
]

# Bucket 4: Groundwater storage
funcs_4 = [
    ExchangeFlux((S=S4,), (p1=vcond, p2=dlev), Val(3), exchange=SEEP),
    RunoffFlux((in=RUN,), (p1=TRAP,), Val(1), runoff=SRUN),
    ExchangeFlux((S=S4, in=SRUN), (p1=k1, p2=k2, p3=k3), Val(1), exchange=FLOW)
]

dfuncs_4 = [
    StateFlux([REC, GWF] => [SEEP, FLOW], S4)
]

# Bucket 5: Channel storage
funcs_5 = [
    BaseflowFlux((S=S5,), (p1=1,), Val(1), baseflow=Q)
]

dfuncs_5 = [
    StateFlux([SRUN, INT, FLOW] => [Q], S5)
]

# Create buckets and model
modhydrolog_bucket_1 = HydroBucket(name=:modhydrolog_bucket_1, funcs=funcs_1, dfuncs=dfuncs_1)
modhydrolog_bucket_2 = HydroBucket(name=:modhydrolog_bucket_2, funcs=funcs_2, dfuncs=dfuncs_2)
modhydrolog_bucket_3 = HydroBucket(name=:modhydrolog_bucket_3, funcs=funcs_3, dfuncs=dfuncs_3)
modhydrolog_bucket_4 = HydroBucket(name=:modhydrolog_bucket_4, funcs=funcs_4, dfuncs=dfuncs_4)
modhydrolog_bucket_5 = HydroBucket(name=:modhydrolog_bucket_5, funcs=funcs_5, dfuncs=dfuncs_5)

modhydrolog_model = HydroModel(name=:modhydrolog, components=[modhydrolog_bucket_1, modhydrolog_bucket_2, 
                                                             modhydrolog_bucket_3, modhydrolog_bucket_4, 
                                                             modhydrolog_bucket_5])

export modhydrolog_bucket_1, modhydrolog_bucket_2, modhydrolog_bucket_3, modhydrolog_bucket_4, modhydrolog_bucket_5, modhydrolog_model
end
