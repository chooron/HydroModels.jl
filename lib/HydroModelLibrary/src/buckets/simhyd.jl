@reexport module SIMHYD
# https://github.com/hydrogo/LHMP/blob/master/models/simhyd_cemaneige.py
using ..HydroModels

@variables prcp pet U IMAX INT INR RMO IRUN EVAP SRUN REC SMF POT BAS GWF SMS GW
@parameters INSC COEFF SQ SMSC SUB CRAK K etmul DELAY X_m

funcs = [
    HydroFlux([pet] => [IMAX], [INSC], exprs=[min(INSC, pet)]),
    HydroFlux([prcp, IMAX] => [INT], exprs=[min(IMAX, prcp)]),
    HydroFlux([prcp, INT] => [INR], exprs=[prcp - INT]),
    HydroFlux([pet, INT] => [POT], exprs=[pet - INT]),
    HydroFlux([POT, SMS] => [EVAP], [SMSC], exprs=[min(10 * (SMS / SMSC), POT)]),
    HydroFlux([INR, SMS] => [RMO], [COEFF, SQ, SMSC], exprs=[min(COEFF * exp(-SQ * (SMS / SMSC)), INR)]),
    HydroFlux([INR, RMO] => [IRUN], exprs=[INR - RMO]),
    HydroFlux([RMO, SMS] => [SRUN], [SUB, SMSC], exprs=[SUB * (SMS / SMSC) * RMO]),
    HydroFlux([RMO, SRUN, SMS] => [REC], [CRAK, SMSC], exprs=[CRAK * (SMS / SMSC) * (RMO - SRUN)]),
    HydroFlux([RMO, SRUN, REC] => [SMF], exprs=[RMO - SRUN - REC]),
    HydroFlux([SMS, SMF] => [GWF], [SMSC], exprs=[ifelse(SMS == SMSC, SMF, 0.0)]),
    HydroFlux([GW] => [BAS], [K], exprs=[K * GW]),
    HydroFlux([IRUN, SRUN, BAS] => [U], [K], exprs=[IRUN + SRUN + BAS]),
]

dfuncs = [
    StateFlux([SMF] => [EVAP, GWF], SMS),
    StateFlux([REC, GWF] => [BAS], GW)
]

soil_bucket = HydroBucket(name=:simhyd_soil, funcs=funcs, dfuncs=dfuncs)

simhyd_model = HydroModel(name=:simhyd, components=[soil_bucket])
export soil_bucket, simhyd_model
end