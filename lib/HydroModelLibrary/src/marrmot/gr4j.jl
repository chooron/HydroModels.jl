@reexport module GR4J
using ..HydroModels
using ..HydroModelLibrary: EvaporationFlux, SaturationFlux, PercolationFlux, RechargeFlux, BaseflowFlux

@variables S1 S2 pet prcp pn en ef ps es perc q9 q9_routed q1 q1_routed fr qr qt ex
@parameters x1 x2 x3 x4

funcs_1 = [
    HydroFlux([prcp, pet] => [pn, en], exprs=[max(0.0, prcp - pet), max(0.0, pet - prcp)]),
    HydroFlux([prcp, pn] => [ef], exprs=[prcp - pn]),
    SaturationFlux((S=S1, in=pn), (p1=x1,), Val(4), saturation=ps),
    EvaporationFlux((S=S1, in=en), (p1=x1,), Val(11), evaporation=es),
    PercolationFlux((S=S1,), (p1=x1,), Val(3), percolation=perc),
]

dfuncs_1 = [StateFlux([ps] => [es, perc], S1)]

uh_q9 = UnitHydrograph([q9] => [q9_routed], [x4], uhfunc=UHFunction(:UH_1_HALF), solvetype=:SPARSE)
uh_q1 = UnitHydrograph([q1] => [q1_routed], [x4], uhfunc=UHFunction(:UH_2_FULL), solvetype=:SPARSE)

funcs_2 = [
    RechargeFlux((S=S2,), (Smax=x3, p1=x2), Val(2), recharge=fr),
    BaseflowFlux((S=S2,), (Smax=x3,), Val(3), saturation=qr),
    HydroFlux([qr, q1_routed, fr] => [qt], exprs=[qr + max(q1_routed + fq, 0)]),
    HydroFlux([qt, q1_routed] => [ex], exprs=[qt - q1_routed]),
]
dfuncs_2 = [StateFlux([q9_routed, fr] => [qr], S2)]

gr4j_soil_bucket = HydroBucket(name=:gr4j_soil_bucket, funcs=funcs_1, dfuncs=dfuncs_1)
gr4j_zone_bucket = HydroBucket(name=:gr4j_zone_bucket, funcs=funcs_2, dfuncs=dfuncs_2)
gr4j_model = HydroModel(name=:gr4j, components=[gr4j_soil_bucket, gr4j_zone_bucket])

export gr4j_soil_bucket, gr4j_zone_bucket, gr4j_model
end