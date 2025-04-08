@reexport module GR4J

using ..HydroModels

#* parameters in the GR4J model
@parameters x1 x2 x3 x4 lag area_coef
#* hydrological flux in the production store
@variables prcp ep soilwater pn en ps es perc pr slowflow fastflow
#* hydrological flux in the unit hydrograph
@variables slowflow_routed fastflow_routed
#* hydrological flux in the routing store
@variables routingstore exch routedflow flow q q_routed s_river

#* define the production store
prod_funcs = [
    HydroFlux([prcp, ep] => [pn, en], exprs=[prcp - min(prcp, ep), ep - min(prcp, ep)]),
    HydroFlux([pn, soilwater] => [ps], [x1], exprs=[max(0.0, pn * (1 - (soilwater / x1)^2))]),
    HydroFlux([en, soilwater] => [es], [x1], exprs=[en * (2 * soilwater / x1 - (soilwater / x1)^2)]),
    HydroFlux([soilwater] => [perc], [x1], exprs=[((x1)^(-4)) / 4 * ((4 / 9)^(4)) * (soilwater^5)]),
    HydroFlux([pn, ps, perc] => [pr], exprs=[pn - ps + perc]),
    HydroFlux([pr] => [slowflow, fastflow], exprs=[0.9 * pr, 0.1 * pr]),
]
prod_dfuncs = [StateFlux([ps] => [es, perc], soilwater)]
prod_ele = HydroBucket(name=:gr4j_prod, funcs=prod_funcs, dfuncs=prod_dfuncs)

#* define the routing store
rst_funcs = [
    HydroFlux([routingstore] => [exch], [x2, x3], exprs=[x2 * abs(routingstore / x3)^3.5]),
    HydroFlux([routingstore, slowflow, exch] => [routedflow], [x3], exprs=[x3^(-4) / 4 * (routingstore + slowflow + exch)^5]),
    HydroFlux([routedflow, fastflow_routed, exch] => [flow], exprs=[routedflow + max(fastflow_routed + exch, 0.0)]),
]
rst_dfuncs = [StateFlux([exch, slowflow] => [routedflow], routingstore)]
rst_ele = HydroBucket(name=:gr4j_rst, funcs=rst_funcs, dfuncs=rst_dfuncs)


uh_1 = UnitHydrograph(slowflow, slowflow_routed, x4, uhfunc=HydroModels.UHFunction(:UH_1_HALF), solvetype=:SPARSE)
uh_2 = UnitHydrograph(fastflow, fastflow_routed, x4, uhfunc=HydroModels.UHFunction(:UH_2_FULL), solvetype=:SPARSE)

gr4j_model = HydroModel(name=:gr4j, components=[prod_ele, uh_1, uh_2, rst_ele])

export gr4j_model, prod_ele, rst_ele, uh_1, uh_2
end