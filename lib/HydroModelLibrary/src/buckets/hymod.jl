@reexport module HyMOD
# https://github.com/KMarkert/hymod
using ..HydroModels

@variables soilwater prcp pet ct_prev raineff1 dummy tmp_soilwater raineff2 evap new_soilwater raineff
@parameters cmax bexp

soil_funcs = [
    HydroFlux([soilwater, prcp, pet] => [ct_prev], [cmax, bexp], exprs=[cmax * (1 - abs((1 - ((bexp + 1) * (soilwater) / cmax)))^(1 / (bexp + 1)))]),
    HydroFlux([prcp, ct_prev] => [raineff1], [cmax], exprs=[max((prcp - cmax + ct_prev), 0.0)]),
    HydroFlux([prcp, ct_prev, raineff1] => [dummy], [cmax], exprs=[min(((ct_prev + prcp - raineff1) / cmax), 1.0)]),
    HydroFlux([dummy] => [tmp_soilwater], [cmax, bexp], exprs=[(cmax / (bexp + 1)) * (1 - abs(1 - dummy)^(bexp + 1))]),
    HydroFlux([soilwater, tmp_soilwater, prcp, raineff1] => [raineff2], exprs=[max(prcp - raineff1 - (tmp_soilwater - soilwater), 0.0)]),
    HydroFlux([tmp_soilwater, pet] => [evap], [cmax, bexp], exprs=[(1 - (((cmax / (bexp + 1)) - tmp_soilwater) / (cmax / (bexp + 1)))) * pet]),
    HydroFlux([tmp_soilwater, evap] => [new_soilwater], exprs=[max(tmp_soilwater - evap, 0.0)]),
    HydroFlux([raineff1, raineff2] => [raineff], exprs=[raineff1 + raineff2])
]

soil_dfuncs = [
    StateFlux(new_soilwater => soilwater)
]

soil_bucket = HydroBucket(name=Symbol(:hymod_soil), funcs=soil_funcs, dfuncs=soil_dfuncs)

@variables raineff slowwater new_slowwater slow_q0 slow_q1 fastwater1 fastwater2 fastwater3
@variables new_fastwater1 new_fastwater2 new_fastwater3 fast_q0 fast_q1 fast_q2 fast_q3 flow

@parameters alpha ks kf

zone_funcs = [
    HydroFlux([raineff] => [fast_q0, slow_q0], [alpha], exprs=[alpha * raineff, (1 - alpha) * raineff]),
    # slow reservoir route
    HydroFlux([slowwater, slow_q0] => [new_slowwater], [ks], exprs=[(1 - ks) * (slowwater + slow_q0)]),
    HydroFlux([new_slowwater] => [slow_q1], [ks], exprs=[(ks / (1 - ks)) * new_slowwater]),
    # fast reservoir route
    HydroFlux([fastwater1, fast_q0] => [new_fastwater1], [kf], exprs=[(1 - kf) * (fastwater1 + fast_q0)]),
    HydroFlux([new_fastwater1] => [fast_q1], [kf], exprs=[(kf / (1 - kf)) * new_fastwater1]),
    HydroFlux([fastwater2, fast_q1] => [new_fastwater2], [kf], exprs=[(1 - kf) * (fastwater2 + fast_q1)]),
    HydroFlux([new_fastwater2] => [fast_q2], [kf], exprs=[(kf / (1 - kf)) * new_fastwater2]),
    HydroFlux([fastwater3, fast_q2] => [new_fastwater3], [kf], exprs=[(1 - kf) * (fastwater3 + fast_q2)]),
    HydroFlux([new_fastwater3] => [fast_q3], [kf], exprs=[(kf / (1 - kf)) * new_fastwater3]),
    # get final output
    HydroFlux([slow_q1, fast_q3] => [flow], exprs=[max(0.0, slow_q1 + fast_q3)])
]

zone_dfuncs = [
    StateFlux(new_slowwater => slowwater),
    StateFlux(new_fastwater1 => fastwater1),
    StateFlux(new_fastwater2 => fastwater2),
    StateFlux(new_fastwater3 => fastwater3),
]

zone_bucket = HydroBucket(name=:hymod_zone, funcs=zone_funcs, dfuncs=zone_dfuncs)

hymod_model = HydroModel(name=:hymod, components=[soil_bucket, zone_bucket])

export soil_bucket, zone_bucket, hymod_model
end
