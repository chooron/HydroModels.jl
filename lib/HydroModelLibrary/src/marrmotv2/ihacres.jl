@reexport module Ihacres
using ..HydroModels
using ..HydroModelLibrary: EvaporationFlux, SaturationFlux, BaseflowFlux, InterflowFlux

@variables S pet prcp ea u uq us xq xs Qt
@parameters Smax lp d p alpha tau_q tau_s tau_d

funcs = [
    EvaporationFlux((S=S, in=pet), (Smax=Smax, p1=lp,), Val(12), evaporation=ea),
    SaturationFlux((S=S, in=prcp), (p1=d, p2=p), Val(5), saturation=u),
    HydroFlux([u] => [uq, us], [alpha], exprs=[alpha * u, (1 - alpha) * u]),
]

# TODO MARRMoT代码中输入输出方向是反的
dfuncs = [StateFlux([prcp] => [ea, u], S)]

ihacres_bucket = HydroBucket(name=:ihacres_bucket, funcs=funcs, dfuncs=dfuncs)

#* route module
uh_1 = UnitHydrograph([uq] => [xq], [tau_q], uhfunc=UHFunction(:UH_5_HALF), solvetype=:SPARSE)
uh_2 = UnitHydrograph([us] => [xs], [tau_s], uhfunc=UHFunction(:UH_5_HALF), solvetype=:SPARSE)
sum_flux = HydroFlux([xq, xs] => [xx], exprs=[xq + xs])
uh_3 = UnitHydrograph([xx] => [xt], [tau_d], uhfunc=UHFunction(:UH_8_DELAY), solvetype=:SPARSE)

ihacres = HydroModel(name=:ihacres, components=[ihacres_bucket, uh_1, uh_2, sum_flux, uh_3])

export ihacres_bucket, ihacres
end