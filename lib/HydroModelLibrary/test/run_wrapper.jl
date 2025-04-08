@parameters Tmin Tmax Df Smax f Qmax
@variables prcp temp lday pet snowpack soilwater rainfall snowfall evap melt baseflow surfaceflow flow

#! load data
ts = collect(1:100)
df = DataFrame(CSV.File("data/exphydro/01013500.csv"))
input_ntp = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
input = Matrix(reduce(hcat, collect(input_ntp[[:temp, :lday, :prcp]]))')

initstates = ComponentVector(snowpack=0.0, soilwater=1303.00)
params = ComponentVector(f=0.0167, Smax=1709.46, Qmax=18.47, Df=2.674, Tmax=0.17, Tmin=-2.09)
pas = ComponentVector(initstates=initstates, params=params)

#! define the snow pack reservoir
snow_fluxes = [
    HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
]
snow_dfluxes = [StateFlux([snowfall] => [melt], snowpack)]
snow_ele = HydroBucket(fluxes=snow_fluxes, dfluxes=snow_dfluxes)

#! define the soil water reservoir
soil_fluxes = [
    HydroFlux([soilwater, pet] => [evap], [Smax], exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
    HydroFlux([soilwater] => [baseflow], [Smax, Qmax, f], exprs=[step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))]),
    HydroFlux([soilwater] => [surfaceflow], [Smax], exprs=[max(0.0, soilwater - Smax)]),
    HydroFlux([baseflow, surfaceflow] => [flow], exprs=[baseflow + surfaceflow]),
]
soil_dfluxes = [StateFlux([rainfall, melt] => [evap, flow], soilwater)]
soil_ele = HydroBucket(fluxes=soil_fluxes, dfluxes=soil_dfluxes)

#! define the Exp-Hydro model
model = HydroModel(name=:exphydro, components=[snow_ele, soil_ele])

@testset "RecordComponentState wrapper" begin
    wrapper = RecordComponentState(model, initstates)
    for i in 1:10
        println(wrapper.states)
        output1 = wrapper(input, pas)
    end
end

@testset "EstimateComponentParams wrapper" begin
    @parameters S0 Q0 ks kq

    est_func1 = HydroFlux([S0] => [Smax, Df], [ks], exprs=[step_func(S0) * S0 * ks, ks * 2.674])
    est_func2 = HydroFlux([Q0] => [Qmax], [kq], exprs=[step_func(Q0) * Q0 * kq])
    est_func1.func
    wrapper1 = EstimateComponentParams(model, [est_func1, est_func2])
    # @test Set(get_param_names(wrapper1.component)) == Set([:Tmin, :Tmax, :Df, :f, :S0, :Q0, :ks, :kq])
    new_params = ComponentVector(f=0.0167, Df=0.0, Tmax=0.17, Tmin=-2.09, S0=1709.46, Q0=18.47, ks=1.0, kq=1.0, Smax=0.0, Qmax=0.0)
    new_pas = ComponentVector(params=new_params, initstates=initstates)
    output = wrapper1(input, new_pas)
    
    pkeys = [:ptype1, :ptype2, :ptype3]
    ptypes = [:ptype1, :ptype2, :ptype3, :ptype1, :ptype2, :ptype3, :ptype1, :ptype2, :ptype3, :ptype1]
    stypes = [Symbol("hru_$i") for i in 1:10]
    wrapper2 = EstimateComponentParams(model, [est_func1, est_func2], pkeys)
    multi_input = fill(input_ntp, 10)
    multi_pas = ComponentVector(params=NamedTuple{Tuple(pkeys)}(fill(new_params, 3)), initstates=NamedTuple{Tuple(stypes)}(fill(initstates, 10)))
    cfg = (ptypes=ptypes, stypes=stypes)
    output = wrapper2(multi_input, multi_pas, config=cfg)
end
