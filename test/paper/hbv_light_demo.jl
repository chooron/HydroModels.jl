using CSV
using Lux
using LuxCore
using Random
using DataFrames
using Symbolics
using ComponentArrays
using OrdinaryDiffEq
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using BenchmarkTools
using StableRNGs
using Optimization
using OptimizationBBO
using CairoMakie
include("../../src/LumpedHydro.jl")

# load data
df = DataFrame(CSV.File("data/exphydro/01013500.csv"));
ts = collect(1:10000)
prcp_vec = df[ts, "prcp(mm/day)"]
temp_vec = df[ts, "tmean(C)"]
dayl_vec = df[ts, "dayl(day)"]
qobs_vec = df[ts, "flow(mm)"]

#! parameters in the HBV-light model
@parameters TT CFMAX CFR CWH FC beta LP PERC k0 k1 k2 UZL

#! hydrological flux in the Exp-Hydro model
@variables prcp(t) = 0.0 [description = "precipitation", unit = "mm"]
@variables temp(t) = 0.0 [description = "precipitation", unit = "°C"]
@variables lday(t) = 0.0 [description = "length of day", unit = "-"]
@variables pet(t) = 0.0 [description = "potential evapotranspiration", unit = "mm"]
@variables rainfall(t) = 0.0 [description = "rain splitted from precipitation", unit = "mm"]
@variables snowfall(t) = 0.0 [description = "snow splitted from precipitation", unit = "mm"]
@variables refreeze(t) = 0.0 [description = "Refreeze of ponding water"]
@variables melt(t) = 0.0 [description = "snow melt", unit = "mm"]
@variables snowinfil(t) = 0.0 [description = " Snowmelt infiltration", unit = "mm"]
@variables snowpack(t) = 0.0 [description = " Snowmelt infiltration", unit = "mm"]
@variables meltwater(t) = 0.0 [description = " Snowmelt infiltration", unit = "mm"]
@variables soilwater(t) = 0.0 [description = " Snowmelt infiltration", unit = "mm"]

@variables soilwetfrac(t) = 0.0
@variables recharge(t) = 0.0
@variables excess(t) = 0.0
@variables evapfrac(t) = 0.0
@variables evap(t) = 0.0

@variables upperzone(t) = 0.0
@variables lowerzone(t) = 0.0
@variables perc(t) = 0.0
@variables q0(t) = 0.0
@variables q1(t) = 0.0
@variables q2(t) = 0.0
@variables flow(t) = 0.0

SimpleFlux = LumpedHydro.SimpleFlux
StdMeanNormFlux = LumpedHydro.StdMeanNormFlux
NeuralFlux = LumpedHydro.NeuralFlux
LagFlux = LumpedHydro.LagFlux
StateFlux = LumpedHydro.StateFlux
HydroElement = LumpedHydro.HydroElement
HydroUnit = LumpedHydro.HydroUnit
step_func = LumpedHydro.step_func

#! define the snow pack reservoir
snow_funcs = [
    SimpleFlux([prcp, temp] => [snowfall, rainfall], [TT],
        flux_exprs=[step_func(TT - temp) * prcp, step_func(temp - TT) * prcp]),
    SimpleFlux([temp, lday] => [pet],
        flux_exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    SimpleFlux([snowpack, temp] => [melt], [TT, CFMAX],
        flux_exprs=[min(snowpack, CFMAX * max(0.0, temp - TT))]),
    SimpleFlux([meltwater, temp] => [refreeze], [TT, CFMAX, CFR],
        flux_exprs=[min(meltwater, CFR * CFMAX * max(0.0, TT - temp))]),
    SimpleFlux([meltwater, snowpack] => [snowinfil], [CWH],
        flux_exprs=[max(0.0, meltwater - CWH * snowpack)]),
]
snow_dfuncs = [StateFlux([snowfall, refreeze] => [melt], snowpack),
    StateFlux([melt] => [refreeze, snowinfil], meltwater)]
snow_ele = HydroElement(:hbv_snow, funcs=snow_funcs, dfuncs=snow_dfuncs)
#! define the soil water reservoir
soil_funcs = [
    SimpleFlux([soilwater] => [soilwetfrac], [FC, beta],
        flux_exprs=[clamp(abs(soilwater / FC)^beta, 0.0, 1.0)]),
    SimpleFlux([rainfall, snowinfil, soilwetfrac] => [recharge],
        flux_exprs=[(rainfall + snowinfil) * soilwetfrac]),
    SimpleFlux([soilwater] => [excess], [FC],
        flux_exprs=[max(0.0, soilwater - FC)]),
    SimpleFlux([soilwater] => [evapfrac], [LP, FC],
        flux_exprs=[soilwater / (LP * FC)]), SimpleFlux([soilwater, pet, evapfrac] => [evap],
        flux_exprs=[min(soilwater, pet * evapfrac)]),
]
soil_dfuncs = [StateFlux([rainfall, snowinfil] => [evap, excess, recharge], soilwater)]
soil_ele = HydroElement(:hbv_soil, funcs=soil_funcs, dfuncs=soil_dfuncs)
#! define the upper and lower subsurface zone 
zone_funcs = [
    SimpleFlux([soilwater, upperzone, evapfrac] => [perc], [PERC],
        flux_exprs=[min(soilwater, PERC)]),
    SimpleFlux([upperzone] => [q0], [k0, UZL],
        flux_exprs=[max(0.0, (upperzone - UZL) * k0)]),
    SimpleFlux([upperzone] => [q1], [k1],
        flux_exprs=[upperzone * k1]),
    SimpleFlux([lowerzone] => [q2], [k2],
        flux_exprs=[lowerzone * k2]),
    SimpleFlux([q0, q1, q2] => [flow],
        flux_exprs=[q0 + q1 + q2]),
]
zone_dfuncs = [StateFlux([recharge, excess] => [perc, q0, q1], upperzone), StateFlux([perc] => [q2], lowerzone)]
zone_ele = HydroElement(:hbv_zone, funcs=zone_funcs, dfuncs=zone_dfuncs)
#! define the HBV-light model
hbv_model = HydroUnit(:hbv, components=[snow_ele, soil_ele, zone_ele]);

params = ComponentVector(TT=0.0, CFMAX=5.0, CWH=0.1, CFR=0.05, FC=200.0, LP=0.6, beta=3.0, k0=0.06, k1=0.2, k2=0.1, PERC=2, UZL=10)
init_states = ComponentVector(upperzone=0.0, lowerzone=0.0, soilwater=0.0, meltwater=0.0, snowpack=0.0)
pas = ComponentVector(params=params, initstates=init_states)
input = (prcp=prcp_vec, lday=dayl_vec, temp=temp_vec)
result = hbv_model(input, pas, timeidx=ts)


fig = Figure(size=(400, 300))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
lines!(ax, ts, result.flow, color=:red)
lines!(ax, ts, qobs_vec, color=:blue)
fig

#! set the tunable parameters boundary
lower_bounds = [-1.5, 1, 0.0, 0.0, 50.0, 0.3, 1.0, 0.05, 0.01, 0.001, 0.0, 0.0]
upper_bounds = [1.2, 8.0, 0.2, 0.1, 500.0, 1.0, 6.0, 0.5, 0.3, 0.15, 3.0, 70.0]
#! prepare flow
output = (flow=qobs_vec,)
#* ComponentVector{Float64}(params = (TT = -1.2223657527438707, CFMAX = 2.201359793941345, CWH = 0.022749518921432663, CFR = 0.058335602629828544, FC = 160.01327559173077, LP = 0.7042581781418978, 
#* beta = 5.580695551758287, k0 = 0.0500023960318018, k1 = 0.04573064980956475, k2 = 0.14881856483902567, PERC = 1.3367222956722589, UZL = 44.059927907190016))
#! model calibration
best_pas = LumpedHydro.param_box_optim(
    hbv_model,
    tunable_pas=ComponentVector(params=params),
    const_pas=ComponentVector(initstates=init_states),
    input=input,
    target=output,
    timeidx=ts,
    lb=lower_bounds,
    ub=upper_bounds,
    solve_alg=BBO_adaptive_de_rand_1_bin_radiuslimited(),
    maxiters=10000,
    loss_func=LumpedHydro.mse,
)

result = hbv_model(input, LumpedHydro.merge_ca(pas, best_pas), timeidx=ts)
LumpedHydro.nse(result.flow, qobs_vec)