include("../src/HydroModels.jl")
using ComponentArrays
using ModelingToolkit
using CSV
using DataFrames

@variables temp lday pet pe
@parameters p1 p2 p3

file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])

pet_flux = HydroModels.HydroFlux(
    [temp, lday] => [pet, pe], [p1, p2, p3],
    exprs=[p1 * lday * 24 * p2 * exp((p3 * temp) / (temp + 237.3)) / (temp + 273.2),
           p1 * lday * 24 * p2 * exp((p3 * temp) / (temp + 237.3)) / (temp + 273.2)],
)

params = (p1=29.8, p2=0.611, p3=17.3)
pas = ComponentVector(params=params)
input_arr = reduce(hcat, collect(input[HydroModels.get_input_names(pet_flux)])) |> permutedims
result1 = pet_flux(input_arr, pas)

# pas = (;params=ComponentVector(p1=ones(10) .* 29.8, p2=ones(10) .* 0.611, p3=ones(10) .* 17.3))
# input_arr3 = repeat(reshape(input_arr, size(input_arr, 1), 1, size(input_arr, 2)), 1, 10, 1)
# result2 = pet_flux(input_arr3, pas, config=(; ptyidx=1:10))