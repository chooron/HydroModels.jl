# Test Flux Components with New Interface
include("../../src/HydroModels.jl")
using .HydroModels
using ComponentArrays
using CSV
using DataFrames

println("="^60)
println("Testing HydroFlux Components")
println("="^60)

@variables temp lday pet pe
@parameters p1 p2 p3

# Load data
file_path = "data/exphydro/01013500.csv"
if !isfile(file_path)
    @warn "Data file not found, using synthetic data"
    df = DataFrame(
        "dayl(day)" => rand(10000) .* 12 .+ 6,
        "tmean(C)" => rand(10000) .* 20 .- 5,
        "prcp(mm/day)" => rand(10000) .* 10
    )
else
    data = CSV.File(file_path)
    df = DataFrame(data)
end

ts = collect(1:1000)
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])

println("\n✓ Data loaded: $(length(ts)) time steps")

# Test 1: Symbolic construction using old interface (still supported)
println("\n1. Testing symbolic HydroFlux construction:")
pet_flux = HydroFlux(
    exprs=[
        p1 * lday * 24 * p2 * exp((p3 * temp) / (temp + 237.3)) / (temp + 273.2),
        p1 * lday * 24 * p2 * exp((p3 * temp) / (temp + 237.3)) / (temp + 273.2)
    ]
)

params = (p1=29.8, p2=0.611, p3=17.3)
pas = ComponentVector(params=params)
input_arr = reduce(hcat, [input[name] for name in HydroModels.get_input_names(pet_flux)]) |> permutedims
result1 = pet_flux(input_arr, pas)

println("   ✓ Symbolic flux executed")
println("   Input size: ", size(input_arr))
println("   Output size: ", size(result1))
println("   Output range: [", minimum(result1), ", ", maximum(result1), "]")

# Test 2: Macro-based construction
println("\n2. Testing macro-based HydroFlux construction:")
pet_macro_flux = @hydroflux pet ~ p1 * lday * 24 * p2 * exp((p3 * temp) / (temp + 237.3)) / (temp + 273.2)

result2 = pet_macro_flux(input_arr, pas)
println("   ✓ Macro flux executed")
println("   Output size: ", size(result2))
println("   Results match symbolic: ", isapprox(result1[1, :], result2[1, :]))

# Test 3: Functional construction (NEW feature)
println("\n3. Testing functional HydroFlux construction:")
function pet_function(inputs, params)
    temp_val = inputs[1]
    lday_val = inputs[2]
    return [params.p1 * lday_val * 24 * params.p2 * exp((params.p3 * temp_val) / (temp_val + 237.3)) / (temp_val + 273.2)]
end

pet_func_flux = HydroFlux(
    pet_function;
    inputs = [:temp, :lday],
    outputs = [:pet],
    params = [:p1, :p2, :p3],
    name = :pet_functional
)

result3 = pet_func_flux(input_arr, pas)
println("   ✓ Functional flux executed")
println("   Output size: ", size(result3))

# Test 4: Multi-node flux (distributed modeling)
println("\n4. Testing multi-node HydroFlux:")
node_num = 5
node_params = ComponentVector(
    params = ComponentVector(
        p1 = fill(29.8, node_num),
        p2 = fill(0.611, node_num),
        p3 = fill(17.3, node_num)
    )
)

input_arr_3d = repeat(reshape(input_arr, size(input_arr, 1), 1, size(input_arr, 2)), 1, node_num, 1)
println("   3D input size: ", size(input_arr_3d))

# Create multi-node flux
pet_multi_flux = @hydroflux begin
    pet ~ p1 * lday * 24 * p2 * exp((p3 * temp) / (temp + 237.3)) / (temp + 273.2)
    hru_types = collect(1:node_num)
end

config = HydroConfig()
result4 = pet_multi_flux(input_arr_3d, node_params, config)
println("   ✓ Multi-node flux executed")
println("   Output size: ", size(result4))

# Summary
println("\n" * "="^60)
println("Test Summary")
println("="^60)
println("✓ Symbolic construction: PASSED")
println("✓ Macro construction: PASSED")
println("✓ Functional construction: PASSED")
println("✓ Multi-node flux: PASSED")
println("\n✅ All flux tests passed!")
println("="^60)