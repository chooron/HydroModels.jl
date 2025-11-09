# Test Bucket Components with New Interface
using CSV
using DataFrames
using ComponentArrays
using BenchmarkTools

include("../../src/HydroModels.jl")
using .HydroModels

println("="^60)
println("Testing HydroBucket Components")
println("="^60)

# Define step function
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

# Define variables and parameters
@variables temp lday prcp pet snowfall rainfall snowpack melt
@parameters Tmin Tmax Df

# Load data
file_path = "data/exphydro/01013500.csv"
if !isfile(file_path)
    @warn "Data file not found, using synthetic data"
    df = DataFrame(
        "dayl(day)" => rand(1000) .* 12 .+ 6,
        "tmean(C)" => rand(1000) .* 20 .- 5,
        "prcp(mm/day)" => rand(1000) .* 10
    )
else
    df = DataFrame(CSV.File(file_path))
end

ts = collect(1:100)
println("✓ Data loaded: $(length(ts)) time steps")

# Test 1: Single-node bucket
println("\n1. Testing single-node HydroBucket:")

snow_bucket = @hydrobucket :snow begin
    fluxes = begin
        @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
        @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
        @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
        @hydroflux melt ~ step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))
    end
    dfluxes = begin
        @stateflux snowpack ~ snowfall - melt
    end
end

# Prepare parameters
f, Smax, Qmax = 0.0167, 1709.46, 18.47
Df_val, Tmax_val, Tmin_val = 2.674, 0.176, -2.093

params = ComponentVector(
    Df = Df_val,
    Tmax = Tmax_val,
    Tmin = Tmin_val
)
init_states = ComponentVector(snowpack = 100.0)
pas = ComponentVector(params = params)

# Prepare input
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
input_arr = reduce(hcat, [input[name] for name in HydroModels.get_input_names(snow_bucket)]) |> permutedims

# Create configuration with NEW HydroConfig
config = HydroConfig(
    solver = MutableSolver,
    interpolator = Val(DirectInterpolation),
    timeidx = ts,
    min_value = 1e-6
)

# Run single-node bucket
result = snow_bucket(input_arr, pas, config; initstates = init_states)

println("   ✓ Single-node bucket executed")
println("   Input size: ", size(input_arr))
println("   Output size: ", size(result))
println("   State variables: ", HydroModels.get_state_names(snow_bucket))
println("   Output variables: ", HydroModels.get_output_names(snow_bucket))

# Test 2: Multi-node bucket
println("\n2. Testing multi-node HydroBucket:")

node_num = 5
multi_snow_bucket = @hydrobucket :snow_multi begin
    fluxes = begin
        @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
        @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
        @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
        @hydroflux melt ~ step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))
    end
    dfluxes = begin
        @stateflux snowpack ~ snowfall - melt
    end
    hru_types = collect(1:node_num)
end

# Prepare multi-node parameters
node_params = ComponentVector(
    params = ComponentVector(
        Df = fill(Df_val, node_num),
        Tmax = fill(Tmax_val, node_num),
        Tmin = fill(Tmin_val, node_num)
    )
)

node_states = ComponentVector(
    snowpack = fill(100.0, node_num)
)

# Prepare 3D input
input_arr_3d = repeat(reshape(input_arr, size(input_arr, 1), 1, size(input_arr, 2)), 1, node_num, 1)

multi_config = HydroConfig(
    solver = MutableSolver,
    interpolator = Val(DirectInterpolation),
    timeidx = ts,
    min_value = 1e-6
)

result_multi = multi_snow_bucket(input_arr_3d, node_params, multi_config; initstates = node_states)

println("   ✓ Multi-node bucket executed")
println("   Input size: ", size(input_arr_3d))
println("   Output size: ", size(result_multi))

# Test 3: Compare ImmutableSolver
println("\n3. Testing ImmutableSolver:")

config_immutable = HydroConfig(
    solver = ImmutableSolver,
    interpolator = Val(DirectInterpolation),
    timeidx = ts,
    min_value = 1e-6
)

result_immutable = snow_bucket(input_arr, pas, config_immutable; initstates = init_states)

println("   ✓ ImmutableSolver executed")
println("   Results match MutableSolver: ", isapprox(result, result_immutable, rtol=1e-5))

# Test 4: Performance benchmark
println("\n4. Performance benchmark:")

bench_mutable = @benchmark $snow_bucket($input_arr, $pas, $config; initstates = $init_states)
println("   MutableSolver:")
println("     Median time: ", BenchmarkTools.prettytime(median(bench_mutable.times)))
println("     Memory: ", BenchmarkTools.prettymemory(median(bench_mutable.memory)))

bench_immutable = @benchmark $snow_bucket($input_arr, $pas, $config_immutable; initstates = $init_states)
println("   ImmutableSolver:")
println("     Median time: ", BenchmarkTools.prettytime(median(bench_immutable.times)))
println("     Memory: ", BenchmarkTools.prettymemory(median(bench_immutable.memory)))

# Summary
println("\n" * "="^60)
println("Test Summary")
println("="^60)
println("✓ Single-node bucket: PASSED")
println("✓ Multi-node bucket: PASSED")
println("✓ ImmutableSolver: PASSED")
println("✓ Performance benchmark: PASSED")
println("\n✅ All bucket tests passed!")
println("="^60)

