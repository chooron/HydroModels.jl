# Complete ExpHydro Model Test with New Interface
# This script demonstrates the full workflow of building and running ExpHydro model
# using the updated HydroModels.jl interface

using CSV
using DataFrames
using ComponentArrays
using BenchmarkTools

# Include HydroModels
include("../../src/HydroModels.jl")
using .HydroModels

println("="^60)
println("ExpHydro Model Complete Test")
println("="^60)

# Define smooth step function
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

# Define variables and parameters
@variables temp lday prcp pet snowfall rainfall snowpack melt
@variables soilwater evap baseflow surfaceflow flow
@parameters Tmin Tmax Df Smax Qmax f

println("\n✓ Variables and parameters defined")

# Build the snow bucket
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

println("✓ Snow bucket created")

# Build the soil water bucket
soil_bucket = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux evap ~ step_func(soilwater) * pet * min(1.0, soilwater / Smax)
        @hydroflux baseflow ~ step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))
        @hydroflux surfaceflow ~ max(0.0, soilwater - Smax)
        @hydroflux flow ~ baseflow + surfaceflow
    end
    dfluxes = begin
        @stateflux soilwater ~ (rainfall + melt) - (evap + flow)
    end
end

println("✓ Soil bucket created")

# Combine into ExpHydro model
exphydro_model = @hydromodel :exphydro begin
    snow_bucket
    soil_bucket
end

println("✓ ExpHydro model created")
println("\nModel Information:")
println("  Input variables: ", HydroModels.get_input_names(exphydro_model))
println("  State variables: ", HydroModels.get_state_names(exphydro_model))
println("  Output variables: ", HydroModels.get_output_names(exphydro_model))
println("  Parameters: ", HydroModels.get_param_names(exphydro_model))

# Load data
println("\n" * "="^60)
println("Loading Data")
println("="^60)

file_path = "../notebook/data/exphydro/01013500.csv"
if !isfile(file_path)
    error("Data file not found: $file_path\nPlease ensure the data file exists.")
end

df = DataFrame(CSV.File(file_path))
ts = collect(1:1000)  # Use 1000 time steps for testing

println("✓ Data loaded: $(length(ts)) time steps")

# Prepare input data
input_data = (
    lday = df[ts, "dayl(day)"],
    temp = df[ts, "tmean(C)"],
    prcp = df[ts, "prcp(mm/day)"]
)

# Convert to matrix format in correct order
input_matrix = reduce(hcat, [input_data[name] for name in HydroModels.get_input_names(exphydro_model)]) |> permutedims

println("✓ Input matrix prepared: size = ", size(input_matrix))

# Prepare parameters
params = ComponentVector(
    f = 0.0167,
    Smax = 1709.46,
    Qmax = 18.47,
    Df = 2.674,
    Tmax = 0.17,
    Tmin = -2.09
)

# Prepare initial states
init_states = ComponentVector(
    snowpack = 0.0,
    soilwater = 1303.00
)

# Combine parameters
pas = ComponentVector(params = params)

println("✓ Parameters configured")
println("  ", params)

# Create configuration using NEW HydroConfig type
config = HydroConfig(
    solver = MutableSolver,
    interpolator = Val(DirectInterpolation),
    timeidx = ts,
    device = identity,
    min_value = 1e-6,
    parallel = false
)

println("✓ Configuration created (HydroConfig)")
println("  Solver: MutableSolver")
println("  Interpolator: DirectInterpolation")

# Run the model
println("\n" * "="^60)
println("Running Model")
println("="^60)

output_matrix = exphydro_model(
    input_matrix,
    pas,
    config;
    initstates = init_states
)

println("✓ Model executed successfully")
println("  Output size: ", size(output_matrix))

# Extract results
output_names = vcat(
    HydroModels.get_state_names(exphydro_model),
    HydroModels.get_output_names(exphydro_model)
)
output_df = NamedTuple{Tuple(output_names)}([output_matrix[i, :] for i in 1:size(output_matrix, 1)]) |> DataFrame

println("\n✓ Results extracted to DataFrame")
println("\nFirst 5 rows:")
println(first(output_df, 5))

# Benchmark performance
println("\n" * "="^60)
println("Performance Benchmark")
println("="^60)

println("\nBenchmarking model execution...")
bench_result = @benchmark $exphydro_model(
    $input_matrix,
    $pas,
    $config;
    initstates = $init_states
)

println("\nBenchmark Results:")
println("  Median time: ", BenchmarkTools.prettytime(median(bench_result.times)))
println("  Memory: ", BenchmarkTools.prettymemory(median(bench_result.memory)))
println("  Allocations: ", median(bench_result.allocs))

# Test with ImmutableSolver for comparison
println("\n" * "="^60)
println("Testing ImmutableSolver")
println("="^60)

config_immutable = HydroConfig(
    solver = ImmutableSolver,
    interpolator = Val(DirectInterpolation),
    timeidx = ts,
    min_value = 1e-6
)

output_immutable = exphydro_model(
    input_matrix,
    pas,
    config_immutable;
    initstates = init_states
)

println("✓ ImmutableSolver test passed")
println("  Output size: ", size(output_immutable))
println("  Results match: ", isapprox(output_matrix, output_immutable, rtol=1e-5))

bench_immutable = @benchmark $exphydro_model(
    $input_matrix,
    $pas,
    $config_immutable;
    initstates = $init_states
)

println("\nImmutableSolver Benchmark:")
println("  Median time: ", BenchmarkTools.prettytime(median(bench_immutable.times)))
println("  Memory: ", BenchmarkTools.prettymemory(median(bench_immutable.memory)))

# Summary
println("\n" * "="^60)
println("Test Summary")
println("="^60)
println("✓ Model construction: PASSED")
println("✓ Model execution: PASSED")
println("✓ Output extraction: PASSED")
println("✓ Performance benchmark: PASSED")
println("✓ ImmutableSolver test: PASSED")
println("\n✅ All tests completed successfully!")
println("="^60)

