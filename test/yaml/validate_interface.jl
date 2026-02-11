#!/usr/bin/env julia

"""
Quick validation script to check if the updated interface works correctly.

This script tests basic functionality without running the full test suite.
"""

println("=" ^ 70)
println("HydroModels.jl Interface Validation")
println("=" ^ 70)

using Pkg
Pkg.activate(".")

println("\n[1/6] Loading packages...")
using HydroModels
using ComponentArrays
using Statistics

println("✓ Packages loaded successfully")

println("\n[2/6] Testing interpolation aliases...")
@assert isdefined(HydroModels, :ConstantInterpolation)
@assert isdefined(HydroModels, :LinearInterpolation)
@assert isdefined(HydroModels, :DirectInterpolation)
@assert HydroModels.DirectInterpolation === HydroModels.ConstantInterpolation
println("✓ Interpolation types available")

println("\n[3/6] Testing single-node flux...")
@variables a b c
@parameters p1 p2

simple_flux = @hydroflux c ~ a * p1 + b * p2
input_2d = [2.0 3.0 1.0; 3.0 2.0 2.0]
params = ComponentVector(params=(p1=3.0, p2=4.0))
output = simple_flux(input_2d, params)

@assert size(output) == (1, 3)
@assert output ≈ [18.0 17.0 11.0]
println("✓ Single-node flux works correctly")

println("\n[4/6] Testing multi-node flux...")
multi_flux = @hydroflux begin
    c ~ a * p1 + b * p2
    htypes = [1, 2, 3]
end

input_3d = repeat(reshape(input_2d, 2, 1, 3), 1, 3, 1)
multi_params = ComponentVector(params=(p1=[3.0, 3.0, 3.0], p2=[4.0, 4.0, 4.0]))
multi_output = multi_flux(input_3d, multi_params)

@assert size(multi_output) == (1, 3, 3)
println("✓ Multi-node flux works correctly")

println("\n[5/6] Testing single-node bucket...")
@variables temp prcp snowfall rainfall snowpack
@parameters Tmin

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

snow_bucket = @hydrobucket :snow begin
    fluxes = begin
        @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
        @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
    end
    dfluxes = begin
        @stateflux snowpack ~ snowfall - rainfall
    end
end

input_2d = [5.0 -2.0 10.0; 10.0 15.0 5.0]  # temp, prcp
params = ComponentVector(params=(Tmin=0.0,))
config = HydroConfig(solver=MutableSolver, timeidx=1:3)
output = snow_bucket(input_2d, params, config; initstates=ComponentVector(snowpack=0.0))

@assert size(output) == (3, 3)  # snowpack + snowfall + rainfall
println("✓ Single-node bucket works correctly")

println("\n[6/6] Testing multi-node bucket...")
multi_bucket = @hydrobucket :snow begin
    fluxes = begin
        @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
        @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
    end
    dfluxes = begin
        @stateflux snowpack ~ snowfall - rainfall
    end
    htypes = [1, 1, 2]
end

input_3d = repeat(reshape(input_2d, 2, 1, 3), 1, 3, 1)
multi_params = ComponentVector(params=(Tmin=[0.0, -5.0],))
multi_states = [0.0, 0.0, 0.0]  # Flattened states for 3 nodes
multi_output = multi_bucket(input_3d, multi_params, config; initstates=multi_states)

@assert size(multi_output) == (3, 3, 3)  # (states+outputs) × nodes × time
println("✓ Multi-node bucket works correctly")

println("\n" * "=" ^ 70)
println("✓ All validation tests passed!")
println("=" ^ 70)
println("\nThe updated interface is working correctly.")
println("You can now run the full test suite with:")
println("  julia test/run_tests.jl")
println("  or")
println("  using Pkg; Pkg.test(\"HydroModels\")")
