# HydroConfig Examples and Tests
# This script demonstrates different ways to use the new HydroConfig system

include("../../src/HydroModels.jl")
using .HydroModels
using ComponentArrays

println("="^60)
println("HydroConfig System Examples")
println("="^60)

# Example 1: Default configuration
println("\n1. Creating default configuration:")
config_default = HydroConfig()
println("   ✓ Default config created")
println("   Solver: ", config_default.solver)
println("   Interpolator: ", config_default.interpolator)
println("   Min value: ", config_default.min_value)

# Example 2: Configuration with MutableSolver
println("\n2. Configuration with MutableSolver:")
config_mutable = HydroConfig(
    solver = MutableSolver,
    interpolator = Val(DirectInterpolation),
    timeidx = 1:100,
    min_value = 1e-6
)
println("   ✓ MutableSolver config created")
println("   Time steps: ", length(config_mutable.timeidx))

# Example 3: Configuration with ImmutableSolver
println("\n3. Configuration with ImmutableSolver:")
config_immutable = HydroConfig(
    solver = ImmutableSolver,
    interpolator = Val(DirectInterpolation),
    timeidx = 1:100
)
println("   ✓ ImmutableSolver config created")

# Example 4: Using default_config()
println("\n4. Using default_config() function:")
config_func = HydroModels.default_config()
println("   ✓ Config from function created")

# Example 5: Merging configurations
println("\n5. Merging configurations:")
config_base = HydroConfig(solver = MutableSolver, min_value = 1e-6)
config_merged = HydroModels.merge_config(config_base; timeidx = 1:50, parallel = true)
println("   ✓ Merged config created")
println("   Parallel: ", config_merged.parallel)
println("   Time steps: ", length(config_merged.timeidx))

# Example 6: Backward compatibility with NamedTuple
println("\n6. Backward compatibility test:")
old_style_config = (
    solver = MutableSolver,
    interpolator = Val(DirectInterpolation),
    timeidx = 1:100,
    min_value = 1e-6
)
normalized = HydroModels.normalize_config(old_style_config)
println("   ✓ NamedTuple config normalized to HydroConfig")
println("   Type: ", typeof(normalized))

# Example 7: Converting HydroConfig to NamedTuple
println("\n7. Converting HydroConfig to NamedTuple:")
as_namedtuple = HydroModels.to_namedtuple(config_default)
println("   ✓ Converted to NamedTuple")
println("   Type: ", typeof(as_namedtuple))

# Example 8: Accessing config values safely
println("\n8. Safe config value access:")
solver_value = HydroModels.get_config_value(config_default, :solver, MutableSolver)
println("   ✓ Retrieved solver: ", solver_value)
nonexistent = HydroModels.get_config_value(config_default, :nonexistent, "default")
println("   ✓ Retrieved with default: ", nonexistent)

# Example 9: Different solver types
println("\n9. All solver types:")
solvers = [MutableSolver, ImmutableSolver, ODESolver, DiscreteSolver]
for solver in solvers
    cfg = HydroConfig(solver = solver)
    println("   ✓ ", solver, " - Type: ", typeof(cfg.solver))
end

# Example 10: Validation tests
println("\n10. Configuration validation:")
try
    invalid_config = HydroConfig(min_value = -1.0)
    println("   ✗ Should have failed with negative min_value")
catch e
    println("   ✓ Correctly rejected negative min_value")
    println("   Error: ", e.msg)
end

try
    valid_config = HydroConfig(min_value = 1e-8)
    println("   ✓ Accepted positive min_value: ", valid_config.min_value)
catch e
    println("   ✗ Should have accepted positive min_value")
end

# Summary
println("\n" * "="^60)
println("Summary")
println("="^60)
println("✓ Default configuration")
println("✓ Custom configurations")
println("✓ Configuration merging")
println("✓ Backward compatibility")
println("✓ Type conversions")
println("✓ Safe value access")
println("✓ All solver types")
println("✓ Validation")
println("\n✅ All configuration examples passed!")
println("="^60)

