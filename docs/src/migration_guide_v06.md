# 📚 Migration Guide: v0.5.x → v0.6.0

## Overview

HydroModels.jl v0.6.0 introduces several improvements focused on clarity, consistency, and performance. This guide helps you migrate your existing code from v0.5.x to v0.6.0.

**Key Changes:**
- Clearer interpolation method naming
- Consistent terminology for multi-node modeling
- Enhanced parameter handling flexibility
- Improved internal code organization

**Good News:** Most changes are backward compatible! Old names still work as aliases, giving you time to update your code gradually.

## Interpolation Changes

### What Changed

The interpolation method names have been updated to better reflect their actual behavior:

| Old Name (v0.5.x) | New Name (v0.6.0) | Description |
|-------------------|-------------------|-------------|
| `DirectInterpolation` | `ConstantInterpolation` | Ceiling-indexed constant interpolation |
| `EnzymeInterpolation` | `LinearInterpolation` | Linear interpolation between points |
| `EnzymeCompatibleInterpolation` | `LinearInterpolation` | Removed (use `LinearInterpolation`) |

### Why the Change?

The old names were confusing:
- "DirectInterpolation" didn't indicate it was constant interpolation
- "EnzymeInterpolation" suggested it was only for Enzyme, but both methods are Enzyme-compatible
- "EnzymeCompatibleInterpolation" was redundant

The new names clearly describe what each method does.

### Migration Steps

**Before (v0.5.x):**
```julia
using HydroModels

# Old configuration
config = (
    solver = MutableSolver,
    interpolator = Val(DirectInterpolation)
)

# Or with HydroConfig
config = HydroConfig(
    solver = MutableSolver,
    interpolator = Val(DirectInterpolation)
)
```

**After (v0.6.0):**
```julia
using HydroModels

# New configuration
config = HydroConfig(
    solver = MutableSolver,
    interpolator = Val(ConstantInterpolation)
)
```

### Backward Compatibility

⚠️ **Important:** Old names still work as aliases in v0.6.0, so your existing code won't break immediately. However, we recommend updating to the new names for clarity and future compatibility.

```julia
# These all work in v0.6.0 (backward compatible)
Val(DirectInterpolation)              # → ConstantInterpolation
Val(EnzymeInterpolation)              # → LinearInterpolation
Val(EnzymeCompatibleInterpolation)    # → LinearInterpolation
```

### Choosing the Right Interpolation Method

**Use `ConstantInterpolation` (default) when:**
- Working with daily timestep data
- Forcing data represents step functions
- Maximum performance is needed
- Data points represent period averages (e.g., daily rainfall)

**Use `LinearInterpolation` when:**
- Working with sub-daily timesteps
- Smooth transitions between forcing values are needed
- Higher accuracy is more important than speed
- Data points represent instantaneous measurements

See the [Interpolation Methods Guide](tutorials/interpolation_guide.md) for detailed comparison.

## Terminology Updates

### Multi-Node Modeling: `hru_types` → `htypes`

The parameter name for HRU (Hydrological Response Unit) type indices has been shortened for consistency:

**Before (v0.5.x):**
```julia
# Creating multi-node bucket
bucket = @hydrobucket :my_bucket begin
    fluxes = begin
        flux1
        flux2
    end
    dfluxes = begin
        state_flux1
    end
    hru_types = [1, 1, 2, 2, 3]  # Old name
end
```

**After (v0.6.0):**
```julia
# Creating multi-node bucket
bucket = @hydrobucket :my_bucket begin
    fluxes = begin
        flux1
        flux2
    end
    dfluxes = begin
        state_flux1
    end
    htypes = [1, 1, 2, 2, 3]  # New name
end
```

**Note:** This change applies to:
- `HydroBucket` constructor
- `@hydrobucket` macro
- `UnitHydrograph` constructor
- `@unithydro` macro
- `NeuralBucket` constructor

## Parameter Handling Improvements

### AbstractVector Support

Components now accept `AbstractVector` instead of requiring `ComponentVector`, providing more flexibility:

**Before (v0.5.x):**
```julia
# Had to use ComponentVector
params = ComponentVector(params=(k=0.5, n=2.0), nns=())
output = bucket(input, params, config)
```

**After (v0.6.0):**
```julia
# Can use ComponentVector (recommended)
params = ComponentVector(params=(k=0.5, n=2.0), nns=())
output = bucket(input, params, config)

# Or plain Vector with proper structure
params_vec = [0.5, 2.0]  # Must match expected parameter order
output = bucket(input, params_vec, config)
```

**Recommendation:** Continue using `ComponentVector` for clarity and named parameter access. Plain vectors are supported for advanced use cases and optimization frameworks.

## File Structure Changes

### Internal Reorganization

The internal file structure has been reorganized for better maintainability:

| Old File | New Files | Impact |
|----------|-----------|--------|
| `src/tools.jl` | `src/interpolate.jl` + `src/solve.jl` | Internal only - no user impact |

**User Impact:** None. These are internal changes that don't affect the public API. All exported functions remain available through `using HydroModels`.

## Testing Your Migration

Follow these steps to verify your code works with v0.6.0:

### 1. Update Package

```julia
using Pkg
Pkg.update("HydroModels")
```

### 2. Check Version

```julia
using HydroModels
println(pkgversion(HydroModels))  # Should show 0.6.0 or higher
```

### 3. Run Your Tests

```julia
# Run your existing test suite
using Pkg
Pkg.test("YourPackage")
```

### 4. Update Terminology (Recommended)

Search and replace in your codebase:
- `DirectInterpolation` → `ConstantInterpolation`
- `EnzymeInterpolation` → `LinearInterpolation`
- `EnzymeCompatibleInterpolation` → `LinearInterpolation`
- `hru_types` → `htypes` (in component constructors)

### 5. Verify Functionality

Run a simple test to ensure everything works:

```julia
using HydroModels

# Define a simple model
@variables temp prcp flow
@parameters k

flux = @hydroflux flow ~ k * prcp
bucket = @hydrobucket :test begin
    fluxes = begin
        flux
    end
end

# Test with new configuration
config = HydroConfig(
    solver = MutableSolver,
    interpolator = Val(ConstantInterpolation)
)

# Create test data
input = rand(1, 100)  # 1 variable, 100 timesteps
params = ComponentVector(params=(k=0.5,), nns=())

# Run model
output = bucket(input, params, config)

println("✅ Migration successful! Output size: ", size(output))
```

## Common Issues and Solutions

### Issue: "DirectInterpolation not found"

**Cause:** You're using the old name without the backward compatibility aliases.

**Solution:** Update to `ConstantInterpolation` or ensure you're using v0.6.0+ which includes the aliases.

### Issue: "hru_types not recognized"

**Cause:** Using old parameter name in v0.6.0+.

**Solution:** Change `hru_types` to `htypes` in your component constructors.

### Issue: Performance regression after update

**Cause:** Unlikely, but if you notice performance changes, check your interpolation method.

**Solution:** Ensure you're using the appropriate interpolation method for your use case. `ConstantInterpolation` is faster than `LinearInterpolation`.

## Getting Help

If you encounter issues during migration:

1. **Check the documentation:** [HydroModels.jl Docs](https://chooron.github.io/HydroModels.jl/)
2. **Review examples:** See the `examples/` directory in the repository
3. **Ask for help:** Open an issue on [GitHub](https://github.com/chooron/HydroModels.jl/issues)

## Summary

**Required Changes:**
- None! Backward compatibility is maintained.

**Recommended Changes:**
- Update `DirectInterpolation` → `ConstantInterpolation`
- Update `EnzymeInterpolation` → `LinearInterpolation`
- Update `hru_types` → `htypes`

**Benefits of Updating:**
- Clearer, more descriptive code
- Better alignment with future versions
- Improved code readability

Take your time with the migration - the backward compatibility ensures your code will continue working while you update at your own pace.
