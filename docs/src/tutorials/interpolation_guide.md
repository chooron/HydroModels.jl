# 🔄 Interpolation Methods Guide

## Introduction

Interpolation is a critical component in hydrological modeling that determines how input forcing data (precipitation, temperature, etc.) is handled between discrete timesteps. The choice of interpolation method affects both model accuracy and computational performance.

HydroModels.jl provides two built-in interpolation methods, both fully compatible with automatic differentiation frameworks (Enzyme and Zygote), making them suitable for parameter optimization and sensitivity analysis.

## Available Methods

### ConstantInterpolation

**Algorithm:** Ceiling-indexed constant interpolation

**How it works:**
- For any time `t` between timesteps, uses the value at `ceil(t)`
- Creates a step-function representation of the forcing data
- No interpolation between points - values change instantaneously at timestep boundaries

**Mathematical representation:**
```
f(t) = data[ceil(t)]  for t ∈ [1, n]
```

**Use when:**
- ✅ Working with daily timestep data
- ✅ Forcing data represents period averages (e.g., daily rainfall totals)
- ✅ Maximum computational performance is needed
- ✅ Step-function behavior is physically appropriate

**Example:**
```julia
using HydroModels

config = HydroConfig(
    solver = MutableSolver,
    interpolator = Val(ConstantInterpolation)
)

# Run model with constant interpolation
output = model(input, params, config)
```

**Performance characteristics:**
- **Speed:** Fastest (simple array indexing)
- **Memory:** Minimal overhead
- **AD compatibility:** Full support for Enzyme and Zygote

### LinearInterpolation

**Algorithm:** Linear interpolation between adjacent data points

**How it works:**
- For any time `t` between timesteps `i` and `i+1`, computes weighted average
- Creates smooth transitions between forcing values
- Interpolation weight based on fractional position between points

**Mathematical representation:**
```
f(t) = data[i] * (1 - α) + data[i+1] * α
where i = floor(t), α = t - i
```

**Use when:**
- ✅ Working with sub-daily timesteps (hourly, 3-hourly)
- ✅ Forcing data represents instantaneous measurements
- ✅ Smooth transitions between values are physically important
- ✅ Higher accuracy is more important than speed

**Example:**
```julia
using HydroModels

config = HydroConfig(
    solver = ImmutableSolver,  # Often paired with LinearInterpolation for AD
    interpolator = Val(LinearInterpolation)
)

# Run model with linear interpolation
output = model(input, params, config)
```

**Performance characteristics:**
- **Speed:** Slightly slower than ConstantInterpolation (requires arithmetic operations)
- **Memory:** Minimal overhead
- **AD compatibility:** Full support for Enzyme and Zygote

## Comparison

### Performance Benchmark

Typical performance comparison on a standard hydrological model (1000 timesteps, single node):

| Method | Relative Speed | Memory Usage | AD Overhead |
|--------|---------------|--------------|-------------|
| ConstantInterpolation | 1.0x (baseline) | Minimal | ~1.5x |
| LinearInterpolation | 0.95x | Minimal | ~1.6x |

**Note:** Performance differences are typically <5% for most applications. Choose based on physical appropriateness rather than performance alone.

### Accuracy Comparison

**For daily timestep models:**
- Both methods produce similar results
- ConstantInterpolation is physically appropriate for daily rainfall
- Differences typically <1% in final outputs

**For sub-daily timestep models:**
- LinearInterpolation provides smoother forcing transitions
- Can improve accuracy by 5-15% for sub-daily simulations
- More physically realistic for temperature and radiation data

### When to Use Each Method

```
Decision Tree:

Is your timestep daily or longer?
├─ Yes → Use ConstantInterpolation
└─ No (sub-daily)
   ├─ Is smooth forcing transition important?
   │  ├─ Yes → Use LinearInterpolation
   │  └─ No → Use ConstantInterpolation
   └─ Are you doing parameter optimization?
      ├─ Yes → Consider LinearInterpolation (smoother gradients)
      └─ No → Either method works
```

## Advanced: DataInterpolations.jl Integration

For more sophisticated interpolation needs, HydroModels.jl integrates with [DataInterpolations.jl](https://github.com/SciML/DataInterpolations.jl) through an extension.

**Available methods:**
- Cubic spline interpolation
- Akima interpolation
- B-spline interpolation
- And more...

**Example:**
```julia
using HydroModels
using DataInterpolations

# Create custom interpolator
custom_interp = CubicSpline(data, time_points)

# Use in configuration
config = HydroConfig(
    solver = MutableSolver,
    interpolator = custom_interp
)

output = model(input, params, config)
```

**Note:** Advanced interpolators may have limited AD support. Check DataInterpolations.jl documentation for compatibility.

## Enzyme Compatibility

Both `ConstantInterpolation` and `LinearInterpolation` are designed to be fully compatible with Enzyme.jl for automatic differentiation:

**Key features:**
- No dynamic memory allocation in hot loops
- Type-stable implementations
- Explicit handling of interpolation logic

**Example with Enzyme:**
```julia
using HydroModels
using Enzyme

# Define model with interpolation
function run_model(params_vec)
    config = HydroConfig(interpolator = Val(LinearInterpolation))
    output = model(input, params_vec, config)
    return sum(output)  # Scalar objective
end

# Compute gradient
grad = Enzyme.gradient(Reverse, run_model, params)
```

## Practical Examples

### Example 1: Daily Rainfall-Runoff Model

```julia
using HydroModels

# Daily data - use ConstantInterpolation
@variables prcp temp pet flow
@parameters k n

# Define model components
flux = @hydroflux flow ~ k * prcp
bucket = @hydrobucket :daily_model begin
    fluxes = begin
        flux
    end
end

# Configuration for daily timesteps
config = HydroConfig(
    solver = MutableSolver,
    interpolator = Val(ConstantInterpolation),  # Appropriate for daily data
    timeidx = 1:365
)

# Run simulation
output = bucket(daily_input, params, config)
```

### Example 2: Hourly Energy Balance Model

```julia
using HydroModels

# Hourly data - use LinearInterpolation for smooth transitions
@variables temp radiation snowmelt
@parameters melt_factor

flux = @hydroflux snowmelt ~ melt_factor * max(temp, 0) * radiation
bucket = @hydrobucket :hourly_model begin
    fluxes = begin
        flux
    end
end

# Configuration for hourly timesteps
config = HydroConfig(
    solver = ImmutableSolver,  # Better for AD
    interpolator = Val(LinearInterpolation),  # Smooth forcing transitions
    timeidx = 1:8760  # One year of hourly data
)

# Run simulation
output = bucket(hourly_input, params, config)
```

### Example 3: Parameter Optimization with Interpolation

```julia
using HydroModels
using Optimization, OptimizationOptimJL

# Define objective function
function objective(params_vec, p)
    config = HydroConfig(
        solver = ImmutableSolver,
        interpolator = Val(LinearInterpolation)  # Smoother gradients
    )

    simulated = model(input, params_vec, config)
    observed = p.observed

    # Nash-Sutcliffe Efficiency
    return -nse(simulated, observed)
end

# Setup optimization
prob = OptimizationProblem(objective, initial_params, observed_data)
sol = solve(prob, BFGS())
```

## Troubleshooting

### Issue: Unexpected discontinuities in output

**Cause:** Using ConstantInterpolation with sub-daily timesteps

**Solution:** Switch to LinearInterpolation for smoother forcing transitions

```julia
# Before (discontinuous)
config = HydroConfig(interpolator = Val(ConstantInterpolation))

# After (smooth)
config = HydroConfig(interpolator = Val(LinearInterpolation))
```

### Issue: Slow gradient computation

**Cause:** Using complex interpolator from DataInterpolations.jl

**Solution:** Use built-in LinearInterpolation for better AD performance

```julia
# Before (slow)
config = HydroConfig(interpolator = CubicSpline(data, times))

# After (fast)
config = HydroConfig(interpolator = Val(LinearInterpolation))
```

### Issue: "Interpolator not found" error

**Cause:** Typo in interpolator name or using old names

**Solution:** Use correct names: `ConstantInterpolation` or `LinearInterpolation`

```julia
# Wrong (old name)
config = HydroConfig(interpolator = Val(DirectInterpolation))

# Correct (new name)
config = HydroConfig(interpolator = Val(ConstantInterpolation))
```

## Summary

**Quick Reference:**

| Scenario | Recommended Method | Reason |
|----------|-------------------|--------|
| Daily rainfall-runoff | ConstantInterpolation | Period averages, step function appropriate |
| Hourly temperature | LinearInterpolation | Smooth transitions, instantaneous values |
| Parameter optimization | LinearInterpolation | Smoother gradients |
| Maximum performance | ConstantInterpolation | Fastest execution |
| Sub-daily energy balance | LinearInterpolation | Physical realism |

**Key Takeaways:**
- Both methods are fast and AD-compatible
- Choose based on physical appropriateness, not just performance
- ConstantInterpolation is default and works well for daily models
- LinearInterpolation provides smoother behavior for sub-daily models
- Performance difference is typically <5%

## Further Reading

- [Configuration Guide](../concepts_en.md#configuration) - Complete configuration options
- [Migration Guide](../migration_guide_v06.md) - Upgrading from v0.5.x
- [Getting Started](../get_start_en.md) - Basic usage examples
- [DataInterpolations.jl](https://github.com/SciML/DataInterpolations.jl) - Advanced interpolation methods
