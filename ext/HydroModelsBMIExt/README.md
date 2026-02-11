# HydroModels BMI Extension

This extension provides Basic Model Interface (BMI) support for HydroModels.jl, enabling interoperability with other modeling frameworks.

## Overview

The BMI extension wraps HydroModels components (HydroBucket, HydroModel) with a standardized interface that allows:
- Time-stepping control
- Variable access and modification
- State management
- Grid information queries

## Installation

The BMI extension is automatically loaded when both HydroModels and BasicModelInterface are available:

```julia
using HydroModels
using BasicModelInterface
using Dates
```

## Quick Start

### 1. Create a Model

```julia
using HydroModels
using Symbolics
using ComponentArrays
using Dates

# Define symbolic variables
@variables temp lday prcp pet snowfall rainfall snowpack melt
@variables soilwater evap baseflow surfaceflow flow
@parameters Tmin Tmax Df Smax Qmax f

# Define step function
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

# Create snow bucket
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

# Create soil bucket
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

# Create complete model
model = @hydromodel :exphydro begin
    snow_bucket
    soil_bucket
end
```

### 2. Initialize BMI

```julia
using BasicModelInterface
const BMI = BasicModelInterface

# Prepare inputs
n_steps = 365
forcing = rand(3, n_steps)  # temp, lday, prcp
params = ComponentVector(
    Df = 2.67, Tmax = 0.18, Tmin = -2.09,
    Smax = 1709.0, f = 0.017, Qmax = 18.5
)
initstates = ComponentVector(snowpack = 0.0, soilwater = 1303.0)
timeidx = Date(2020, 1, 1):Day(1):Date(2020, 12, 31)
config = (solver = MutableSolver,)

# Initialize BMI model
bmi_model = BMI.initialize(
    model, forcing, params, initstates,
    collect(timeidx), config,
    lat=45.0, lon=-75.0  # Optional location
)
```

### 3. Run Model

```julia
# Single time step
BMI.update(bmi_model)

# Run to specific date
target_date = Date(2020, 6, 1)
output = BMI.update_until(bmi_model, target_date)

# Complete simulation
final_output = BMI.finalize(bmi_model)
```

### 4. Access Variables

```julia
# Get variable values
snowpack_values = zeros(100)
get_value(bmi_model, "snowpack", snowpack_values)

# Get pointer to values (no copy)
flow_ptr = get_value_ptr(bmi_model, "flow")

# Set variable values
new_snowpack = ones(100) * 50.0
set_value(bmi_model, "snowpack", new_snowpack)
```

## BMI Functions Reference

### Core Control Functions

- `BMI.initialize(component, forcing, pas, initstates, timeidx, config; lat, lon)` - Initialize model
- `BMI.update(model)` - Advance by one time step
- `BMI.update_until(model, time)` - Advance to specific time
- `BMI.finalize(model)` - Complete simulation and clean up

### Information Functions

- `get_component_name(model)` - Get model name
- `get_input_item_count(model)` - Number of input variables
- `get_output_item_count(model)` - Number of output variables
- `get_input_var_names(model)` - List of input variable names
- `get_output_var_names(model)` - List of output variable names

### Variable Information

- `get_var_grid(model, name)` - Grid ID for variable
- `get_var_type(model, name)` - Data type of variable
- `get_var_units(model, name)` - Units of variable
- `get_var_itemsize(model, name)` - Size of single item
- `get_var_nbytes(model, name)` - Total bytes for variable
- `get_var_location(model, name)` - Location (node/edge/face)

### Time Functions

- `get_current_time(model)` - Current simulation time
- `get_start_time(model)` - Start time
- `get_end_time(model)` - End time
- `get_time_units(model)` - Time units
- `get_time_step(model)` - Time step size

### Variable Access

- `get_value(model, name, dest)` - Get variable values
- `get_value_ptr(model, name)` - Get pointer to values
- `get_value_at_indices(model, name, dest, indices)` - Get values at indices
- `set_value(model, name, src)` - Set variable values
- `set_value_at_indices(model, name, indices, src)` - Set values at indices

### Grid Functions

For lumped models, grid functions return fixed values:
- `get_grid_rank(model, grid)` - Returns 2
- `get_grid_size(model, grid)` - Returns 1
- `get_grid_type(model, grid)` - Returns "uniform_rectilinear"
- `get_grid_shape(model, grid, shape)` - Returns (1, 1)
- `get_grid_origin(model, grid, origin)` - Returns (lat, lon) if provided

## Helper Functions

### State Management

```julia
# Save current state
saved_state = save_state(bmi_model)

# Run some steps
BMI.update(bmi_model)
BMI.update(bmi_model)

# Restore previous state
restore_state!(bmi_model, saved_state)
```

### Time Stepping Utilities

```julia
# Step forward multiple times
step_forward!(bmi_model, 10)  # Advance 10 steps

# Get time information
total_steps = get_time_step_count(bmi_model)
remaining = get_remaining_steps(bmi_model)
```

## Complete Example

```julia
using HydroModels
using BasicModelInterface
using ComponentArrays
using Dates

# Load model (from YAML or programmatic)
model = load_model_from_yaml("exphydro.yaml")

# Prepare data
forcing = rand(3, 365)
params = ComponentVector(Df=2.67, Tmax=0.18, Tmin=-2.09,
                         Smax=1709.0, f=0.017, Qmax=18.5)
initstates = ComponentVector(snowpack=0.0, soilwater=1303.0)
timeidx = Date(2020,1,1):Day(1):Date(2020,12,31)
config = (solver=MutableSolver,)

# Initialize BMI
bmi_model = BMI.initialize(model, forcing, params, initstates,
                           collect(timeidx), config)

# Run simulation with state saving
states_history = []
for i in 1:365
    # Save state every 30 days
    if i % 30 == 0
        push!(states_history, save_state(bmi_model))
    end

    # Update model
    BMI.update(bmi_model)

    # Access current values
    if i % 10 == 0
        flow = get_value_ptr(bmi_model, "flow")
        println("Day $i: Flow = $(flow[end])")
    end
end

# Finalize
final_output = BMI.finalize(bmi_model)
```

## Notes

### Current Limitations

1. **Lumped Models Only**: Current implementation supports only lumped (0-dimensional) models
2. **Date-based Time**: Time indices must be Date objects
3. **Memory Management**: Large forcing arrays are stored in memory

### Future Enhancements

1. **Distributed Models**: Support for 2D grid-based and network-based routing
2. **Flexible Time**: Support for different time representations
3. **Streaming Data**: Support for on-demand data loading

## Troubleshooting

### Error: "No outputs available"

This occurs when trying to access variables before running the model:

```julia
# Wrong
bmi_model = BMI.initialize(...)
flow = get_value_ptr(bmi_model, "flow")  # Error!

# Correct
bmi_model = BMI.initialize(...)
BMI.update(bmi_model)  # Run at least one step
flow = get_value_ptr(bmi_model, "flow")  # OK
```

### Error: "Cannot update: already at end time"

The model has reached the end of the time series:

```julia
# Check remaining steps
if get_remaining_steps(bmi_model) > 0
    BMI.update(bmi_model)
end
```

### Error: "Missing required parameter"

Ensure all model parameters are provided in the ComponentVector:

```julia
# Get required parameters
param_names = get_param_names(model)
println("Required parameters: ", param_names)

# Provide all parameters
params = ComponentVector(
    Df = 2.67,
    Tmax = 0.18,
    # ... all required parameters
)
```

## References

- [Basic Model Interface Specification](https://bmi.readthedocs.io/)
- [HydroModels.jl Documentation](https://github.com/chooron/HydroModels.jl)
- [ComponentArrays.jl](https://github.com/jonniedie/ComponentArrays.jl)
