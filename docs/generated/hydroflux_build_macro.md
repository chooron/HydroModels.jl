# `@hydroflux` Macro Documentation

## Overview

The `@hydroflux` macro provides a convenient way to create `HydroFlux` objects directly from mathematical equations. This approach is more intuitive and less error-prone than manually specifying inputs, outputs, parameters, and expressions separately.

## Syntax

```julia
# With a name
@hydroflux name equations

# Without a name (auto-generated name)
@hydroflux equations
```

Where:
- `name` (optional): A symbol representing the name of the flux (e.g., `:my_flux`)
- `equations`: One or more equations where:
  - Left side: Output variables
  - Right side: Expressions containing input variables and parameters

## Features

- Automatically identifies inputs, outputs, and parameters from the equations
- Supports multiple equations in a single macro call
- Automatically detects parameters using ModelingToolkit's `isparameter` function
- Handles variables that appear on both sides of equations correctly
- Supports optional naming of the flux

## Examples

### Single Equation with Named Flux

```julia
using HydroModels
using ModelingToolkit

@variables x, y, z
@parameters a, b

# Create a simple flux that calculates z = a*x + b*y
flux = @hydroflux :linear_flux z = a*x + b*y
```

### Multiple Equations with Named Flux

```julia
using HydroModels
using ModelingToolkit

@variables x, y, z₁, z₂
@parameters a, b, c

# Create a flux with multiple outputs
flux = @hydroflux :multi_output_flux begin
    z₁ = a*x + b*y
    z₂ = c*x*y
end
```

### Single Equation without a Name

```julia
using HydroModels
using ModelingToolkit

@variables x, y, z
@parameters a, b

# Create a flux with an auto-generated name
flux = @hydroflux z = a*x + b*y
```

### Multiple Equations without a Name

```julia
using HydroModels
using ModelingToolkit

@variables x, y, z₁, z₂
@parameters a, b, c

# Create a flux with multiple outputs and an auto-generated name
flux = @hydroflux begin
    z₁ = a*x + b*y
    z₂ = c*x*y
end
```

## Using the Created Flux

After creating a flux with the `@hydroflux` macro, you can use it like any other `HydroFlux` object:

```julia
# Define input values
input_values = [2.0, 3.0]  # Values for x and y

# Define parameter values using ComponentArrays
using ComponentArrays
params = ComponentVector(a=0.5, b=1.5, c=2.0)

# Calculate the flux
result = flux(input_values, params)
```

## Advanced Usage

### Variables on Both Sides

If a variable appears on both sides of an equation, it will be correctly handled:

```julia
@variables x, y, z
@parameters a, b

# z appears on both sides, but will be treated as an output
flux = @hydroflux z = a*x + b*y + 0.1*z
```

### Complex Expressions

The macro supports complex mathematical expressions:

```julia
@variables x, y, z
@parameters a, b, c

# Complex mathematical expression
flux = @hydroflux z = a*sin(x) + b*exp(c*y)
```

## Implementation Details

The macro works by:
1. Parsing the provided equations
2. Extracting variables from both sides
3. Using ModelingToolkit's `isparameter` to identify parameters
4. Removing output variables from the input list (if they appear on both sides)
5. Constructing a `HydroFlux` object with the extracted information
