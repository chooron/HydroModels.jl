# `@stateflux_build` Macro Documentation

## Overview

The `@stateflux_build` macro provides a convenient way to create `StateFlux` objects directly from mathematical equations. This approach simplifies the definition of state variables and their dynamics in hydrological models.

## Syntax

```julia
# With a name
@stateflux_build name equation

# Without a name (auto-generated name)
@stateflux_build equation
```

Where:
- `name` (optional): A symbol representing the name of the flux (e.g., `:my_flux`)
- `equation`: An equation where:
  - Left side: The state variable
  - Right side: Expression describing the rate of change of the state variable

## Features

- Automatically identifies inputs and parameters from the equation
- Supports optional naming of the flux
- Creates a `StateFlux` object that can be incorporated into larger models

## Examples

### State Flux with a Name

```julia
using HydroModels
using ModelingToolkit

@variables S, P, ET, Q
@parameters a, b

# Create a state flux for a simple bucket model
# S is the state variable (e.g., water storage)
# P is precipitation input
# ET is evapotranspiration output
# Q is runoff output
flux = @stateflux_build :bucket_flux S = P - ET - Q
```

### State Flux with Parameters

```julia
using HydroModels
using ModelingToolkit

@variables S, P
@parameters a, b, c

# Create a state flux with parameters
# S is the state variable
# P is precipitation input
# a*S represents evapotranspiration as a function of storage
# b*S^c represents runoff as a function of storage
flux = @stateflux_build :param_bucket S = P - a*S - b*S^c
```

### State Flux without a Name

```julia
using HydroModels
using ModelingToolkit

@variables S, P, ET, Q

# Create a state flux with an auto-generated name
flux = @stateflux_build S = P - ET - Q
```

## Using State Fluxes

State fluxes represent the rate of change of state variables and are typically used in differential equations or as part of larger models. They cannot be run directly like `HydroFlux` or `NeuralFlux` objects.

```julia
# This would result in an error
# flux([1.0, 0.2, 0.3], ComponentVector())

# Instead, state fluxes are incorporated into larger models
# For example, in a differential equation:
# dS/dt = P - ET - Q
```

## Implementation Details

The macro works by:
1. Parsing the provided equation
2. Extracting the state variable from the left side
3. Extracting variables from the right side
4. Using ModelingToolkit's `isparameter` to identify parameters
5. Constructing a `StateFlux` object with the extracted information

## Comparison with Other Flux Macros

- `@hydroflux_build`: Creates a simple flux with explicit mathematical formulas for multiple outputs
- `@neuralflux`: Creates a neural network-based flux where the transformation is learned
- `@stateflux_build`: Creates a state flux that defines the rate of change of a state variable

Each macro serves a different purpose in hydrological modeling, allowing for flexible and expressive model construction.
