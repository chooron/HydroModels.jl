# `@neuralflux` Macro Documentation

## Overview

The `@neuralflux` macro provides a convenient way to create `NeuralFlux` objects using a clean, intuitive syntax. It simplifies the integration of neural networks into hydrological models by allowing a declarative syntax that clearly shows the relationship between inputs and outputs.

## Syntax

```julia
@neuralflux output ~ chain([inputs])
```

Where:
- `output`: A single output variable or a vector of output variables
- `~`: The separator (tilde) indicating a neural network transformation
- `chain`: A Lux neural network chain with a name property
- `inputs`: A vector of input variables enclosed in square brackets

## Features

- Clean, intuitive syntax that clearly shows the transformation from inputs to outputs
- Supports both single and multiple output variables
- Integrates seamlessly with Lux neural network chains
- Automatically handles variable type conversions

## Examples

### Single Output

```julia
using HydroModels
using ModelingToolkit
using Lux

@variables x, y, z
chain = Chain(Dense(2 => 10, relu), Dense(10 => 1), name=:my_net)

# Create a neural flux with a single output
flux = @neuralflux z ~ chain([x, y])
```

### Multiple Outputs

```julia
using HydroModels
using ModelingToolkit
using Lux

@variables x, y, z₁, z₂
chain = Chain(Dense(2 => 16, relu), Dense(16 => 2), name=:multi_net)

# Create a neural flux with multiple outputs
flux = @neuralflux [z₁, z₂] ~ chain([x, y])
```

## Using the Created Flux

After creating a flux with the `@neuralflux` macro, you can use it like any other `NeuralFlux` object:

```julia
# Define input values
input_values = [2.0, 3.0]  # Values for x and y

# Initialize parameters for the neural network
using Random
rng = StableRNG(42)
params = ComponentVector(Lux.initialparameters(rng, chain))

# Create ComponentVector for parameters
using ComponentArrays
params_cv = ComponentVector(my_net=params)

# Calculate the flux
result = flux(input_values, params_cv)
```

## Implementation Details

The macro works by:
1. Parsing the provided expression with the tilde (`~`) operator
2. Extracting the output variables from the left side
3. Extracting the neural network chain and input variables from the right side
4. Constructing a `NeuralFlux` object with the extracted information

## Important Notes

1. The neural network chain must have a name property set (e.g., `name=:my_net`)
2. The number of inputs to the neural network should match the number of input variables provided
3. The number of outputs from the neural network should match the number of output variables provided
4. When using the flux, the parameters ComponentVector must have a field with the same name as the neural network chain
5. Input variables must be provided as a vector enclosed in square brackets: `chain([x, y])`
