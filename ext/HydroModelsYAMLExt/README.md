# HydroModelsYAMLExt

Extension module for loading HydroModels from YAML configuration files.

## Overview

This extension provides functionality to define hydrological models using YAML configuration files instead of Julia code. It automatically loads when both `HydroModels` and `YAML` packages are available.

## Installation

```julia
using Pkg
Pkg.add("HydroModels")
Pkg.add("YAML")
```

The extension will be automatically loaded when you use both packages:

```julia
using HydroModels
using YAML  # This triggers the extension
```

## Features

- **Declarative Model Definition**: Define models using human-readable YAML syntax
- **Parameter Metadata**: Include descriptions, units, bounds, and default values
- **Automatic Variable Detection**: Variables and parameters are automatically extracted from formulas
- **Configuration Support**: Specify solver, interpolator, and other settings in YAML
- **Type Safety**: Full integration with HydroModels' type system

## Usage

### Basic Example

```julia
using HydroModels
using YAML

# Load model from YAML file
model = load_model_from_yaml("my_model.yaml")

# Load configuration
config = load_config_from_yaml("my_model.yaml")

# Load parameter metadata
params_info = load_parameters_from_yaml("my_model.yaml")

# Run model
output = model(forcing_data, parameters, config)
```

## YAML Format

### Complete Example

```yaml
version: "1.0"
schema: "hydromodels"

# Parameter definitions with metadata
parameters:
  Smax:
    description: "Maximum soil moisture storage"
    units: "mm"
    default: 100.0
    bounds: [50.0, 500.0]

  f:
    description: "Exponential decay parameter"
    units: "dimensionless"
    default: 2.0
    bounds: [0.5, 5.0]

  Qmax:
    description: "Maximum baseflow rate"
    units: "mm/day"
    default: 10.0
    bounds: [1.0, 50.0]

# Component definitions
components:
  - type: HydroBucket
    name: soil_bucket

    # Flux calculations
    fluxes:
      - name: baseflow
        formula: "baseflow ~ Qmax * exp(-f * max(0.0, Smax - soilwater) / Smax)"

      - name: surface_runoff
        formula: "runoff ~ max(0.0, rainfall - infiltration)"

    # State variable derivatives
    state_fluxes:
      - name: soil_balance
        formula: "soilwater ~ rainfall - evap - baseflow - runoff"

    # Optional: HRU types for spatial models
    hru_types: [1, 1, 2, 2]

# Model composition
model:
  type: HydroModel
  name: simple_hydro_model
  components: [soil_bucket]

# Configuration (optional)
config:
  solver: MutableSolver
  interpolator: DirectInterpolation
  min_value: 1.0e-6
  parallel: false
```

### YAML Structure

#### 1. Parameters Section

Define model parameters with metadata:

```yaml
parameters:
  parameter_name:
    description: "Human-readable description"
    units: "Physical units"
    default: 100.0
    bounds: [min_value, max_value]
```

#### 2. Components Section

Define model components (currently supports `HydroBucket`):

```yaml
components:
  - type: HydroBucket
    name: component_name

    fluxes:
      - name: flux_name
        formula: "output_var ~ mathematical_expression"

    state_fluxes:
      - name: state_name
        formula: "state_var ~ rate_of_change_expression"

    hru_types: [1, 1, 2]  # Optional, for multi-node models
```

#### 3. Model Section

Compose components into a complete model:

```yaml
model:
  type: HydroModel
  name: model_name
  components: [component1, component2]
```

#### 4. Config Section (Optional)

Specify runtime configuration:

```yaml
config:
  solver: MutableSolver  # or ImmutableSolver, ODESolver, DiscreteSolver
  interpolator: DirectInterpolation  # or EnzymeCompatibleInterpolation
  min_value: 1.0e-6
  parallel: false
```

## Formula Syntax

Formulas use standard mathematical notation:

### Operators
- Arithmetic: `+`, `-`, `*`, `/`, `^`
- Comparison: `<`, `>`, `<=`, `>=`, `==`

### Functions
- Basic: `min`, `max`, `abs`, `sqrt`
- Exponential: `exp`, `log`
- Trigonometric: `sin`, `cos`, `tan`, `tanh`
- Rounding: `ceil`, `floor`, `round`
- Custom: `step_func`, `smoothlogistic_func` (from HydroModels)

### Examples

```yaml
# Simple linear relationship
formula: "Q ~ k * S"

# Nonlinear storage-discharge
formula: "Q ~ Qmax * (S / Smax)^f"

# Conditional logic using max/min
formula: "infiltration ~ min(rainfall, Smax - S)"

# Complex expression
formula: "ET ~ PET * (1.0 - exp(-S / Smax))"
```

## API Reference

### Main Functions

#### `load_model_from_yaml(yaml_file::String)`

Load a complete hydrological model from YAML file.

**Arguments:**
- `yaml_file`: Path to YAML configuration file

**Returns:**
- `HydroModel` object ready for simulation

**Example:**
```julia
model = load_model_from_yaml("exphydro.yaml")
```

#### `load_config_from_yaml(yaml_file::String)`

Load configuration from YAML file.

**Arguments:**
- `yaml_file`: Path to YAML configuration file

**Returns:**
- `HydroConfig` object

**Example:**
```julia
config = load_config_from_yaml("exphydro.yaml")
```

#### `load_parameters_from_yaml(yaml_file::String)`

Load parameter metadata from YAML file.

**Arguments:**
- `yaml_file`: Path to YAML configuration file

**Returns:**
- Dictionary mapping parameter symbols to metadata dictionaries

**Example:**
```julia
params_info = load_parameters_from_yaml("exphydro.yaml")
println(params_info[:Smax][:description])  # "Maximum soil moisture storage"
println(params_info[:Smax][:bounds])       # (50.0, 500.0)
```

## Advantages of YAML Configuration

1. **Readability**: Non-programmers can understand and modify models
2. **Version Control**: Easy to track changes in model structure
3. **Documentation**: Built-in parameter descriptions and metadata
4. **Portability**: Models can be shared across different projects
5. **Validation**: Parameter bounds and units are explicitly defined
6. **Separation of Concerns**: Model structure separate from implementation

## Limitations

- Currently only supports `HydroBucket` components
- Formula parsing is limited to standard mathematical expressions
- No support for custom Julia functions in formulas (use functional interface instead)
- Variables and parameters must be valid Julia identifiers

## Future Extensions

Planned features for future versions:

- [ ] Support for `HydroRoute` components
- [ ] Support for `UnitHydrograph` components
- [ ] Support for `NeuralFlux` and `NeuralBucket` components
- [ ] Schema validation
- [ ] YAML to Julia code generation
- [ ] Interactive YAML editor/validator

## Examples

See the `examples/` directory for complete YAML model examples:

- `exphydro.yaml`: ExpHydro conceptual rainfall-runoff model
- (Add more examples as they are created)

## Technical Details

### Extension Loading

This module uses Julia's package extension system (Julia 1.9+). The extension is automatically loaded when both parent packages are available:

```julia
# In Project.toml
[extensions]
HydroModelsYAMLExt = "YAML"
```

### Variable and Parameter Detection

The extension automatically:
1. Extracts all symbols from formulas using regex
2. Distinguishes between parameters (defined in `parameters` section) and variables
3. Creates Symbolics variables and parameters
4. Builds symbolic expressions for model compilation

### Type Safety

All components created from YAML are fully type-stable and compatible with:
- Automatic differentiation (Zygote, Enzyme)
- GPU computation
- Optimization frameworks
- All HydroModels features

## Troubleshooting

### Extension Not Loading

**Problem**: Functions like `load_model_from_yaml` are not available.

**Solution**: Make sure both packages are loaded:
```julia
using HydroModels
using YAML
```

### Formula Parsing Errors

**Problem**: Error message "Failed to parse expression".

**Solution**: Check that:
- All variables are used in formulas
- Function names are spelled correctly
- Parentheses are balanced
- The `~` separator is present in formulas

### Unknown Variable Errors

**Problem**: "Unknown variable in formula LHS".

**Solution**: Make sure the left-hand side variable is either:
- Defined in another flux formula (for flux outputs)
- A state variable (for state flux formulas)

## Contributing

To contribute to this extension:

1. Add new component builders in `component_builder.jl`
2. Extend formula parser in `expression_parser.jl` for new syntax
3. Update YAML loader in `yaml_loader.jl` for new features
4. Add tests and examples
5. Update this README

## License

This extension is part of HydroModels.jl and uses the same license.
