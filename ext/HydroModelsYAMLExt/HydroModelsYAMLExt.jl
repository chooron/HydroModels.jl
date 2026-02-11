"""
    HydroModelsYAMLExt

Extension module for loading HydroModels from YAML configuration files.

This extension is automatically loaded when both HydroModels and YAML packages are available.

# Main Functions
- `load_model_from_yaml(yaml_file::String)`: Load a complete model from YAML file
- `load_config_from_yaml(yaml_file::String)`: Load configuration from YAML file
- `load_parameters_from_yaml(yaml_file::String)`: Load parameter metadata from YAML file

# Example
```julia
using HydroModels: toparam
using YAML  # This triggers the extension

# Load model from YAML
model = load_model_from_yaml("exphydro.yaml")

# Run model
output = model(forcing, params, config)
```

# YAML Format
The YAML file should contain:
- `parameters`: Parameter definitions with metadata
- `components`: List of component definitions (HydroBucket, etc.)
- `model`: Model composition specification
- `config`: Optional configuration settings

See the examples directory for sample YAML files.
"""
module HydroModelsYAMLExt

using HydroModels
using YAML
using Symbolics
using Symbolics: Num
using ComponentArrays: ComponentVector

# Include submodules
include("expression_parser.jl")
include("component_builder.jl")
include("yaml_loader.jl")
include("executor.jl")

end # module HydroModelsYAMLExt
