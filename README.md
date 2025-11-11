# HydroModels.jl

[![version](https://docs.juliahub.com/HydroModels/version.svg)](https://juliahub.com/ui/Packages/General/HydroModels)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://chooron.github.io/HydroModels.jl/dev/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

> A modern hydrological modeling framework for flexible model construction, deep learning integration, and high-performance computation

## 📖 Overview

**HydroModels.jl** is a comprehensive hydrological modeling framework that extends and enhances SUPERFLEX's design philosophy. Built on the Julia language and the SciML (Scientific Machine Learning) ecosystem, it combines flexible model construction with computational efficiency, particularly supporting deep learning integration in hydrological modeling.

The framework is built on [HydroModelCore.jl](https://github.com/chooron/HydroModelCore.jl), which provides the foundational type system, runtime code generation, and automatic differentiation support.

## ✨ Key Features

- **🧩 Flexible Model Construction**: Support for lumped, semi-distributed, and distributed hydrological models
- **🤖 Deep Learning Integration**: Seamless neural network integration for enhanced flux calculations and dynamic parameter estimation
- **⚡ High Performance**: Leverages Julia's performance and the SciML ecosystem for efficient computation
- **🔍 Gradient-Based Optimization**: Full support for automatic differentiation (Zygote, ForwardDiff) and advanced optimization
- **🌊 Modular Components**: Compose models from reusable buckets, fluxes, and routing components
- **🔬 Scientific Ecosystem**: Built on ComponentArrays, Symbolics, Lux, and OrdinaryDiffEq

## 🚀 Installation

```julia
using Pkg
Pkg.add("HydroModels")
```

Or install from GitHub for the latest development version:

```julia
Pkg.add(url="https://github.com/chooron/HydroModels.jl")
```

## 📚 Quick Start

### Example 1: Building a Simple Bucket Model

```julia
using HydroModels
using ComponentArrays

# Define a helper function
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

# Declare variables and parameters
@variables temp lday prcp pet snowfall rainfall snowpack melt
@parameters Tmin Tmax Df

# Define a snow bucket with fluxes and state
snow_bucket = @hydrobucket :snow begin
    fluxes = begin
        @hydroflux begin
            snowfall ~ step_func(Tmin - temp) * prcp
            rainfall ~ step_func(temp - Tmin) * prcp
        end
        @hydroflux melt ~ step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))
        @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
    end
    dfluxes = begin
        @stateflux snowpack ~ snowfall - melt
    end
end

# Setup parameters and initial states
params = ComponentVector(params = ComponentVector(Tmin = -2.0, Tmax = 0.0, Df = 2.5))
initstates = ComponentVector(snowpack = 0.0)

# Prepare input data (example with 100 time steps)
input_data = (
    temp = rand(100) .* 20 .- 5,    # Temperature: -5 to 15°C
    lday = fill(0.5, 100),           # Day length: 12 hours
    prcp = rand(100) .* 10           # Precipitation: 0-10 mm/day
)
input_mat = reduce(hcat, [input_data.temp, input_data.lday, input_data.prcp])'

# Run the model
config = (solver = HydroModels.ManualSolver(), timeidx = 1:100)
results = snow_bucket(input_mat, params, config; initstates = initstates)

# Results contain: [snowpack, snowfall, rainfall, melt, pet] × timesteps
println("Result dimensions: ", size(results))  # (5, 100)
```

### Example 2: Building a Complete Hydrological Model (ExpHydro)

```julia
using HydroModels
using ComponentArrays

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

# Define all variables and parameters
@variables temp lday pet prcp snowfall rainfall snowpack melt
@variables soilwater evap baseflow surfaceflow flow
@parameters Tmin Tmax Df Smax Qmax f

# Snow and surface processes bucket
snow_bucket = @hydrobucket :surface begin
    fluxes = begin
        @hydroflux begin
            snowfall ~ step_func(Tmin - temp) * prcp
            rainfall ~ step_func(temp - Tmin) * prcp
        end
        @hydroflux melt ~ step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))
        @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
    end
    dfluxes = begin
        @stateflux snowpack ~ snowfall - melt
    end
end

# Soil water and runoff generation bucket
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

# Compose the complete model
exphydro_model = @hydromodel :exphydro begin
    snow_bucket
    soil_bucket
end

# Setup parameters
params = ComponentVector(
    params = ComponentVector(
        Tmin = -2.0, Tmax = 0.0, Df = 2.5,
        Smax = 1500.0, Qmax = 20.0, f = 0.01
    )
)
initstates = ComponentVector(snowpack = 0.0, soilwater = 1000.0)

# Prepare input data
input_data = (
    temp = rand(365) .* 20 .- 5,
    lday = fill(0.5, 365),
    prcp = rand(365) .* 10
)
input_mat = reduce(hcat, [input_data.temp, input_data.lday, input_data.prcp])'

# Run the complete model
config = (solver = HydroModels.ManualSolver(), timeidx = 1:365)
results = exphydro_model(input_mat, params, config; initstates = initstates)

# Extract outputs
flow_series = results[end, :]  # Last row contains flow output
println("Mean annual flow: ", mean(flow_series), " mm/day")
```

### Example 3: Gradient-Based Parameter Optimization

```julia
using HydroModels
using ComponentArrays
using Zygote

# Define a simple model (using exphydro_model from Example 2)
# ... (model definition code here)

# Define loss function
function loss_fn(params, input_mat, observed_flow, initstates)
    config = (solver = HydroModels.ManualSolver(mutable=false), timeidx = 1:length(observed_flow))
    predictions = exphydro_model(input_mat, params, config; initstates = initstates)
    predicted_flow = predictions[end, :]
    return sum((predicted_flow .- observed_flow).^2)  # MSE loss
end

# Compute gradients using Zygote
observed_flow = rand(365) .* 5  # Synthetic observations
gradients = Zygote.gradient(params) do p
    loss_fn(p, input_mat, observed_flow, initstates)
end

println("Parameter gradients computed successfully!")
println("Gradient for Smax: ", gradients[1].params.Smax)
```

## 🏗️ Framework Architecture

HydroModels.jl provides several key component types:

- **`@hydroflux`**: Define individual flux calculations
- **`@stateflux`**: Define state differential equations
- **`@hydrobucket`**: Compose fluxes and states into storage components
- **`@hydromodel`**: Assemble buckets into complete models
- **`@unithydro`**: Add unit hydrograph routing
- **`@hydroroute`**: Add spatial routing (grid-based or vector-based)

### Neural Network Integration

```julia
using Lux

# Define a neural network-enhanced flux
nn_flux = @neuralflux begin
    nn_model = Lux.Chain(
        Dense(3, 10, tanh),
        Dense(10, 1)
    )
    input_names = [:temp, :prcp, :soilwater]
    output_name = :nn_et
end
```

## 📖 Documentation

- **[Online Documentation](https://chooron.github.io/HydroModels.jl/dev/)**: Comprehensive guides and tutorials
- **[API Reference](https://chooron.github.io/HydroModels.jl/dev/)**: Detailed function documentation
- **[Examples](test/)**: Working examples and test cases
- **[HydroModelCore.jl](https://github.com/chooron/HydroModelCore.jl)**: Core framework documentation

## 🌊 Model Library

HydroModels.jl includes implementations of several classic hydrological models:

- **ExpHydro**: Simple lumped conceptual model with snow module
- **GR4J**: Four-parameter daily model with unit hydrographs
- **HBV**: Scandinavian conceptual model
- **M50**: Neural network-enhanced hybrid model
- **XAJ**: Xinanjiang model (Chinese watershed model)

See the [test/models/](test/models/) directory for complete implementations.

## 🔗 Ecosystem Integration

HydroModels.jl integrates seamlessly with the Julia ecosystem:

- **[SciML](https://sciml.ai/)**: ODE solving, sensitivity analysis, optimization
- **[Lux.jl](https://github.com/LuxDL/Lux.jl)**: Neural network integration
- **[ComponentArrays.jl](https://github.com/jonniedie/ComponentArrays.jl)**: Named parameter arrays
- **[Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl)**: Symbolic computation
- **[Optimization.jl](https://github.com/SciML/Optimization.jl)**: Advanced optimization algorithms

## 📊 Performance Features

- **Runtime Code Generation**: Models compiled to efficient Julia code
- **Automatic Differentiation**: Full gradient support for optimization
- **Multiple Solvers**: Manual solvers (mutable/immutable) and ODE solvers
- **Distributed Computing**: Support for multi-node spatial models
- **GPU Support**: CUDA-compatible operations (experimental)

## ⚠️ Known Issues & Future Plans

- When solving ODE problems, `EnzymeVJP` from SciMLSensitivity.jl provides efficient gradients but ~~has issues with multi-node problems~~ (being addressed)
- Future plans include:
  - Enhanced gradient computation using [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl)
  - Code compilation with [Reactant.jl](https://github.com/EnzymeAD/Reactant.jl)
  - Extended support for large-scale model generation

## 🤝 Contributing

HydroModels.jl is an open-source project. We welcome contributions from the community!

- **Report Issues**: [GitHub Issues](https://github.com/chooron/HydroModels.jl/issues)
- **Submit PRs**: Help improve the code and documentation
- **Share Models**: Contribute new model implementations

## 📝 Citation

If you use HydroModels.jl in your research, please cite:

```bibtex
@software{hydromodels_jl,
  author = {Jing Xu},
  title = {HydroModels.jl: A Modern Hydrological Modeling Framework},
  year = {2024},
  url = {https://github.com/chooron/HydroModels.jl}
}
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 👤 Author

**Jing Xu** ([@chooron](https://github.com/chooron))

*Note: I am not a professional full-time software developer. If you find any issues or have suggestions, please feel free to post them in the [issues](https://github.com/chooron/HydroModels.jl/issues) section.*
