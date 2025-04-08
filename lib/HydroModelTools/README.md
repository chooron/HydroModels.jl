# HydroModelTools.jl

[![英文文档](https://img.shields.io/badge/docs-dev-blue.svg)](https://chooron.github.io/HydroModelTools.jl/dev/)

The HydroModelTools.jl package provides a suite of tools for working with hydrological models, including ODE problem solvers and parameters optimizers. The package is built on the SciML ecosystem and serves as the support tools for the [HydroModels.jl](https://github.com/chooron/HydroModels.jl) package.

**Note**: If you don't use HydroModels.jl, we suggest using DifferentialEquations.jl and Optimization.jl directly.

## ODE Problem Solvers

The HydroModelTools.jl package provides a suite of ODE problem solvers which are wrappers of the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) package. For continuous and discrete problems, the package provides `ODESolver` and `DiscreteSolver` types, and provides the sensitivity analysis using [SciMLSensitivity.jl](https://github.com/SciML/SciMLSensitivity.jl).

## Parameter Optimizers

The HydroModelTools.jl package provides a suite of parameter optimizers which are wrappers of the [Optimization.jl](https://github.com/JuliaOpt/Optimization.jl) package. For black box optimization, the package provides `BatchOptimizer` type, and for gradient optimization, the package provides `GradOptimizer` type.

## Installation

```julia
using Pkg
Pkg.add("HydroModelTools")
```

## Contributing

HydroModelTools.jl is an open-source project available on Github. We welcome contributions from the community to help advance deep learning applications in hydrology.