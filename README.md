# HydroModels.jl

[![version](https://docs.juliahub.com/HydroModels/version.svg)](https://juliahub.com/ui/Packages/General/HydroModels)
[![英文文档](https://img.shields.io/badge/docs-dev-blue.svg)](https://chooron.github.io/HydroModels.jl/dev/)
[![中文文档](https://img.shields.io/badge/docs-zh-red.svg)](https://chooron.github.io/HydroModels.jl/dev-zh/)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

## Overview

HydroModels.jl is a modern hydrological modeling framework that extends and enhances SUPERFLEX's design philosophy. Built on Julia language and the SciML (Scientific Machine Learning) ecosystem, it combines flexible model construction with computational efficiency, particularly supporting deep learning integration in hydrological modeling.

## Key Features

- **Flexible Model Construction**: Supports development of lumped, semi-distributed, and distributed hydrological models
- **Deep Learning Integration**: Enables neural network integration for enhanced flux calculations and dynamic parameter estimation
- **Computational Efficiency**: Leverages Julia's high-performance capabilities and the SciML ecosystem
- **Gradient-Based Optimization**: Supports advanced parameter optimization techniques by Optimization.jl
- **Comprehensive Framework**: Provides tools for both traditional hydrological modeling and modern machine learning approaches

## Framework Capabilities

HydroModels.jl offers:
- Easy implementation and customization of hydrological models
- Integration with Julia's scientific computing ecosystem
- Support for various modeling approaches:
  - Traditional conceptual models
  - Neural network enhanced models
  - Distributed hydrological systems

## Installation

```julia
using Pkg
Pkg.add("HydroModels")
```

## Important Notations

- I am not a professional full-time software developer, so the code may have imperfections. If you find any issues, please post them in the [issues](https://github.com/chooron/HydroModels.jl/issues) section, and I will respond promptly.
- The documentation is still being improved. If you want to quickly understand how to use the software, you can refer to the [dev](https://github.com/chooron/HydroModels.jl/tree/main/dev) folder on GitHub. I will integrate this into a notebook as soon as possible.
- The software has some issues. When solving ODE files, using EnzymeVJP from SciMLSensitivity.jl generally provides more efficient gradient calculations, but there are problems with multi-node ODE solving. I will fix this as soon as possible and begin trying to calculate gradients using [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl) and compile code using [Reactant.jl](https://github.com/EnzymeAD/Reactant.jl).
- To be responsible, if you only want to build a hydrological model or a neural network-embedded hydrological model, I would recommend building it manually. My main purpose in writing this package is to reduce the difficulty of model development and to serve future large model generation tasks.

## Contributing

HydroModels.jl is an open-source project available on Github. We welcome contributions from the community to help advance deep learning applications in hydrology.