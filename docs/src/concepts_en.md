# HydroModels.jl Concepts

HydroModels.jl implements a modular and extensible architecture, designed to support diverse hydrological modeling paradigms. The framework's core structure comprises four main classes: Flux, Bucket, Route, Model and Wrapper, each playing a vital role in constructing comprehensive hydrological systems. The framework is implemented in Julia programming language for model development, computation, and parameter optimization. It leverages Lux.jl as the deep learning framework and integrates with SciML for scientific computing, symbolic programming, and parameter optimization capabilities. The conceptual framework of the model design philosophy is illustrated in the following diagram:

![Framework Concept](assets/concept.jpg)
Architecture and ecosystem of the HydroModels.jl framework. (a) illustrates the supporting ecosystem and its functional capabilities, highlighting the key dependencies and their roles in the framework; (b) shows the core architectural design of HydroModels.jl, demonstrating the main components and their interactions.

## Flux class: Water Process Transfer Representation

The Flux class serves as the fundamental building block of HydroModels.jl, drawing inspiration from flux library of the MARRMoT while extending its capabilities. This class encapsulates the physical equations governing water movement throughout the hydrological cycle. The conceptual formulation of the Flux class can be expressed as:


```math
\begin{aligned}
Y(t) = f(X(t),\theta) && (1) \\
\end{aligned}
```

where the functional representation $f$ maps input variables $X(t)$ and parameters $\theta$ to output variables $Y(t)$. Here, $X(t)$ represents input variables like precipitation and temperature, $\theta$ denotes the parameters that can be calibrated, and $Y(t)$ represents output variables like runoff and infiltration. Based on different applications, the Flux class includes:

- HydroFlux: Implements fundamental hydrological processes through mathematical formulations
- StateFlux: Handles mass-balance equations for state variables like soil moisture and snowpack
- NeuralFlux: Uses neural networks to model processes or predict parameters, leveraging Lux.jl capabilities

## Bucket Class: Storage Volume Change Simulation

The Bucket class represents fundamental water storage components within the hydrological system. It consists of multiple HydroFlux components (including NeuralFlux) and StateFlux components, which collectively define the water balance dynamics through coupled differential equations. These components can represent various hydrological stores such as soil moisture, groundwater reservoirs, and surface water bodies.

The Bucket class generates two essential functions:
1. An ODE solver function for state variables
2. A hydrological flux computation function for output calculations

This dual-function architecture enables efficient computation of hydrological processes, following a two-step process of solving ODEs and computing output fluxes.

## Route Class: Spatial Flow Propagation

The Route class simulates lateral water fluxes across landscapes and river networks, supporting both lumped "integral" and distributed "differential" models. Its core components include:

```math
\begin{aligned}
Q_{out}(t) &= f_{rflux}(S_{route}(t),Q_{gen}(t); ps) && (2) \\
\frac{dS_{route}}{dt} &= Q_{in}(t) - Q_{out}(t) && (3) \\
Q_{in}(t+1) &= f_{aggr}(Q_{out}(t)) && (4) \\
\end{aligned}
```

Where: $Q_{in}$ is inflow, $Q_{out}$ is outflow, $S_{route}$ is the routing state variable, and $f_{aggr}$ is the aggregation function.

The class supports various routing methods including:
- Standard state-update approaches
- Unit hydrograph through convolution operations
- Dynamic routing methods (Muskingum, kinematic wave)

## Model Class: Hydrological Process Integration

The Model class serves as the central management component, orchestrating the integration of Flux, Bucket, and Route components. It facilitates:
- Construction of distributed hydrological models
- Integration of vertical water movement with lateral connectivity
- Systematic parameter optimization and sensitivity analysis
- Efficient computation through metadata-driven approaches

## Configuration System

HydroModels.jl v0.5 introduces a type-stable configuration system using the `HydroConfig` struct:

### Key Components

```julia
struct HydroConfig{S<:SolverType, I, D}
    solver::S                    # Solver type
    interpolator::I              # Interpolator (wrapped in Val)
    timeidx::Vector{Int}         # Time indices
    device::D                    # Device function
    min_value::Float64           # Minimum threshold
    parallel::Bool               # Parallel flag
end
```

### Solver Types

- **MutableSolver**: Iterative updates, memory-efficient (default)
- **ImmutableSolver**: Functional accumulate, best for Zygote AD
- **ODESolver**: For DifferentialEquations.jl integration
- **DiscreteSolver**: For algebraic equations only

### Usage Example

```julia
# Create configuration
config = HydroConfig(
    solver = ImmutableSolver,
    interpolator = Val(DirectInterpolation),
    timeidx = 1:1000,
    min_value = 1e-6
)

# Run model
output = model(input, params, config; initstates = states)
```

### Performance Benefits

- **Type stability**: 100% type-stable for compiler optimization
- **Validation**: Automatic parameter checking (e.g., min_value > 0)
- **Flexibility**: Easy configuration merging and modification
- **Backward compatible**: NamedTuple configs still work

## Wrapper Class: Enhanced Component Capabilities

The Wrapper class extends component customization capabilities while maintaining interface uniformity. Key features include:
- NamedTupleIO wrapper for customized input/output specifications
- EstimateParam wrapper for parameter prediction based on basin characteristics
- RecordComponentState wrapper for state variable storage, supporting both batch training and online forecasting