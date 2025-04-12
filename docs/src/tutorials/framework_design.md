# Understanding DeepFlex.jl Core Concepts

This document explains the core concepts of the DeepFlex.jl framework, focusing on the three fundamental components: Flux, Element/Bucket, and Model. We'll use examples from the ExpHydro model to illustrate these concepts in practice.

## Framework Architecture

DeepFlex.jl is built on a modular architecture with three core components:

1. **Flux**: Basic computational units that define mathematical relationships
2. **Element/Bucket**: Storage components that integrate fluxes and manage state variables
3. **Model**: The top-level structure that combines all components into a complete system

This architecture allows you to build models ranging from simple conceptual models to complex physically-based or hybrid models using the same framework.

## Flux Components

Fluxes are the fundamental building blocks in DeepFlex.jl. They define the mathematical relationships between variables in your hydrological system.

### HydroFlux

A `HydroFlux` defines a mathematical relationship between input and output variables using symbolic expressions.

In the ExpHydro model, several HydroFlux components are used to calculate processes like snowfall, rainfall, and potential evapotranspiration:

```julia
@hydroflux snowfall ~ step_func(Tmin - temp) * prcp
@hydroflux rainfall ~ step_func(temp - Tmin) * prcp
@hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
```

These expressions define:

- How precipitation is partitioned into snowfall and rainfall based on temperature
- How potential evapotranspiration is calculated from temperature and day length

HydroFlux components are used for processes that can be expressed as explicit mathematical formulas.

### StateFlux

A `StateFlux` defines how state variables change over time. These are used to create differential equations that are solved during simulation.

In the ExpHydro model, StateFlux components define how the snowpack and soil water storage change:

```julia
@stateflux snowpack ~ snowfall - melt
@stateflux soilwater ~ (rainfall + melt) - (evap + flow)
```

These expressions define:

- The snowpack increases with snowfall and decreases with snowmelt
- The soil water increases with rainfall and snowmelt, and decreases with evapotranspiration and streamflow

StateFlux components are essential for water balance calculations and form the basis of the ordinary differential equations (ODEs) that are solved during simulation.

### NeuralFlux

A `NeuralFlux` integrates neural networks into the hydrological modeling framework. This allows for data-driven components to be seamlessly combined with process-based equations.

While the basic ExpHydro model doesn't use neural networks, more advanced models in DeepFlex.jl can incorporate neural networks for processes that are difficult to model with explicit equations:

```julia
# Example from M100 model
@neuralflux [log_evap_div_lday, log_flow, asinh_melt, asinh_ps, asinh_pr] ~ m100_nn([norm_snw, norm_slw, norm_prcp, norm_temp])
```

NeuralFlux components maintain the same interface as traditional fluxes, allowing for seamless integration of data-driven and process-based approaches.

## Element/Bucket Components

Elements are higher-level components that integrate multiple fluxes and manage state variables. The primary type of element in DeepFlex.jl is the `Bucket`.

### HydroBucket

A `HydroBucket` represents a storage component with state variables that change over time according to input and output fluxes. It encapsulates both water movement processes and state evolution equations.

In the ExpHydro model, two buckets are defined:

1. **Snow Bucket**: Handles snow accumulation and melt processes

```julia
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
```

2. **Soil Water Bucket**: Manages soil moisture, evapotranspiration, and runoff generation

```julia
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
```

Each bucket:

- Integrates multiple flux components
- Manages state variables and their dynamics
- Handles ODE solving for time evolution
- Supports both single-node and multi-nodes simulations

The bucket structure allows for modular model construction, where each bucket represents a distinct hydrological process or storage component.

### HydroRoute

Another type of element is `HydroRoute`, which handles water routing through a network of connected nodes. While not used in the basic ExpHydro model, routing components are essential for distributed hydrological modeling:

```julia
# Example routing component
route = @hydroroute :river_routing begin
    fluxes = begin
        @hydroflux outflow ~ inflow * (1 - exp(-k * travel_time))
    end
    dfluxes = begin
        @stateflux channel_storage ~ inflow - outflow
    end
    aggr_func = network_aggregation_function
end
```

Routing components manage flow calculations and inter-node water transfer in network-based models.

## Specialized Structures

Hydrological models incorporate several specialized structures, with the unit hydrograph being the most typical example. The unit hydrograph is a module used to describe hillslope routing, simulating the process of runoff concentration after runoff generation through convolution calculations.

```math
Q(t) = \sum_{i=1}^{t} P(t-i+1) \cdot UH(i)
```

To represent this concept, HydroModels.jl has designed two types: `UHFunction` and `UnitHydrograph`. These represent the unit hydrograph function and the unit hydrograph model, respectively.

The unit hydrograph function calculates the distribution weights for different time periods in the runoff process, while the unit hydrograph model describes which flux is used for unit hydrograph routing calculations. Here, we use the unit hydrographs for slow flow and fast flow routing calculations in the GR4J model as implementation examples.

```julia
# Unit hydrograph components for GR4J model
uh = @unithydro :maxbas_uh begin
    uh_func = begin
        2lag => (1 - 0.5 * (2 - t / lag)^2.5)
        lag => (0.5 * (t / lag)^2.5)
    end
    uh_vars = [q]
    configs = (solvetype=:DISCRETE, suffix=:_lag)
end
```

## Model: The Integration Framework

The `HydroModel` is the top-level structure that integrates multiple components into a complete simulation system.

In the ExpHydro model, the snow and soil buckets are combined into a complete model:

```julia
exphydro_model = @hydromodel :exphydro begin
    snow_bucket
    soil_bucket
end
```

The model automatically:

- Analyzes dependencies between components
- Creates a computational graph for efficient execution
- Manages state variables across the entire model
- Prepares the model for simulation with various solver options

The model provides a unified interface for simulation, allowing you to run the model with different input data, parameters, and configuration options:

```julia
# Run the model
results = exphydro_model(
    input_matrix,   # Input data matrix
    params,         # Parameters
    initstates = init_states,  # Initial states
    config = config  # Configuration
)
```

## How Components Work Together

The power of DeepFlex.jl comes from how these components work together:

1. **Fluxes** define the mathematical relationships and processes
2. **Buckets/Elements** integrate these fluxes and manage state variables
3. **Models** combine multiple elements into a complete simulation system

In the ExpHydro model:

- HydroFlux components define processes like snowfall, rainfall, and potential evapotranspiration
- StateFlux components define how snowpack and soil water storage change over time
- HydroBucket components integrate these fluxes into snow and soil water components
- HydroModel combines these components into a complete ExpHydro model

The framework automatically handles the connections between components. For example, the rainfall and melt outputs from the snow bucket are used as inputs to the soil bucket without requiring explicit connections.

## Execution Process

When running a model in DeepFlex.jl, the following steps occur:

1. **Data Preparation**:

   - Input data is organized as a matrix (features × time)
   - Parameters are provided as a ComponentVector
   - Initial states are specified
2. **Configuration**:

   - Solver selection (e.g., ManualSolver for explicit Euler method)
   - Interpolation method (e.g., LinearInterpolation)
   - Time indices for output
3. **Execution**:

   - Components are executed in dependency order
   - State variables are integrated over time
   - Fluxes are calculated at each time step
   - Results are collected and organized
4. **Output**:

   - Time series of state variables (e.g., snowpack, soil moisture)
   - Time series of output fluxes (e.g., streamflow)

## Conclusion

Understanding the relationships between Flux, Element/Bucket, and Model components is essential for effectively using the HydroModels.jl framework. These components provide a flexible and powerful system for building hydrological models of varying complexity.

The ExpHydro model demonstrates how these components work together to create a complete hydrological model:

- Flux components define the mathematical relationships
- Bucket components integrate these fluxes and manage state variables
- The Model combines these components into a complete simulation system

By leveraging this modular architecture, you can:

- Build models from reusable components
- Combine process-based and data-driven approaches

This framework design enables a wide range of hydrological modeling applications while maintaining computational efficiency and scientific rigor.
