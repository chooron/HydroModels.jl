# Multi-Node/Distributed Hydrological Modeling

## Introduction to Multi-Node/Distributed Hydrological Modeling

The hydrological models introduced earlier are primarily lumped models, which use a single model to simulate an entire watershed. This approach is typically suitable for smaller watersheds. For large watersheds with significant topographic features and river network topology, a single hydrological model is insufficient to simulate the processes accurately.

Researchers have proposed methods to build hydrological models based on watershed topographic characteristics and then integrate them using routing models. This approach is known as distributed hydrological modeling. The process begins by dividing the watershed into multiple hydrological response units (which we call multi-nodes). This division can be done either by partitioning the watershed into standard grids of fixed size (distributed approach) or by dividing it into irregularly shaped sub-watersheds based on river topology (semi-distributed approach). The runoff results from these multiple nodes are then integrated using geographic elevation information and river network topology.

## Implementation of Multi-Node Models

All types constructed in HydroModels.jl (Flux, Bucket, Route, Model) support multi-node computation without requiring additional declarations during construction. Each module automatically allocates the appropriate computational scenario based on the input data type. This section primarily introduces the construction methods for routing modules. Using the HBV model as the runoff generation model, we have constructed three types of routing models: unit hydrograph-based routing, semi-distributed routing based on river network topology, and distributed routing based on geographic elevation information.

### 1. Unit Hydrograph-Based Routing Model

HBV-maxbas uses HBV as its foundation and a unit hydrograph model for routing. The unit hydrograph parameters are determined by the distance from the computational grid to the outlet node. The implementation is as follows:

```julia
# Define unit hydrograph using macro
uh = @unithydro :maxbas_uh begin
    uh_func = begin
        2lag => (1 - 0.5 * (2 - t / lag)^2.5)
        lag => (0.5 * (t / lag)^2.5)
    end
    uh_vars = q => q_uh  # Input => Output
end
```

Both semi-distributed routing based on river network topology and distributed routing based on geographic elevation information can be constructed using HydroRoute. The difference between these two construction methods lies in the aggregate function.

Semi-distributed hydrological models typically describe the connectivity of sub-watersheds using directed graphs. This allows for the integration of runoff results from different sub-watersheds based on the connection relationships provided by the network:

```julia
# Define routing component with aggregation function
discharge_route = @hydroroute :exphydro_routed begin
    fluxes = begin
        @hydroflux q_routed ~ s_river / (1 + lag) + q
    end
    dfluxes = begin
        @stateflux s_river ~ q - q_routed
    end
    hru_types = collect(1:9)  # 9 nodes
    aggr_func = HydroModels.build_aggr_func(network)
end
```

The construction of `fluxes` and `dfluxes` is similar to that of Bucket, but with some limitations:

1. The number of `fluxes` and `dfluxes` must correspond. `fluxes` represents the flow from the current computational unit to the next, while `dfluxes` represents the storage volume of the routing module in the current computational unit. It is composed of the runoff from the current unit, inflow from other units, and outflow to other units. The complete formula should be:

$$
\frac{dS_{route}}{dt} = Q_{in}(t) - Q_{out}(t) + Q_{gen}(t)
$$

Here, we assume that $Q_{in}(t)$ equals $aggr\_func(Q_{out}(t))$, meaning the inflow for the current time period equals the result of routing the outflow for the current time period (this process is continuous, so we consider this expression valid). To simplify the construction of the Route module, this calculation is relatively fixed, but more modifications can be made in the writing of `fluxes`.

1. Compared to Bucket, the Route module also needs to provide an aggregate function. This function is constructed from the flow direction matrix calculated from the topological relationships of sub-watersheds or geographic elevation data, and is used to calculate $Q_{in}(t)$ for the current time period.

For semi-distributed hydrological models, a directed graph can be constructed based on the upstream-downstream connection relationships of the watershed (using Graphs.jl):

```julia
network = DiGraph(9)
add_edge!(network, 1, 2)
add_edge!(network, 2, 5)
add_edge!(network, 3, 5)
add_edge!(network, 4, 5)
add_edge!(network, 5, 8)
add_edge!(network, 6, 9)
add_edge!(network, 7, 8)
add_edge!(network, 8, 9)

aggr_func = HydroModels.build_aggr_func(network)
```

The code represents the connection relationships of the watershed by constructing a directed graph, where the numbers in the directed graph represent the node numbers.

For distributed hydrological models, a flow direction matrix (d8) and a list of matrix indices representing the input node index order are needed to ensure the correctness of the input integration:

```julia
flwdir = [1 4 8; 1 4 4; 1 1 2]
positions = [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]
```

## Implementing Multi-Node Computation with HydroModels.jl

All types constructed in HydroModels.jl (Flux, Bucket, Route, Model) support multi-node computation without requiring additional declarations during construction. Each module automatically allocates the appropriate computational scenario based on the input data type. Below, I will introduce the execution process of multi-node computation from three aspects: data preparation, parameter settings, and runtime configuration.

### 1. Data Preparation

For single-node input, the data format is variables × time dimensions. For multi-node input, the data format becomes variables × nodes × time dimensions, resulting in a three-dimensional matrix.

### 2. Parameter Preparation

For single-node input, the parameter format is:

```julia
params = ComponentVector(
    params=(p1=0.0, p2=0.0, ...),
    nns=(nn1=[...], )
)
```

As you can see, each parameter name corresponds to a single value.

For multi-node input, the parameter format is:

```julia
params = ComponentVector(
    params=(p1=[0.0, 0.0, ...], p2=[0.0, 0.0, ...], ...),
    nns=(nn1=[...], )
)
```

Each parameter corresponds to an array of values, representing the parameters for each node. Note that the neural network parameters do not show significant changes. This is because constructing a neural network for each node during multi-node computation would make the computational cost relatively large and unnecessary.

Similarly, the initial states follow the same pattern as the parameters:

```julia
initstates = ComponentVector(
    s1=[0.0, 0.0, ...],
    s2=[0.0, 0.0, ...],
    ...
)
```

### 3. Runtime Configuration

#### HydroConfig for Multi-Node Models (NEW in v2.0)

For multi-node models, you use the same `HydroConfig` system as single-node models:

```julia
# Create configuration for multi-node model
config = HydroConfig(
    solver = MutableSolver,          # or ImmutableSolver for AD
    interpolator = Val(DirectInterpolation),
    timeidx = 1:365,                 # Time steps
    min_value = 1e-6,
    parallel = false                  # Future: enable parallel computation
)
```

#### Parameter Type Index (hru_types)

In distributed hydrological models, some computational units have highly similar natural attributes, so it's reasonable to share parameters. This is achieved using `hru_types` during bucket/route construction:

```julia
# Define parameter types for each node
# Suppose we have 10 nodes falling into 4 terrain types
hru_types = [1, 1, 1, 2, 2, 2, 3, 3, 3, 4]

# Define bucket with hru_types
snow_bucket = @hydrobucket :snow begin
    fluxes = begin
        @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
        # ... other fluxes
    end
    dfluxes = begin
        @stateflux snowpack ~ snowfall - melt
    end
    hru_types = [1, 1, 1, 2, 2, 2, 3, 3, 3, 4]  # Specify node types
end
```

In the multi-node input parameters, each parameter corresponds to a parameter array with one value per terrain type:

```julia
# Parameters for 4 terrain types
params = ComponentVector(
    params = (
        Df = [2.5, 3.0, 2.0, 2.8],      # 4 values for 4 types
        Tmax = [0.5, 0.3, 0.7, 0.4],
        # ... other parameters
    )
)
```

The framework automatically expands parameters based on `hru_types` during computation:
```julia
# Internally expands: Df[hru_types] = [2.5, 2.5, 2.5, 3.0, 3.0, 3.0, 2.0, 2.0, 2.0, 2.8]
```

This significantly reduces the number of parameters to calibrate while maintaining spatial heterogeneity.

## Real-World Example (Based on HBV+maxbas/route Model)

To demonstrate the practical application of multi-node modeling, let's consider a real-world example using the HBV model for runoff generation combined with either maxbas or route for flow routing.

First, we'll define our watershed structure using a directed graph to represent the connectivity between sub-watersheds:

```julia
# Define watershed network structure
network = DiGraph(9)  # 9 sub-watersheds
add_edge!(network, 1, 2)  # Sub-watershed 1 flows into 2
add_edge!(network, 2, 5)
add_edge!(network, 3, 5)
add_edge!(network, 4, 5)
add_edge!(network, 5, 8)
add_edge!(network, 6, 9)
add_edge!(network, 7, 8)
add_edge!(network, 8, 9)  # Sub-watershed 8 flows into outlet 9

# Create the routing component
discharge_route = @hydroroute :hbv_routed begin
    fluxes = begin
        @hydroflux q_routed ~ s_river / (1 + lag) + q
    end
    dfluxes = begin
        @stateflux s_river ~ q - q_routed
    end
    aggr_func = HydroModels.build_aggr_func(network)
end

# Combine with HBV model
hbv_distributed = @hydromodel :hbv_distributed begin
    hbv_bucket  # Runoff generation component
    discharge_route  # Routing component
end
```

Next, we prepare our input data and parameters for the multi-node model:

```julia
# Prepare input data (variables × nodes × time)
input_data = zeros(3, 9, 365)  # Example: 3 variables, 9 nodes, 365 days
# Fill with actual precipitation, temperature, etc. data for each node

# Define parameter types (4 types of terrain)
params = ComponentVector(
    params = (
        FC = [150.0, 200.0, 100.0, 180.0],    # Field capacity for each terrain type
        LP = [0.7, 0.8, 0.6, 0.75],
        BETA = [2.0, 3.0, 1.5, 2.5],
        # Other HBV parameters...
        lag = [1.5, 2.0, 1.0, 1.8]             # Routing lag parameter
    )
)

# Initial states for each node (9 nodes total)
initstates = ComponentVector(
    snowpack = zeros(9),
    soilwater = zeros(9),
    s_river = zeros(9)
)

# Configure model execution (NEW in v2.0)
config = HydroConfig(
    solver = MutableSolver,
    interpolator = Val(DirectInterpolation),
    timeidx = 1:365,
    min_value = 1e-6
)

# Run the distributed model
results = hbv_distributed(
    input_data,
    params,
    config;
    initstates = initstates
)
```

> **Note**: The `hru_types` parameter is now specified during model construction (in the bucket definition), not at runtime. This allows for better type stability and compiler optimization.

This example demonstrates how HydroModels.jl can efficiently handle multi-node hydrological modeling with parameter sharing across similar terrain types, significantly reducing the number of parameters that need to be calibrated while maintaining the spatial heterogeneity of the watershed.
