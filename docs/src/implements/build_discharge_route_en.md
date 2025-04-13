# Building a Discharge Routing Model

## Mathematical Formulation

The discharge routing model is based on the following equations:

```math
\begin{aligned}
\frac{dS_{rf,n}}{dt} &= Q_{rf,up} - Q_{rf,n} && (1) \\
Q_{rf,n} &= \frac{S_{rf,n}}{\text{LAG}_{rf}} + Q_{rf,n-1} && (2) \\
Q_{rf,up} &= \sum_{i \in up(n)} Q_{rf,i} && (3) \\
Q_{rf} &= Q_{rf,n} + (R_{sw} + R_{gw}) \cdot A_{gc} && (4)
\end{aligned}
```

where $n$ is the reach index, ranging from ${1,...,c_{max}}$; $\Delta S_{rf,n}$ is the storage change in reach n; $Q_{rf,n}$ is the outflow from reach n; $S_{rf,n}$ is the storage in reach n; $\text{LAG}_{rf}$ is the routing lag parameter; $Q_{rf,up}$ is the inflow from upstream reaches; $R_{sw}$ and $R_{gw}$ are surface and groundwater runoff; $A_{gc}$ is the grid cell area; $c_{max}$ is the maximum number of reaches.

The equations describe: (1) storage change in the reach, (2) reach outflow calculation, (3) upstream inflow accumulation, and (4) final discharge calculation. This approach determines the inflow to each Hydrological Response Unit (HRU) based on its upstream inputs, which can be represented using two different routing functions. The outflow from the current HRU is calculated using the current reach storage $S_{rf,n}$ and combined with the HRU's $R_{sw}$ and $R_{gw}$ (after area conversion) to obtain the final discharge $Q_{rf}$.

## Implementation in HydroModels.jl

### Modified Calculation Formula

In the HydroModels.jl framework, this storage-based discharge calculation method is represented by the HydroRoute type, which solves the problem continuously through ODE equations. The implementation differs from the original formula in that it assumes no instantaneous changes when calculating $Q_{rf,n}$, transforming equation (2) into:

```math
\begin{aligned}
Q_{rf,n} &= S_{rf,n} \cdot \frac{1}{\text{LAG}_{rf} + 1} && (5)
\end{aligned}
```

The $Q_{rf,n-1}$ term is removed because in continuous change, $S_{rf,n}$ is constantly evolving, making it unnecessary to consider $Q_{rf,n-1}$ in calculating $Q_{rf,n}$. Additionally, there are dimensional differences between $S_{rf,n}$ and $Q_{rf,n}$, and including $Q_{rf,n}$ would add an instantaneous state variable, which is unsuitable as a state variable. In HydroRoute construction, the primary focus is expressing the outflow calculation formula. Like HydroBucket, HydroRoute accepts HydroFlux to represent the $Q_{rf,n}$ calculation formula:

```julia
using HydroModels
# Define parameters and variables
@variables q q_routed s_river
@parameters lag
# Define routing flux
rflux = HydroFlux([q, s_river] => [q_routed], [lag], exprs=[s_river / (1 + lag) + q])
```

### Flow Accumulation Functions

The method of obtaining upstream flow inputs is another crucial component of the routing model. Flow accumulation functions are typically used to calculate the upstream inflow for each HRU. This accumulation can be represented using either Grid-Based or Vector-Based approaches, corresponding to two types of HydroRoute construction in HydroModels.jl: GridRoute and VectorRoute.

#### Vector-Based Routing Model

The Vector-Based routing model is based on watershed subdivision and topological relationships, using an adjacency matrix for flow accumulation:

```julia
using Graphs
# Build river network topology
network = DiGraph(9)
add_edge!(network, 1, 2)
add_edge!(network, 2, 5)
add_edge!(network, 3, 5)
add_edge!(network, 4, 5)
add_edge!(network, 5, 8)
add_edge!(network, 6, 9)
add_edge!(network, 7, 8)
add_edge!(network, 8, 9)
vroute = HydroModels.VectorRoute(rfunc=rflux, rstate=s_river, network=network)
```

The code uses Graphs.jl's DiGraph data structure to represent watershed topology, which, along with the constructed rflux and s_river, completes the Vector-Based routing model construction.

#### Grid-Based Routing Model

The Grid-Based routing model represents HRU connectivity using a D8 flow direction matrix. The construction is based on the D8 matrix and corresponding HRU coordinates:

```julia
# Build river network topology
flwdir = [1 4 8; 1 4 4; 1 1 2]
positions = [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]
groute = HydroModels.GridRoute(rfunc=rflux, rstate=s_river, flwdir=flwdir, positions=positions)
```

The code uses matrix data to represent flow direction and HRU coordinates, which, combined with the constructed rflux and s_river, completes the Grid-Based routing model construction.
