# HydroModelsOptimizationExt

Extension for integrating HydroModels with the SciML Optimization ecosystem.

## Features

### 1. Direct Metric Optimization

Use hydrological metrics directly without defining custom loss functions:

```julia
prob = OptimizationProblem(
    component, input, target;
    metric="KGE",  # Supported: KGE, NSE, LogKGE, MSE
    lb_pas=[0.1, 0.1],
    ub_pas=[10.0, 10.0]
)
```

**Supported Metrics:**
- `KGE`: Kling-Gupta Efficiency (good for overall performance)
- `NSE`: Nash-Sutcliffe Efficiency (emphasizes high flows)
- `LogKGE`: Log-transformed KGE (better for low flows)
- `MSE`: Mean Squared Error (default)

### 2. Fixed Parameters with ComponentVector

Calibrate only a subset of parameters while keeping others fixed. Supports hierarchical parameter structures including both traditional parameters and neural network parameters:

```julia
# Fix traditional parameters
prob = OptimizationProblem(
    component, input, target;
    metric="NSE",
    fixed_params=ComponentVector(params=(k2=2.5,)),  # Fix k2, calibrate others
    lb_pas=[0.1, 0.1],  # Bounds for calibratable params only
    ub_pas=[10.0, 10.0]
)

# Fix both traditional and neural network parameters
# Neural network parameters are automatically initialized
prob = OptimizationProblem(
    model, input, target;
    fixed_params=ComponentVector(
        params=(k=0.5, b=2.0),  # Fix traditional parameters
        nns=(et_net=nn_params,)  # Fix neural network parameters
    ),
    lb_pas=[0.1],  # Bounds for remaining calibratable params
    ub_pas=[5.0]
)
```

This is useful for:
- Identifiability studies
- Sensitivity analysis
- Constraining poorly-identified parameters
- Hybrid modeling (fix neural networks, calibrate physics-based params)

### 3. Neural Network Support

Automatic initialization and management of neural network parameters:

```julia
using HydroModels
using Lux

# Create model with neural flux
@variables prcp temp et
nn = Chain(Dense(2 => 10, tanh), Dense(10 => 1), name=:et_net)
neural_flux = @neuralflux et ~ nn([prcp, temp])
model = HydroModel([neural_flux, other_components...])

# Neural network parameters are automatically initialized
# No need to manually extract nn_params!
prob = OptimizationProblem(
    model, input, target;
    fixed_params=ComponentVector(params=(k=0.5,))  # Fix traditional params
)

# Or get initial parameters explicitly if needed
initial_params = get_initial_params(model)
# Returns: ComponentVector(params=(k=1.0, ...), nns=(et_net=ComponentVector(...),))
```

**Supported Components:**
- `NeuralFlux`: Single neural network
- `NeuralBucket`: Three networks (flux, state, output)
- `HydroBucket`: Extracts from embedded neural fluxes
- `HydroModel`: Recursively collects from all sub-components

### 4. Custom Loss Functions (Backward Compatible)

You can still use custom loss functions:

```julia
custom_loss = (obs, sim) -> sum(abs.(obs .- sim))

prob = OptimizationProblem(
    component, input, target;
    loss_func=custom_loss,
    lb_pas=[0.1, 0.1],
    ub_pas=[10.0, 10.0]
)
```

### 5. Automatic Differentiation Support

```julia
using Optimization, OptimizationOptimJL

prob = OptimizationProblem(
    component, input, target;
    metric="KGE",
    adtype=Optimization.AutoForwardDiff(),  # Enable AD
    lb_pas=[0.1, 0.1],
    ub_pas=[10.0, 10.0]
)

sol = solve(prob, BFGS())  # Gradient-based optimizer
```

## Complete Examples

### Example 1: Basic Calibration

```julia
using HydroModels
using Optimization
using OptimizationBBO

# Setup
component = ExpHydro(name=:exphydro)
input = rand(5, 365)  # 5 inputs, 365 days
target = rand(365)    # Observed streamflow

# Create optimization problem
prob = OptimizationProblem(
    component,
    input,
    target;
    metric="KGE",
    warm_up=30,  # Skip first 30 days
    fixed_params=ComponentVector(params=(f=0.5,)),  # Fix one parameter
    lb_pas=[0.1, 0.1, 0.1, 0.1],
    ub_pas=[10.0, 10.0, 10.0, 10.0]
)

# Solve with BBO algorithm
sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=1000)

# Extract results
optimal_params = sol.u
final_loss = sol.objective
```

### Example 2: Hybrid Model with Neural Networks

```julia
using HydroModels
using Lux
using ComponentArrays

# Create hybrid model
@variables prcp temp et q
nn = Chain(Dense(2 => 10, tanh), Dense(10 => 1), name=:et_net)
neural_flux = @neuralflux et ~ nn([prcp, temp])
model = HydroModel([neural_flux, other_components...])

# Get initial parameters (optional, for inspection)
initial_params = get_initial_params(model)

# Get only neural network parameters if needed
nn_params = get_nn_initial_params(model)

# Calibrate traditional parameters, fix neural networks
prob = OptimizationProblem(
    model, input, target;
    fixed_params=ComponentVector(
        nns=nn_params  # Fix all neural network parameters
    ),
    lb_pas=[0.1, 0.1],
    ub_pas=[10.0, 10.0]
)

sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=1000)
```

## API Reference

### OptimizationProblem

```julia
OptimizationProblem(
    component::AbstractComponent,
    input::AbstractArray{T,2},
    target::AbstractVector;
    metric::Union{String,Nothing} = nothing,
    loss_func::Function = (y, y_hat) -> sum((y .- y_hat) .^ 2) / length(y),
    fixed_params::Union{ComponentVector,NamedTuple,Nothing} = nothing,
    initial_params::Union{ComponentVector,NamedTuple,Nothing} = nothing,
    warm_up::Int = 1,
    lb_pas::Union{AbstractVector,Nothing} = nothing,
    ub_pas::Union{AbstractVector,Nothing} = nothing,
    adtype = nothing,
    interpolator = Val(HydroModels.ConstantInterpolation),
    timeidx = collect(1:size(input, 2)),
    solver = HydroModels.MutableSolver,
    solve_alg = nothing,
    sense_alg = nothing,
    solve_cb = nothing,
    default_initstates = ComponentVector(...)
)
```

**Parameters:**
- `component`: HydroModels component to calibrate
- `input`: Input data matrix (n_inputs × n_timesteps)
- `target`: Observed output vector (n_timesteps)
- `metric`: Optimization metric ("KGE", "NSE", "LogKGE", "MSE")
- `loss_func`: Custom loss function (if metric not specified)
- `fixed_params`: ComponentVector of parameters to fix (supports hierarchical structure)
- `initial_params`: ComponentVector of initial parameter values (optional, defaults to 1.0 for params and random for nns)
- `warm_up`: Number of initial timesteps to skip in loss calculation
- `lb_pas`, `ub_pas`: Lower and upper bounds for calibratable parameters only
- `adtype`: Automatic differentiation type for gradient-based optimization

### get_initial_params

```julia
get_initial_params(component::AbstractComponent; rng=Random.default_rng(), eltype=Float64)
```

Get initial parameters for a component based on its metadata.

**Returns:** ComponentVector with `params` (traditional parameters initialized to 1.0) and `nns` (neural network parameters randomly initialized).

**Example:**
```julia
flux = @hydroflux q ~ k * p
params = get_initial_params(flux)
# Returns: ComponentVector(params=(k=1.0,))

neural_model = HydroModel([neural_flux, bucket])
params = get_initial_params(neural_model)
# Returns: ComponentVector(params=(k=1.0, ...), nns=(et_net=ComponentVector(...),))
```

**Note:** Traditional parameters are initialized to 1.0. For custom initialization, provide `initial_params` to `OptimizationProblem`.

### get_nn_initial_params

```julia
get_nn_initial_params(component::AbstractComponent; rng=Random.default_rng())
```

Extract initial neural network parameters from a component.

**Returns:** NamedTuple with neural network names as keys and ComponentVector parameters as values, or `nothing` if no neural networks found.

**Supported Components:**
- `NeuralFlux`: Returns single network parameters
- `NeuralBucket`: Returns flux, state, and output network parameters
- `HydroBucket`: Extracts from embedded neural flux components
- `HydroModel`: Recursively collects from all sub-components

## Notes

- When using `fixed_params`, provide bounds in `lb_pas`/`ub_pas` for calibratable parameters only
- Metrics are minimized (negative values for KGE/NSE to maximize them)
- The extension integrates seamlessly with all Optimization.jl solvers
- `fixed_params` uses `update_ca` for hierarchical parameter updates, supporting nested structures
