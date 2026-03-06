---
name: hydromodels-ext-modeling
description: Build, execute, and calibrate HydroModels.jl models from scratch using HydroModelsYAMLExt, HydroModelsOrdinaryDiffEqExt, HydroModelsDataInterpolationsExt, and HydroModelsOptimizationExt. Use for end-to-end modeling workflows: symbolic model construction, forcing preparation, solver/interpolator setup, YAML-driven execution, optimizer selection, and parameter calibration without GPU. Exclude HydroModelsBMIExt and HydroModelsRastersExt.
---

# HydroModels Extension Modeling (No GPU)

Follow this workflow to build and run hydrological models in this project.

## Apply hard constraints

- Use only these extensions: `HydroModelsYAMLExt`, `HydroModelsOrdinaryDiffEqExt`, `HydroModelsDataInterpolationsExt`, `HydroModelsOptimizationExt`.
- Exclude `HydroModelsBMIExt` and `HydroModelsRastersExt`.
- Keep computation on CPU:
  - Do not import CUDA/GPU toolchains.
  - Keep `device=identity` (or omit `device`).

## Build a model from scratch (symbolic workflow)

1. Define symbolic variables and parameters with `@variables` and `@parameters`.
2. Build fluxes with `@hydroflux` and state equations with `@stateflux`.
3. Compose buckets with `@hydrobucket`.
4. Compose the full model with `@hydromodel`.
5. Build `params` and `initstates` as `ComponentVector`.
6. Build a `HydroConfig` and execute with `model(input, params, config; initstates=initstates)`.

Use this minimal pattern:

```julia
using HydroModels
using ComponentArrays

@variables prcp pet soilwater evap runoff
@parameters k smax

bucket = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux evap ~ min(pet, soilwater)
        @hydroflux runoff ~ k * max(0.0, soilwater - smax)
    end
    dfluxes = begin
        @stateflux soilwater ~ prcp - evap - runoff
    end
end

model = @hydromodel :demo begin
    bucket
end

params = ComponentVector(params=(k=0.4, smax=100.0))
initstates = ComponentVector(soilwater=40.0)
input = rand(2, 365)  # rows = [prcp, pet], cols = time

config = HydroConfig(
    solver=MutableSolver,
    interpolator=Val(ConstantInterpolation),
    timeidx=collect(1:365),
    min_value=1e-6,
    parallel=false,
)

output = model(input, params, config; initstates=initstates)
```

## Execute models correctly

- Keep input shape as `variables × timesteps` (2D) or `variables × nodes × timesteps` (3D).
- Keep variable order aligned with `get_input_names(model)`.
- Prefer explicit `HydroConfig` over loose named tuples.
- Use keyword `initstates=...` when states exist.

## Configure solver behavior

Use solver according to task:

- `MutableSolver`: fastest default for forward runs.
- `ImmutableSolver`: preferred for AD-sensitive calibration workflows.
- `ODESolver`: use DifferentialEquations-based integration via `HydroModelsOrdinaryDiffEqExt`.
- `DiscreteSolver`: use map/discrete integration via SciML discrete problem path.

To trigger ODE/discrete extension methods, load:

```julia
using HydroModels
using OrdinaryDiffEq
using SciMLSensitivity
```

Tune with config keys passed through optimization config or named tuples:

- `solve_alg` (e.g., `Tsit5()`)
- `sense_alg`
- `solve_cb`

## Choose interpolation method

### Built-in interpolation (core)

- `Val(ConstantInterpolation)`: default and fastest for daily/step-like forcings.
- `Val(LinearInterpolation)`: smoother forcings and smoother gradients.

### DataInterpolations extension

Load extension and use:

```julia
using HydroModels
using DataInterpolations

config_linear = HydroConfig(interpolator=Val(DataInterpolations.LinearInterpolation))
config_cubic  = HydroConfig(interpolator=Val(DataInterpolations.CubicSpline))
```

Selection rule:

- Use `ConstantInterpolation` for performance and daily step forcings.
- Use `LinearInterpolation` for smoother calibration gradients.
- Use `DataInterpolations.CubicSpline` only when smooth forcing reconstruction is necessary.

## Run YAML-driven workflows

1. Trigger extension:

```julia
using HydroModels
using YAML
```

2. Load model/config/parameter metadata:

```julia
model = load_model_from_yaml("model.yaml")
config = load_config_from_yaml("model.yaml")
meta = load_parameters_from_yaml("model.yaml")
```

3. Prepare forcing matrix explicitly as `variables × timesteps` and run with explicit params/initstates.

Prefer explicit execution (`load_model_from_yaml` + manual run) over opaque wrappers when debugging dimensions, parameter mapping, or initialization.

## Perform parameter optimization

### Create optimization problem

Use the extension API:

```julia
using HydroModels
using Optimization

prob = OptimizationProblem(
    model,
    input,
    target;
    metric="KGE",              # KGE, NSE, LogKGE, MSE
    warm_up=30,
    solver=ImmutableSolver,
    interpolator=Val(LinearInterpolation),
    lb_pas=[0.1, 10.0],
    ub_pas=[2.0, 500.0],
)
```

### Use different optimizers

Load algorithm packages based on method class:

- Gradient-based:

```julia
using OptimizationOptimJL
sol = solve(prob, BFGS())
```

- Derivative-free/global:

```julia
using OptimizationBBO
sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=1000)
```

### Fix subset of parameters

Use `fixed_params` to calibrate only selected parameters:

```julia
fixed = ComponentVector(params=(k=0.4,))
prob = OptimizationProblem(
    model, input, target;
    metric="NSE",
    fixed_params=fixed,
    lb_pas=[10.0],
    ub_pas=[500.0],
)
```

### Initialize parameter vectors safely

- Use `get_initial_params(model)` for robust default initialization.
- Pass `initial_params=...` into `OptimizationProblem` when custom initialization is required.

## Validate each workflow stage

Run checks in this order:

1. Build model and print `get_input_names/get_param_names/get_state_names`.
2. Run a short forward simulation (e.g., 10-30 steps).
3. Verify output dimensions and absence of NaN/Inf.
4. Run optimization with small `maxiters` as a smoke test.
5. Run full calibration only after smoke tests pass.

## Avoid common mistakes

- Do not use GPU-specific configuration.
- Do not pass wrongly oriented forcing matrices.
- Do not mix old solver naming patterns with current `MutableSolver/ImmutableSolver/ODESolver/DiscreteSolver`.
- Do not provide full-parameter bounds when using `fixed_params`; provide bounds for calibratable parameters only.