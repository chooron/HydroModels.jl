# Using Optimization to Calibrate HydroModels.jl Parameters

This tutorial will guide you on how to use Julia's Optimization ecosystem to calibrate parameters in HydroModels.jl hydrological models. Through parameter optimization, you can improve model outputs to better fit observed data, thereby enhancing the model's predictive capabilities.

## Table of Contents

1. [Overview](#overview)
2. [Required Packages](#required-packages)
3. [Data Preparation](#data-preparation)
4. [Model Initialization](#model-initialization)
5. [Parameter Boundary Setup](#parameter-boundary-setup)
6. [Objective Function Definition](#objective-function-definition)
7. [Optimization Problem Configuration](#optimization-problem-configuration)
8. [Executing Optimization](#executing-optimization)
9. [Result Analysis and Visualization](#result-analysis-and-visualization)
10. [Choosing Different Optimization Algorithms](#choosing-different-optimization-algorithms)
11. [Advanced Applications](#advanced-applications)

## Overview

Parameter optimization is a core step in hydrological model calibration. In HydroModels.jl, we leverage Julia's Optimization ecosystem to implement efficient parameter calibration. This tutorial will demonstrate how to set up and execute the optimization process to find parameter sets that best fit observed data.

## Required Packages

First, we need to import the necessary packages:

```julia
using CSV, DataFrames, Dates, ComponentArrays
using HydroModels, HydroModelTools
using DataInterpolations
using OptimizationEvolutionary  # Evolutionary algorithms
using OptimizationGCMAES       # Covariance Matrix Adaptation Evolution Strategy
using OptimizationBBO          # Biogeography-Based Optimization
using Optimization             # Core optimization framework
using ModelingToolkit
using Distributions
using ProgressMeter
using Plots
using JLD2
```

## Data Preparation

The optimization process requires input data (such as precipitation, evaporation) and observed data (such as streamflow):

```julia
# Load data
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path)
df = DataFrame(data)

# Select time period
ts = collect(1:10000)

# Extract necessary variables
lday_vec = df[ts, "dayl(day)"]           # Day length
prcp_vec = df[ts, "prcp(mm/day)"]        # Precipitation
temp_vec = df[ts, "tmean(C)"]            # Temperature
flow_vec = df[ts, "flow(mm)"]            # Observed streamflow

# Calculate potential evapotranspiration (PET)
pet_vec = @. 29.8 * lday_vec * 24 * 0.611 * exp((17.3 * temp_vec) / (temp_vec + 237.3)) / (temp_vec + 273.2)

# Set warm-up period and maximum iterations
warm_up = 365    # Model warm-up period (days)
max_iter = 1000  # Maximum optimization iterations
```

## Model Initialization

Next, we initialize the hydrological model and set up initial parameters and states:

```julia
# Select model (e.g., HBV model)
model_nm = "hbv"
model = HydroModelLibrary.load_model(Symbol(model_nm))

# Get parameter boundaries
param_bounds = getbounds.(get_params(model))

# Randomly generate initial parameter values (within parameter boundaries)
random_param_values = map(param_bounds) do param_bound
    rand(Uniform(param_bound[1], param_bound[2]))
end

# Create parameter vector
init_params = ComponentVector(params=NamedTuple{Tuple(get_param_names(model))}(random_param_values))
ps_axes = getaxes(init_params)

# Initialize model states (all state variables set to 0)
init_states = NamedTuple{Tuple(get_state_names(model))}(zeros(length(get_state_names(model)))) |> ComponentVector
```

## Preparing Model Inputs and Run Configuration

Before optimization, we need to prepare model input data and run configuration:

```julia
# Prepare input data
input = (P=prcp_vec, Ep=pet_vec, T=temp_vec)
input_matrix = Matrix(reduce(hcat, collect(input[HydroModels.get_input_names(model)])'))

# Set model run configuration
config = (solver=HydroModelTools.DiscreteSolver(), interp=LinearInterpolation)
run_kwargs = (config=config, initstates=init_states)
```

## Objective Function Definition

The core of optimization is defining an objective function to evaluate the performance of model parameters. Here we use KGE (Kling-Gupta Efficiency) as the evaluation metric:

```julia
# Calculate mean for R² computation
y_mean = mean(flow_vec)

# Define R² function (coefficient of determination)
r2_func(y, y_hat) = sum((y .- y_hat) .^ 2) ./ sum((y .- y_mean) .^ 2)

# Define KGE function (Kling-Gupta Efficiency)
kge_func(y, y_hat) = sqrt((r2_func(y, y_hat))^2 + (std(y_hat) / std(y) - 1)^2 + (mean(y_hat) / mean(y) - 1)^2)

# Define optimization objective function
function obj_func(p, _)
    return kge_func(
        flow_vec[warm_up:end],  # Observed values after warm-up period
        model(input_matrix, ComponentVector(p, ps_axes); run_kwargs...)[end, warm_up:end]  # Model predictions
    )
end
```

## Optimization Progress Tracking

To monitor the optimization process, we set up a progress bar and callback function:

```julia
# Create progress bar
progress = Progress(max_iter, desc="Optimization")

# Create recorder to save results from each iteration
recorder = []

# Define callback function
callback_func!(state, l) = begin
    push!(recorder, (iter=state.iter, loss=l, time=now(), params=state.u))
    next!(progress)  # Update progress bar
    false  # Return false to continue optimization
end
```

## Parameter Boundary Setup

To ensure the optimization algorithm searches within a reasonable parameter space, we need to set upper and lower bounds for parameters:

```julia
# Extract parameter lower bounds
lb_list = first.(param_bounds) .|> eltype(input_matrix)

# Extract parameter upper bounds
ub_list = last.(param_bounds) .|> eltype(input_matrix)
```

## Optimization Problem Configuration

Configure the optimization problem using the Optimization.jl framework:

```julia
# Create optimization function object
optf = Optimization.OptimizationFunction(obj_func, Optimization.AutoForwardDiff())

# Create optimization problem
optprob = Optimization.OptimizationProblem(optf, Vector(init_params), lb=lb_list, ub=ub_list)
```

## Executing Optimization

Select an optimization algorithm and execute the optimization process:

```julia
# Solve optimization problem
sol = Optimization.solve(
    optprob,
    # CMAES(),  # Covariance Matrix Adaptation Evolution Strategy
    BBO_adaptive_de_rand_1_bin_radiuslimited(),  # Adaptive Differential Evolution algorithm
    # GCMAESOpt(),  # Global Covariance Matrix Adaptation Evolution Strategy
    maxiters=max_iter,  # Maximum iteration count
    callback=callback_func!  # Callback function
)

# Save optimization process records
recorder_df = DataFrame(recorder)
CSV.write("cache/recorder_$(model_nm).csv", recorder_df)
```

## Result Analysis and Visualization

After optimization is complete, we can run the model with the optimal parameters and analyze the results:

```julia
# Convert optimization results to ComponentVector format
params = ComponentVector(sol.u, ps_axes)

# Run model with optimal parameters
output = model(input_matrix, params; run_kwargs...)

# Extract model outputs
model_output_names = vcat(get_state_names(model), get_output_names(model)) |> Tuple
output_df = DataFrame(NamedTuple{model_output_names}(eachslice(output, dims=1)))

# Calculate final KGE value
@info loss = kge_func(flow_vec, output_df[!, :q_routed])

# Visualize results
plot(output_df[!, :q_routed], label="Simulated Flow")
plot!(flow_vec, label="Observed Flow")
```

## Choosing Different Optimization Algorithms

HydroModels.jl supports multiple optimization algorithms. You can choose the most suitable algorithm based on your specific problem:

1. **Differential Evolution Algorithm** (BBO_adaptive_de_rand_1_bin_radiuslimited): Suitable for most hydrological model parameter optimization problems, with good global search capabilities.

   ```julia
   sol = Optimization.solve(optprob, BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters=max_iter)
   ```

2. **Covariance Matrix Adaptation Evolution Strategy** (CMAES): Works well for complex nonlinear problems.

   ```julia
   sol = Optimization.solve(optprob, CMAES(), maxiters=max_iter)
   ```

3. **Global Covariance Matrix Adaptation Evolution Strategy** (GCMAESOpt): An improved version of CMAES with better global search capabilities.

   ```julia
   sol = Optimization.solve(optprob, GCMAESOpt(), maxiters=max_iter)
   ```

## Advanced Applications

### Multi-Objective Optimization

For situations requiring simultaneous optimization of multiple objectives, you can define a composite objective function:

```julia
function multi_obj_func(p, _)
    model_output = model(input_matrix, ComponentVector(p, ps_axes); run_kwargs...)[end, warm_up:end]
    
    # Calculate multiple objectives
    kge_value = kge_func(flow_vec[warm_up:end], model_output)
    bias_value = abs(mean(model_output) / mean(flow_vec[warm_up:end]) - 1)
    
    # Combine multiple objectives
    return 0.7 * kge_value + 0.3 * bias_value
end
```

### Parameter Sensitivity Analysis

After optimization, you can conduct parameter sensitivity analysis to understand which parameters have the greatest impact on model performance:

```julia
function sensitivity_analysis(model, optimal_params, param_names)
    results = []
    
    for (i, param_name) in enumerate(param_names)
        # Create parameter perturbation sequence
        param_values = range(lb_list[i], ub_list[i], length=10)
        
        # Evaluate model performance for each perturbation value
        for val in param_values
            # Copy optimal parameters and modify current parameter
            test_params = copy(optimal_params)
            test_params[i] = val
            
            # Run model and evaluate performance
            model_output = model(input_matrix, ComponentVector(test_params, ps_axes); run_kwargs...)[end, warm_up:end]
            kge_val = kge_func(flow_vec[warm_up:end], model_output)
            
            push!(results, (param=param_name, value=val, kge=kge_val))
        end
    end
    
    return DataFrame(results)
end
```

### Ensemble of Multiple Optimization Algorithms

You can try multiple optimization algorithms and select the best result:

```julia
function ensemble_optimization(optprob, algorithms, max_iter)
    best_sol = nothing
    best_loss = Inf
    
    for alg in algorithms
        sol = Optimization.solve(optprob, alg, maxiters=max_iter)
        loss = sol.minimum
        
        if loss < best_loss
            best_loss = loss
            best_sol = sol
        end
    end
    
    return best_sol
end

# Usage example
algorithms = [BBO_adaptive_de_rand_1_bin_radiuslimited(), CMAES(), GCMAESOpt()]
best_sol = ensemble_optimization(optprob, algorithms, max_iter)
```

## Summary

This tutorial demonstrates how to use Julia's Optimization ecosystem to calibrate hydrological models in HydroModels.jl. By defining appropriate objective functions, selecting suitable optimization algorithms, and setting reasonable parameter boundaries, you can effectively calibrate hydrological models to better fit observed data.

Optimization is an iterative process that may require trying different objective functions, optimization algorithms, and initial parameter values to achieve the best results. Using the methods in this tutorial, you can flexibly customize the optimization process to meet specific hydrological modeling requirements.