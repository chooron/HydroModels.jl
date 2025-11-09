# Test NeuralBucket - ODE-like RNN Structure
# 
# This notebook demonstrates the NEW NeuralBucket design, which mimics
# the ODE structure of HydroBucket using neural networks and Lux.jl's
# recurrence mechanism.

using HydroModels
using Lux
using ComponentArrays
using Random
using Statistics

# Set random seed for reproducibility
Random.seed!(42)

println("=" ^ 70)
println("Testing NeuralBucket - ODE-like Structure")
println("=" ^ 70)

# ============================================================================
# Concept: HydroBucket vs NeuralBucket
# ============================================================================

println("\n[Concept] Structural Similarity between ODE and RNN")
println()
println("Traditional HydroBucket (ODE):")
println("  dS/dt = f(S, input, params)")
println("  S_t = S_{t-1} + dt * f(S_{t-1}, input_t, params)")
println("  output_t = g(S_t, input_t, params)")
println()
println("NeuralBucket (RNN-like):")
println("  flux_t = flux_network(S_{t-1}, input_t)")
println("  S_t = state_network(S_{t-1}, flux_t)")
println("  output_t = output_network(flux_t)")
println()
println("Both maintain states across time and iterate over sequences!")

# ============================================================================
# Example 1: Simple NeuralBucket with linear transformations
# ============================================================================

println("\n" * "=" ^ 70)
println("[Example 1] Simple NeuralBucket (Linear)")
println("=" ^ 70)

# Create a simple bucket
simple_bucket = create_simple_neural_bucket(
    name = :simple_bucket,
    n_inputs = 2,       # precipitation, temperature
    n_states = 1,       # storage
    n_outputs = 1,      # runoff
    inputs = [:prcp, :temp],
    states = [:storage],
    outputs = [:runoff]
)

println("\n✓ Created simple NeuralBucket: $(simple_bucket.name)")
println("  Inputs: $(get_input_names(simple_bucket))")
println("  States: $(get_state_names(simple_bucket))")
println("  Outputs: $(get_output_names(simple_bucket))")
println("  Neural networks: $(get_nn_names(simple_bucket))")

# Initialize parameters
rng = Random.default_rng()
ps_lux = LuxCore.initialparameters(rng, simple_bucket)
ps_vec = ComponentVector(ps_lux)

println("\n  Total parameters: $(length(ps_vec))")
println("  Flux network params: $(length(ComponentVector(ps_lux.flux)))")
println("  State network params: $(length(ComponentVector(ps_lux.state)))")
println("  Output network params: $(length(ComponentVector(ps_lux.output)))")

# Package parameters for HydroModels interface
params = ComponentVector(
    nns = (
        simple_bucket_flux = ps_lux.flux,
        simple_bucket_state = ps_lux.state,
        simple_bucket_output = ps_lux.output
    )
)

# Prepare test data
n_timesteps = 100
test_input = rand(Float32, 2, n_timesteps)  # (2 × 100)

# Run the bucket
config = HydroConfig()
output = simple_bucket(test_input, params, config)

println("\n  Input shape: $(size(test_input))")
println("  Output shape: $(size(output))")
println("  Output range: [$(minimum(output)), $(maximum(output))]")
println("  Output mean: $(mean(output))")

# ============================================================================
# Example 2: NeuralBucket with hidden layers
# ============================================================================

println("\n" * "=" ^ 70)
println("[Example 2] NeuralBucket with Hidden Layers")
println("=" ^ 70)

# Create a bucket with hidden layers
bucket = create_neural_bucket(
    name = :rainfall_runoff,
    n_inputs = 2,           # prcp, temp
    n_states = 2,           # snowpack, soilwater
    n_outputs = 1,          # runoff
    n_fluxes = 3,           # intermediate fluxes
    hidden_size = 16,
    inputs = [:prcp, :temp],
    states = [:snowpack, :soilwater],
    outputs = [:runoff],
    flux_activation = tanh,
    state_activation = identity,
    output_activation = identity
)

println("\n✓ Created NeuralBucket: $(bucket.name)")
println("  n_inputs: $(bucket.n_inputs)")
println("  n_states: $(bucket.n_states)")
println("  n_outputs: $(bucket.n_outputs)")

# Initialize
ps_bucket = LuxCore.initialparameters(rng, bucket)
ps_bucket_vec = ComponentVector(ps_bucket)

println("\n  Total parameters: $(length(ps_bucket_vec))")

params_bucket = ComponentVector(
    nns = (
        rainfall_runoff_flux = ps_bucket.flux,
        rainfall_runoff_state = ps_bucket.state,
        rainfall_runoff_output = ps_bucket.output
    )
)

# Run
output_bucket = bucket(test_input, params_bucket, config)

println("\n  Output shape: $(size(output_bucket))")
println("  Output statistics: mean=$(mean(output_bucket)), std=$(std(output_bucket))")

# ============================================================================
# Example 3: Single time step (Lux interface)
# ============================================================================

println("\n" * "=" ^ 70)
println("[Example 3] Single Time Step Computation (Lux Interface)")
println("=" ^ 70)

# Initialize Lux states
st = LuxCore.initialstates(rng, bucket)

# Prepare single time step input
x_t = test_input[:, 1]  # First time step

# Lux interface: (x_t, ps, st) → (y_t, st_new)
y_t, st_new = bucket(x_t, ps_bucket, st)

println("\n✓ Single time step computation:")
println("  Input shape: $(size(x_t))")
println("  Output shape: $(size(y_t))")
println("  Output value: $y_t")
println("\n  Initial hydro_state: $(st.hydro_state)")
println("  Updated hydro_state: $(st_new.hydro_state)")

# Run another step to see state evolution
y_t2, st_new2 = bucket(test_input[:, 2], ps_bucket, st_new)

println("\n  After second step:")
println("  Output value: $y_t2")
println("  Updated hydro_state: $(st_new2.hydro_state)")

# ============================================================================
# Example 4: With initial states
# ============================================================================

println("\n" * "=" ^ 70)
println("[Example 4] Running with Initial States")
println("=" ^ 70)

# Prepare initial states
initstates = ComponentVector(
    states = (
        snowpack = [0.5f0],
        soilwater = [10.0f0]
    )
)

println("\n  Initial states: snowpack=$(initstates[:states].snowpack), soilwater=$(initstates[:states].soilwater)")

# Run with initial states
output_with_init = bucket(test_input, params_bucket, config; initstates=initstates)

println("\n  Output shape: $(size(output_with_init))")
println("  Output statistics: mean=$(mean(output_with_init)), std=$(std(output_with_init))")

# Compare with zero initial states
output_zero_init = bucket(test_input, params_bucket, config)

println("\n  Difference from zero initial states:")
println("  Mean absolute difference: $(mean(abs.(output_with_init .- output_zero_init)))")

# ============================================================================
# Example 5: Multi-node computation
# ============================================================================

println("\n" * "=" ^ 70)
println("[Example 5] Multi-node Computation")
println("=" ^ 70)

# Prepare multi-node data
n_nodes = 5
test_input_multi = rand(Float32, 2, n_nodes, 100)  # (2 × 5 × 100)

println("\n  Multi-node input shape: $(size(test_input_multi))")

# Run multi-node computation
output_multi = bucket(test_input_multi, params_bucket, config)

println("\n✓ Multi-node computation completed")
println("  Output shape: $(size(output_multi))")
println("  Output statistics per node:")
for i in 1:n_nodes
    node_output = output_multi[1, i, :]
    println("    Node $i: mean=$(round(mean(node_output), digits=4)), std=$(round(std(node_output), digits=4))")
end

# ============================================================================
# Example 6: Custom network architecture
# ============================================================================

println("\n" * "=" ^ 70)
println("[Example 6] Custom Network Architecture")
println("=" ^ 70)

# Define custom networks
custom_flux_net = Chain(
    Dense(3 => 32, relu),    # 1 state + 2 inputs
    Dense(32 => 16, tanh),
    Dense(16 => 2),          # 2 fluxes
    name = :custom_flux
)

custom_state_net = Chain(
    Dense(3 => 16, relu),    # 1 state + 2 fluxes
    Dense(16 => 1),          # 1 new state
    name = :custom_state
)

custom_output_net = Chain(
    Dense(2 => 8, relu),     # 2 fluxes
    Dense(8 => 1),           # 1 output
    name = :custom_output
)

custom_bucket = NeuralBucket(;
    name = :custom,
    flux_network = custom_flux_net,
    state_network = custom_state_net,
    output_network = custom_output_net,
    n_inputs = 2,
    n_states = 1,
    n_outputs = 1,
    inputs = [:prcp, :temp],
    states = [:storage],
    outputs = [:flow]
)

println("\n✓ Created custom NeuralBucket")

# Initialize and run
ps_custom = LuxCore.initialparameters(rng, custom_bucket)
params_custom = ComponentVector(
    nns = (
        custom_flux = ps_custom.flux,
        custom_state = ps_custom.state,
        custom_output = ps_custom.output
    )
)

output_custom = custom_bucket(test_input, params_custom, config)

println("  Output shape: $(size(output_custom))")
println("  Total parameters: $(length(ComponentVector(ps_custom)))")

# ============================================================================
# Example 7: Comparison with traditional patterns
# ============================================================================

println("\n" * "=" ^ 70)
println("[Example 7] ODE vs RNN Structure Comparison")
println("=" ^ 70)

println("\nTraditional HydroBucket:")
println("  ✓ Uses symbolic expressions (e.g., @hydroflux)")
println("  ✓ Requires ODE solver (MutableSolver, ImmutableSolver, etc.)")
println("  ✓ Explicit physical equations")
println("  ✓ States have physical meaning")
println("  ✓ High interpretability")

println("\nNeuralBucket:")
println("  ✓ Uses neural networks (Lux layers)")
println("  ✓ No ODE solver needed (discrete time stepping)")
println("  ✓ Learned transformations")
println("  ✓ States can be learned representations")
println("  ✓ More flexible, less interpretable")

println("\nCommon features:")
println("  ✓ Both maintain states across time")
println("  ✓ Both iterate over input sequences")
println("  ✓ Both support multi-node computation")
println("  ✓ Both compatible with HydroModels interface")
println("  ✓ Both support Zygote automatic differentiation")

# ============================================================================
# Example 8: State evolution visualization
# ============================================================================

println("\n" * "=" ^ 70)
println("[Example 8] State Evolution Analysis")
println("=" ^ 70)

# Track state evolution manually
st_evolve = LuxCore.initialstates(rng, simple_bucket)
state_history = Float32[]

for t in 1:20
    x_t = test_input[:, t]
    y_t, st_evolve = simple_bucket(x_t, ps_lux, st_evolve)
    
    if !isnothing(st_evolve.hydro_state)
        push!(state_history, st_evolve.hydro_state[1])
    end
end

println("\n  State evolution (first 20 steps):")
for (i, s) in enumerate(state_history[1:min(10, length(state_history))])
    println("    Step $i: state = $(round(s, digits=4))")
end

println("\n  State statistics:")
println("    Mean: $(round(mean(state_history), digits=4))")
println("    Std: $(round(std(state_history), digits=4))")
println("    Min: $(round(minimum(state_history), digits=4))")
println("    Max: $(round(maximum(state_history), digits=4))")

# ============================================================================
# Summary
# ============================================================================

println("\n" * "=" ^ 70)
println("Summary")
println("=" ^ 70)

println("\n✓ NeuralBucket successfully mimics HydroBucket structure")
println("✓ Uses Lux.jl for parameter management and differentiation")
println("✓ Implements three-stage computation:")
println("  1. Flux computation: (states, inputs) → fluxes")
println("  2. State update: (states, fluxes) → new_states")
println("  3. Output computation: fluxes → outputs")

println("\nKey advantages:")
println("  ✓ RNN-like recurrence structure")
println("  ✓ No explicit ODE solver required")
println("  ✓ Flexible neural network architectures")
println("  ✓ Compatible with Lux ecosystem")
println("  ✓ Maintains hydrological states like traditional buckets")

println("\nUse cases:")
println("  • Pure data-driven hydrological modeling")
println("  • Hybrid physics-ML models (combine with HydroBucket)")
println("  • Learning complex state update rules")
println("  • Transfer learning from other domains")

println("\n" * "=" ^ 70)

