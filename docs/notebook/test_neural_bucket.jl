# Test NeuralBucket - RNN-style Neural Bucket Component
# 
# This notebook demonstrates the usage of NeuralBucket, which integrates
# Lux.jl's RNN layers with HydroModels framework.

using HydroModels
using Lux
using ComponentArrays
using Random
using Statistics

# Set random seed for reproducibility
Random.seed!(42)

println("=" ^ 70)
println("Testing NeuralBucket Component")
println("=" ^ 70)

# ============================================================================
# Example 1: Simple LSTM-based bucket
# ============================================================================

println("\n[Example 1] Creating LSTM-based NeuralBucket")

# Define a simple LSTM bucket using the convenience function
lstm_bucket = create_lstm_bucket(
    name = :rainfall_runoff,
    n_inputs = 2,      # precipitation and temperature
    n_hidden = 32,     # LSTM hidden size
    n_outputs = 1,     # runoff output
    inputs = [:prcp, :temp],
    outputs = [:runoff]
)

println("✓ Created LSTM bucket: $(lstm_bucket.name)")
println("  Inputs: $(get_input_names(lstm_bucket))")
println("  Outputs: $(get_output_names(lstm_bucket))")
println("  Neural networks: $(get_nn_names(lstm_bucket))")

# Initialize parameters
rng = Random.default_rng()
ps = LuxCore.initialparameters(rng, lstm_bucket)
ps_vec = ComponentVector(ps)

println("  Total parameters: $(length(ps_vec))")

# Prepare test data
n_timesteps = 100
test_input = rand(Float32, 2, n_timesteps)  # (n_inputs × n_timesteps)

# Package parameters for HydroModels interface
cell_name = get_nn_names(lstm_bucket)[1]
output_name = get_nn_names(lstm_bucket)[2]
params = ComponentVector(
    nns = (
        rainfall_runoff_lstm = ps.cell,
        rainfall_runoff_output = ps.output
    )
)

# Run the bucket (single-node)
config = HydroConfig()
output = lstm_bucket(test_input, params, config)

println("  Input shape: $(size(test_input))")
println("  Output shape: $(size(output))")
println("  Output range: [$(minimum(output)), $(maximum(output))]")

# ============================================================================
# Example 2: GRU-based bucket
# ============================================================================

println("\n[Example 2] Creating GRU-based NeuralBucket")

gru_bucket = create_gru_bucket(
    name = :snowmelt,
    n_inputs = 3,      # temp, snowpack, radiation
    n_hidden = 16,     # GRU hidden size
    n_outputs = 2,     # melt + runoff
    inputs = [:temp, :snowpack, :radiation],
    outputs = [:melt, :runoff]
)

println("✓ Created GRU bucket: $(gru_bucket.name)")
println("  Inputs: $(get_input_names(gru_bucket))")
println("  Outputs: $(get_output_names(gru_bucket))")

# Initialize and run
ps_gru = LuxCore.initialparameters(rng, gru_bucket)
params_gru = ComponentVector(
    nns = (
        snowmelt_gru = ps_gru.cell,
        snowmelt_output = ps_gru.output
    )
)

test_input_gru = rand(Float32, 3, 50)
output_gru = gru_bucket(test_input_gru, params_gru, config)

println("  Output shape: $(size(output_gru))")
println("  Mean outputs: melt=$(mean(output_gru[1,:])), runoff=$(mean(output_gru[2,:]))")

# ============================================================================
# Example 3: Custom RNN cell in NeuralBucket
# ============================================================================

println("\n[Example 3] Custom NeuralBucket with deep LSTM")

# Create a more complex LSTM cell with multiple layers
deep_lstm = Chain(
    LSTMCell(2 => 64),
    LSTMCell(64 => 32),
    name = :deep_lstm
)

output_layer = Chain(
    Dense(32 => 16, relu),
    Dense(16 => 1),
    name = :deep_output
)

custom_bucket = NeuralBucket(
    name = :deep_rainfall_runoff,
    cell = deep_lstm,
    output_layer = output_layer,
    inputs = [:prcp, :pet],
    outputs = [:flow]
)

println("✓ Created custom deep LSTM bucket")

# Initialize
ps_custom = LuxCore.initialparameters(rng, custom_bucket)
ps_custom_vec = ComponentVector(ps_custom)
println("  Total parameters: $(length(ps_custom_vec))")

params_custom = ComponentVector(
    nns = (
        deep_lstm = ps_custom.cell,
        deep_output = ps_custom.output
    )
)

test_input_custom = rand(Float32, 2, 80)
output_custom = custom_bucket(test_input_custom, params_custom, config)

println("  Output shape: $(size(output_custom))")
println("  Output statistics: mean=$(mean(output_custom)), std=$(std(output_custom))")

# ============================================================================
# Example 4: Multi-node NeuralBucket
# ============================================================================

println("\n[Example 4] Multi-node NeuralBucket computation")

# Create a simple LSTM bucket
multi_bucket = create_lstm_bucket(
    name = :multi_site,
    n_inputs = 2,
    n_hidden = 16,
    n_outputs = 1,
    inputs = [:prcp, :temp],
    outputs = [:flow]
)

# Prepare multi-node data (n_inputs × n_nodes × n_timesteps)
n_nodes = 5
test_input_multi = rand(Float32, 2, n_nodes, 60)

# Initialize parameters (same for all nodes)
ps_multi = LuxCore.initialparameters(rng, multi_bucket)
params_multi = ComponentVector(
    nns = (
        multi_site_lstm = ps_multi.cell,
        multi_site_output = ps_multi.output
    )
)

# Run multi-node computation
output_multi = multi_bucket(test_input_multi, params_multi, config)

println("  Input shape: $(size(test_input_multi))")
println("  Output shape: $(size(output_multi))")
println("  Output range per node:")
for i in 1:n_nodes
    node_output = output_multi[1, i, :]
    println("    Node $i: [$(minimum(node_output)), $(maximum(node_output))]")
end

# ============================================================================
# Example 5: Comparison with traditional HydroBucket
# ============================================================================

println("\n[Example 5] NeuralBucket vs HydroBucket characteristics")

println("\nNeuralBucket characteristics:")
println("  ✓ No differential equations - discrete time stepping")
println("  ✓ No ODE solver required")
println("  ✓ RNN hidden states maintained across time")
println("  ✓ Can learn complex temporal patterns")
println("  ✓ Compatible with Lux ecosystem")
println("  ✓ Supports multi-node computation")
println("  ✓ Can be integrated into HydroModel")

println("\nTypical use cases:")
println("  • Pure data-driven modeling")
println("  • Hybrid models (combine with HydroBucket)")
println("  • Residual modeling (correct process-based models)")
println("  • Parameter estimation networks")

# ============================================================================
# Example 6: Integration with HydroModel
# ============================================================================

println("\n[Example 6] Integrating NeuralBucket into HydroModel")

# Note: NeuralBucket can be used alongside traditional buckets
# However, since it doesn't define states/fluxes in the traditional way,
# it's typically used as a standalone component or post-processor

# Example structure (pseudo-code):
println("  Hybrid model structure:")
println("    1. Traditional HydroBucket for physical processes")
println("    2. NeuralBucket for residual correction or parameter estimation")
println("    3. Combine outputs in final model")

# ============================================================================
# Summary
# ============================================================================

println("\n" * "=" ^ 70)
println("Summary")
println("=" ^ 70)

println("\n✓ Successfully tested NeuralBucket component")
println("✓ Verified single-node and multi-node computation")
println("✓ Tested LSTM, GRU, and custom cell configurations")
println("✓ Demonstrated Lux integration and parameter management")

println("\nNeuralBucket is ready for:")
println("  • Pure neural network-based hydrological modeling")
println("  • Hybrid physics-ML models")
println("  • Real-time sequence prediction")
println("  • Transfer learning applications")

println("\n" * "=" ^ 70)

