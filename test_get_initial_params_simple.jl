using HydroModels
using ComponentArrays
using Lux
using Random

println("=" ^ 60)
println("Testing get_initial_params")
println("=" ^ 60)

# Test 1: Simple flux with parameters
println("\n[Test 1] HydroFlux with parameters")
@variables prcp q
@parameters k b

flux = @hydroflux q ~ k * prcp + b
params = get_initial_params(flux; rng=Random.MersenneTwister(42))

println("  Params: ", params)
println("  k = ", params.params.k, " (default: 1.0)")
println("  b = ", params.params.b, " (default: 1.0)")
println("  ✓ Parameters initialized to 1.0")

# Test 2: Neural flux
println("\n[Test 2] NeuralFlux with neural network")
@variables prcp temp et
nn = Chain(Dense(2 => 10, tanh), Dense(10 => 1), name=:et_net)
neural_flux = @neuralflux et ~ nn([prcp, temp])

params = get_initial_params(neural_flux; rng=Random.MersenneTwister(42))
println("  Params structure: ", keys(params))
if :nns in keys(params)
    println("  Neural network names: ", keys(params.nns))
    println("  NN param count: ", length(params.nns.et_net))
    println("  ✓ Has neural network parameters")
end

# Test 3: get_nn_initial_params
println("\n[Test 3] get_nn_initial_params")
nn_params = get_nn_initial_params(neural_flux; rng=Random.MersenneTwister(42))
println("  NN params: ", keys(nn_params))
println("  ✓ Extracted NN params: ", !isnothing(nn_params))

# Test 4: Reproducibility
println("\n[Test 4] Reproducibility with RNG")
@variables p q2
@parameters k2
flux2 = @hydroflux q2 ~ k2 * p

rng1 = Random.MersenneTwister(123)
rng2 = Random.MersenneTwister(123)

# For traditional params, should be same (1.0)
params1 = get_initial_params(flux2; rng=rng1)
params2 = get_initial_params(flux2; rng=rng2)

println("  Params1 k2: ", params1.params.k2)
println("  Params2 k2: ", params2.params.k2)
println("  ✓ Reproducible: ", params1.params.k2 == params2.params.k2)

# Test 5: Neural network reproducibility
println("\n[Test 5] Neural network reproducibility")
rng3 = Random.MersenneTwister(456)
rng4 = Random.MersenneTwister(456)

nn_params1 = get_nn_initial_params(neural_flux; rng=rng3)
nn_params2 = get_nn_initial_params(neural_flux; rng=rng4)

println("  NN params equal: ", nn_params1.et_net == nn_params2.et_net)
println("  ✓ Neural network params reproducible")

println("\n" * "=" ^ 60)
println("All tests completed!")
println("=" ^ 60)
