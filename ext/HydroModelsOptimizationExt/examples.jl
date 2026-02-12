# Examples for HydroModelsOptimizationExt

using HydroModels
using Optimization
using OptimizationBBO
using ComponentArrays
using Lux

# Example 1: Basic calibration with KGE metric
function example_basic_kge()
    # Load your component and data
    component = YourHydroComponent()
    input = rand(5, 100)  # 5 inputs, 100 timesteps
    target = rand(100)    # Observed output

    # Create optimization problem with KGE metric
    prob = OptimizationProblem(
        component,
        input,
        target;
        metric="KGE",  # Use KGE instead of MSE
        warm_up=10,
        lb_pas=[0.1, 0.1, 0.1],
        ub_pas=[10.0, 10.0, 10.0]
    )

    # Solve
    sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=1000)
    return sol
end

# Example 2: Calibration with fixed parameters (ComponentVector interface)
function example_fixed_params()
    component = YourHydroComponent()
    input = rand(5, 100)
    target = rand(100)

    # Fix some parameters using ComponentVector
    # Suppose component has parameters: [:k1, :k2, :k3]
    prob = OptimizationProblem(
        component,
        input,
        target;
        metric="NSE",
        fixed_params=ComponentVector(params=(k2=2.5,)),  # Fix k2, calibrate k1 and k3
        lb_pas=[0.1, 0.1],  # Bounds for k1 and k3 only
        ub_pas=[10.0, 10.0]
    )

    sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=1000)
    return sol
end

# Example 3: Model with NeuralFlux - Fix traditional parameters
function example_neural_flux()
    @variables prcp temp et q

    # Create neural network
    nn = Chain(Dense(2 => 10, tanh), Dense(10 => 1), name=:et_net)
    neural_flux = @neuralflux et ~ nn([prcp, temp])

    # Build model with neural flux (add other components as needed)
    # model = HydroModel([neural_flux, other_components...])
    model = neural_flux  # Simplified for example

    input = rand(2, 100)
    target = rand(100)

    # Neural network parameters are automatically initialized
    # No need to manually extract nn_params!
    prob = OptimizationProblem(
        model,
        input,
        target;
        fixed_params=ComponentVector(
            params=(k=0.5, b=2.0)  # Fix traditional parameters
        ),
        lb_pas=[0.1, 0.0],  # Bounds for calibratable params
        ub_pas=[5.0, 10.0]
    )

    sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=1000)
    return sol
end

# Example 4: Fix both traditional and neural network parameters
function example_fix_both()
    model = YourHybridModel()
    input = rand(5, 100)
    target = rand(100)

    # Get neural network parameters to fix them
    nn_params = get_nn_initial_params(model)

    # Fix some traditional params AND some neural network params
    prob = OptimizationProblem(
        model,
        input,
        target;
        fixed_params=ComponentVector(
            params=(k=0.5, b=2.0),           # Fix traditional parameters
            nns=(et_net=nn_params.et_net,)   # Fix neural network parameters
        ),
        lb_pas=[0.1],  # Only bounds for remaining calibratable params
        ub_pas=[5.0]
    )

    sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=1000)
    return sol
end

# Example 5: HydroBucket with NeuralFlux
function example_bucket_with_neural()
    @variables prcp temp snowpack soilwater et q

    # Create neural flux for ET
    et_nn = Chain(Dense(3 => 16, relu), Dense(16 => 1), name=:et_nn)
    et_flux = @neuralflux et ~ et_nn([prcp, temp, soilwater])

    # Create other fluxes (simplified for example)
    # Add your other flux and dflux components here

    # Create bucket with neural flux
    # bucket = @hydrobucket :my_bucket begin
    #     fluxes = [et_flux, other_fluxes...]
    #     dfluxes = [state_fluxes...]
    # end
    bucket = et_flux  # Simplified for example

    input = rand(3, 100)
    target = rand(100)

    # Neural network parameters are automatically initialized
    prob = OptimizationProblem(
        bucket,
        input,
        target;
        fixed_params=ComponentVector(
            params=(snowmelt_factor=1.5,)  # Fix some traditional params
        )
    )

    sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=1000)
    return sol
end

# Example 6: NeuralBucket calibration
function example_neural_bucket()
    # Create neural bucket
    neural_bucket = create_neural_bucket(
        name=:nb,
        n_inputs=3,
        n_states=2,
        n_outputs=1,
        inputs=[:prcp, :temp, :pet],
        states=[:snowpack, :soilwater],
        outputs=[:q]
    )

    input = rand(3, 100)
    target = rand(100)

    # Get neural network parameters to fix them
    nn_params = get_nn_initial_params(neural_bucket)

    # Calibrate only traditional params, fix neural networks
    prob = OptimizationProblem(
        neural_bucket,
        input,
        target;
        fixed_params=ComponentVector(
            nns=nn_params  # Fix all neural network parameters
        )
    )

    sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=1000)
    return sol
end

# Example 7: Using LogKGE for low-flow calibration
function example_logkge()
    component = YourHydroComponent()
    input = rand(5, 100)
    target = rand(100)

    prob = OptimizationProblem(
        component,
        input,
        target;
        metric="LogKGE",  # Better for low flows
        warm_up=30,
        lb_pas=[0.1, 0.1, 0.1],
        ub_pas=[10.0, 10.0, 10.0]
    )

    sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=1000)
    return sol
end

# Example 8: Custom loss function (backward compatible)
function example_custom_loss()
    component = YourHydroComponent()
    input = rand(5, 100)
    target = rand(100)

    # Define custom loss
    custom_loss = function(obs, sim)
        # Weighted combination of NSE and volume error
        nse = 1 - sum((obs .- sim) .^ 2) / sum((obs .- mean(obs)) .^ 2)
        vol_error = abs(sum(obs) - sum(sim)) / sum(obs)
        return -nse + 0.5 * vol_error
    end

    prob = OptimizationProblem(
        component,
        input,
        target;
        loss_func=custom_loss,  # Still works!
        lb_pas=[0.1, 0.1, 0.1],
        ub_pas=[10.0, 10.0, 10.0]
    )

    sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=1000)
    return sol
end

# Example 9: Hierarchical parameter fixing
function example_hierarchical_fixing()
    model = YourHybridModel()
    input = rand(5, 100)
    target = rand(100)

    # Get neural network parameters
    nn_params = get_nn_initial_params(model)

    # Fix only specific neural network layers using hierarchical structure
    fixed_params = ComponentVector(
        params=(k=0.5,),
        nns=(
            et_net=ComponentVector(
                layer_1=nn_params.et_net.layer_1  # Fix only first layer
                # layer_2 will be calibrated
            ),
        )
    )

    prob = OptimizationProblem(
        model,
        input,
        target;
        fixed_params=fixed_params
    )

    sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); maxiters=1000)
    return sol
end

# Example 10: Get initial parameters for inspection
function example_get_initial_params()
    @variables prcp temp et q
    nn = Chain(Dense(2 => 10, tanh), Dense(10 => 1), name=:et_net)
    neural_flux = @neuralflux et ~ nn([prcp, temp])

    # Get complete initial parameters
    params = get_initial_params(neural_flux)
    # Returns: ComponentVector(params=(), nns=(et_net=ComponentVector(...),))

    # Inspect parameter structure
    println("Parameter names: ", keys(params))
    if haskey(getaxes(params)[1], :nns)
        println("Neural network names: ", keys(params.nns))
    end

    return params
end

