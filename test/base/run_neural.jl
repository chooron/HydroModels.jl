@testset "NeuralBucket" begin
    # Define symbolic variables for the bucket
    @variables lday temp prcp snowfall rainfall melt snowpack snowpack2
    @parameters f Smax Qmax Df Tmax Tmin

    # Define the step function used in the model
    step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

    # Build a neural bucket using the macro syntax
    @testset "NeuralBucket Construction" begin
        # Test bucket construction with macro
        bucket_macro = @neuralbucket :bucket1 begin
            fluxes = begin
                @hydroflux begin
                    snowfall ~ step_func(Tmin - temp) * prcp
                    rainfall ~ step_func(temp - Tmin) * prcp
                end
                @hydroflux melt ~ step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))
            end
            dfluxes = begin
                @stateflux snowpack ~ snowfall - melt
                @stateflux snowpack2 ~ snowfall - melt
            end
        end
        
        # Test bucket construction with direct API
        fluxes = [
            @hydroflux begin
                snowfall ~ step_func(Tmin - temp) * prcp
                rainfall ~ step_func(temp - Tmin) * prcp
            end
            @hydroflux melt ~ step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))
        ]
        
        dfluxes = [
            @stateflux snowpack ~ snowfall - melt
            @stateflux snowpack2 ~ snowfall - melt
        ]
        
        bucket_direct = NeuralBucket(fluxes=fluxes, dfluxes=dfluxes)
        
        # Test that both buckets have the same structure
        @test length(get_input_names(bucket_macro)) == length(get_input_names(bucket_direct))
        @test length(get_state_names(bucket_macro)) == length(get_state_names(bucket_direct))
        @test length(get_param_names(bucket_macro)) == length(get_param_names(bucket_direct))
        
        # Check that the expected names are present
        @test :temp ∈ get_input_names(bucket_macro)
        @test :prcp ∈ get_input_names(bucket_macro)
        @test :snowpack ∈ get_state_names(bucket_macro)
        @test :snowpack2 ∈ get_state_names(bucket_macro)
    end
    
    @testset "NeuralBucket Simulation" begin
        # Create a bucket for simulation
        bucket = @neuralbucket :test_bucket begin
            fluxes = begin
                @hydroflux begin
                    snowfall ~ step_func(Tmin - temp) * prcp
                    rainfall ~ step_func(temp - Tmin) * prcp
                end
                @hydroflux melt ~ step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))
            end
            dfluxes = begin
                @stateflux snowpack ~ snowfall - melt
                @stateflux snowpack2 ~ snowfall - melt
            end
        end
        
        # Set parameter values
        f_, Smax_, Qmax_, Df_, Tmax_, Tmin_ = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
        
        
        # Test multi-node simulation
        @testset "Multi-Node Simulation" begin
            # Number of nodes
            node_num = 5
            
            # Create parameters and initial states for multiple nodes
            node_params = ComponentVector(
                f=fill(f_, node_num), 
                Smax=fill(Smax_, node_num), 
                Qmax=fill(Qmax_, node_num),
                Df=fill(Df_, node_num), 
                Tmax=fill(Tmax_, node_num), 
                Tmin=fill(Tmin_, node_num)
            )
            init_states = ComponentVector(
                snowpack=fill(100.0, node_num), 
                snowpack2=fill(100.0, node_num)
            )
            pas = ComponentVector(params=node_params)
            
            # Create synthetic input data
            rng = StableRNG(456)
            ts_length = 30
            
            # Create base input data
            input_data = (
                lday = rand(rng, ts_length),
                temp = rand(rng, ts_length) .* 20 .- 5,
                prcp = rand(rng, ts_length) .* 10
            )
            
            # Convert to array format for single node
            input_names = get_input_names(bucket)
            input_mat = reduce(hcat, [input_data[name] for name in input_names]) |> permutedims
            
            # Expand to multi-node format (variables, nodes, timesteps)
            input_arr = permutedims(repeat(input_mat, 1, 1, node_num), (1, 3, 2))
            
            # Run simulation
            result = bucket(input_arr, pas; initstates=init_states)
            
            # Test that result has the expected shape
            @test size(result, 1) == length(get_output_names(bucket)) + length(get_state_names(bucket))
            @test size(result, 2) == node_num
            @test size(result, 3) == ts_length
        end
    end
end
