@testset "test unit hydro flux" begin
    # Define the variables and parameters
    @variables q1 q1_lag t
    @parameters lag
    # Create a UnitHydroRouteFlux object
    # Input: q1 (flow)
    # Parameter: lag (routing parameter)
    # Using uh_1_half as the unit hydrograph function
    # Solve type: unithydro1 (convolution method)
    uh1 = HydroModels.@unithydro :maxbas_uh begin
        uh_vars = q1 => q1_lag
        uh_func = begin
            lag => (t / lag)^2.5
        end
    end

    # Test the input names of the router
    @test HydroModels.get_input_names(uh1) == [:q1]
    # Test the parameter names of the router
    @test HydroModels.get_param_names(uh1) == [:lag]
    # Test the output names of the router
    @test HydroModels.get_output_names(uh1) == [:q1_lag]
    # Test the routing function with sample input
    input_flow = Float32[2 3 4 2 3 1]
    params = ComponentVector(params=(lag=3.5,))
    expected_output = [0.0899066 0.643448 2.3442 3.20934 3.44646 2.20934]
    # [0.08726897695099571 0.5373023715895905 1.6508571480656808 2.839759323622619 3.2301609643779736 2.7991762465729138]
    # @test uh1(input_flow, params) ≈ expected_output atol = 1e-3
    input_arr = repeat(reshape(input_flow, 1, 1, length(input_flow)), 1, 10, 1)
    ndtypes = [Symbol("node_$i") for i in 1:10]
    input_pas = ComponentVector(params=(lag=fill(3.5, 10),),)
    config = Dict(:ptypes => ndtypes, :solver => ManualSolver(mutable=true))
    expected_output_arr = repeat(reshape(expected_output, 1, 1, length(expected_output)), 1, 10, 1)
    # @test uh1(input_arr, input_pas; config...) ≈ expected_output_arr atol = 1e-3
end

# @testset "test multi unit hydro flux" begin
#     # Define the variables and parameters
#     @variables q1 q1_lag t
#     @parameters lag
#     # Create a UnitHydroRouteFlux object
#     # Input: q1 (flow)
#     # Parameter: lag (routing parameter)
#     # Using uh_1_half as the unit hydrograph function
#     # Solve type: unithydro1 (convolution method)
#     uh1 = HydroModels.@unithydro :maxbas_uh begin
#         uh_vars = q1 => q1_lag
#         uh_func = begin
#             lag => (t / lag)^2.5
#         end
#         hru_types = collect(1:10)
#     end

#     # Test the input names of the router
#     @test HydroModels.get_input_names(uh1) == [:q1]
#     # Test the parameter names of the router
#     @test HydroModels.get_param_names(uh1) == [:lag]
#     # Test the output names of the router
#     @test HydroModels.get_output_names(uh1) == [:q1_lag]
#     # Test the routing function with sample input
#     input_flow = Float32[2 3 4 2 3 1]
#     params = ComponentVector(params=(lag=3.5,))
#     expected_output = [0.0899066 0.643448 2.3442 3.20934 3.44646 2.20934]
#     # [0.08726897695099571 0.5373023715895905 1.6508571480656808 2.839759323622619 3.2301609643779736 2.7991762465729138]
#     @test uh1(input_flow, params) ≈ expected_output atol = 1e-3
#     input_arr = repeat(reshape(input_flow, 1, 1, length(input_flow)), 1, 10, 1)
#     input_pas = ComponentVector(params=(lag=fill(3.5, 10),),)
#     expected_output_arr = repeat(reshape(expected_output, 1, 1, length(expected_output)), 1, 10, 1)
#     @test uh1(input_arr, input_pas) ≈ expected_output_arr atol = 1e-3
# end