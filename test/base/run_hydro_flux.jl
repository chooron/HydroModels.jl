@testset "test single hydro flux (build by Symbolics.jl)" begin
    @variables a b c
    @parameters p1 p2
    simple_flux_1 = @hydroflux c ~ a * p1 + b * p2
    @test Set(HydroModels.get_input_names(simple_flux_1)) == Set([:a, :b])
    @test Set(HydroModels.get_param_names(simple_flux_1)) == Set([:p1, :p2])
    @test Set(HydroModels.get_output_names(simple_flux_1)) == Set([:c])
    output_mat = [18.0 17.0 11.0]
    @test simple_flux_1([2.0 3.0 1.0; 3.0 2.0 2.0], ComponentVector(params=(p1=3.0, p2=4.0))) == output_mat
end

@testset "test single hydro flux with multiple output" begin
    @variables a b c d
    @parameters p1 p2
    simple_flux_1 = @hydroflux begin
        c ~ a * p1 + b * p2
        d ~ a + p1 + b + p2
    end
    @test Set(HydroModels.get_input_names(simple_flux_1)) == Set([:a, :b])
    @test Set(HydroModels.get_param_names(simple_flux_1)) == Set([:p1, :p2])
    @test Set(HydroModels.get_output_names(simple_flux_1)) == Set([:c, :d])
    output_mat = [18.0 12.0; 17.0 12.0; 16.0 12.0]'
    @test simple_flux_1([2.0 3.0; 3.0 2.0; 4.0 1.0]', ComponentVector(params=(p1=3.0, p2=4.0))) == output_mat
end

@testset "test multi hydro flux (build by Symbolics.jl)" begin
    @variables a b c
    @parameters p1 p2
    simple_flux_1 = @hydroflux begin
        c ~ a * p1 + b * p2
        hru_types = collect(1:10)
    end
    @test Set(HydroModels.get_input_names(simple_flux_1)) == Set([:a, :b])
    @test Set(HydroModels.get_param_names(simple_flux_1)) == Set([:p1, :p2])
    @test Set(HydroModels.get_output_names(simple_flux_1)) == Set([:c])
    output_mat = [18.0 17.0 11.0]
    @test simple_flux_1([2.0 3.0 1.0; 3.0 2.0 2.0], ComponentVector(params=(p1=3.0, p2=4.0))) == output_mat
    # test with multiple nodes
    input_arr = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), [2.0 3.0 1.0; 3.0 2.0 2.0] for _ in 1:10), (1, 3, 2))
    input_pas = ComponentVector(params=(p1=fill(3.0, 10), p2=fill(4.0, 10)))
    @test simple_flux_1(input_arr, input_pas) == permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), output_mat for _ in 1:10), (1, 3, 2))
end

@testset "test multi hydro flux with multiple output" begin
    @variables a b c d
    @parameters p1 p2
    simple_flux_1 = @hydroflux begin
        c ~ a * p1 + b * p2
        d ~ a + p1 + b + p2
        hru_types = collect(1:10)
    end
    @test Set(HydroModels.get_input_names(simple_flux_1)) == Set([:a, :b])
    @test Set(HydroModels.get_param_names(simple_flux_1)) == Set([:p1, :p2])
    @test Set(HydroModels.get_output_names(simple_flux_1)) == Set([:c, :d])
    output_mat = [18.0 12.0; 17.0 12.0; 16.0 12.0]'
    @test simple_flux_1([2.0 3.0; 3.0 2.0; 4.0 1.0]', ComponentVector(params=(p1=3.0, p2=4.0))) == output_mat
    input_arr = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), [2.0 3.0; 3.0 2.0; 4.0 1.0]' for _ in 1:10), (1, 3, 2))
    input_pas = ComponentVector(params=(p1=fill(3.0, 10), p2=fill(4.0, 10)))
    @test simple_flux_1(input_arr, input_pas) == permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), output_mat for _ in 1:10), (1, 3, 2))
end

@testset "test state flux (build by Symbolics.jl)" begin
    # Define variables and parameters
    @variables a b c d e
    @parameters p1 p2

    # Test StateFlux with input and output, but no parameters
    state_flux_1 = @stateflux e ~ (a + b) - (c + d)
    @test Set(HydroModels.get_input_names(state_flux_1)) == Set([:a, :b, :c, :d])
    @test Set(HydroModels.get_param_names(state_flux_1)) == Set(Symbol[])
    @test Set(HydroModels.get_output_names(state_flux_1)) == Set(Symbol[])
    @test Set(HydroModels.get_state_names(state_flux_1)) == Set([:e])

    # Test StateFlux with input, parameters, and a custom expression
    state_flux_2 = @stateflux e ~ (a * p1 + b * p2) - (c + d)
    @test Set(HydroModels.get_input_names(state_flux_2)) == Set([:a, :b, :c, :d])
    @test Set(HydroModels.get_param_names(state_flux_2)) == Set([:p1, :p2])
    @test Set(HydroModels.get_output_names(state_flux_2)) == Set(Symbol[])
    @test Set(HydroModels.get_state_names(state_flux_2)) == Set([:e])
end
