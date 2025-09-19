@testset "test grid route based on discharge route flux" begin
    flwdir = [1 4 8; 1 4 4; 1 1 2]
    positions = [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]

    @variables q1 q1_routed s_river
    @parameters lag

    route = @hydroroute begin
        fluxes = begin
            @hydroflux q1_routed ~ s_river / (1 + lag)
        end
        dfluxes = begin
            @stateflux s_river ~ q1 - q1_routed
        end
        aggr_func = HydroModels.build_aggr_func(flwdir, positions)
        hru_types = [1, 1, 1, 2, 2, 2, 3, 3, 3]
    end

    @test HydroModels.get_input_names(route) == [:q1]
    @test HydroModels.get_output_names(route) == [:q1_routed]
    @test HydroModels.get_param_names(route) == [:lag]
    @test HydroModels.get_state_names(route) == [:s_river]

    input_arr = ones(1, 9, 20)
    timeidx = collect(1:20)
    params = ComponentVector(lag=fill(0.2, 3))
    initstates = ComponentVector(s_river=fill(0.0, length(positions)))
    pas = ComponentVector(params=params)
    config = (solver=ManualSolver(mutable=true), interp=LinearInterpolation)
    output_arr = route(input_arr, pas, config, initstates=initstates, timeidx=timeidx)
    @test size(output_arr) == size(ones(2, 9, 20))
end

@testset "test vector route based on discharge route flux" begin
    @variables q1 q1_routed s_river q_gen
    @parameters lag

    network = DiGraph(9)
    add_edge!(network, 1, 2)
    add_edge!(network, 2, 5)
    add_edge!(network, 3, 5)
    add_edge!(network, 4, 5)
    add_edge!(network, 5, 8)
    add_edge!(network, 6, 9)
    add_edge!(network, 7, 8)
    add_edge!(network, 8, 9)

    params = ComponentVector(lag=fill(0.2, 3))
    initstates = ComponentVector(s_river=fill(0.0, 9))
    pas = ComponentVector(params=params)

    route = @hydroroute begin
        fluxes = begin
            @hydroflux q1_routed ~ s_river / (1 + lag)
        end
        dfluxes = begin
            @stateflux s_river ~ q1 - q1_routed
        end
        aggr_func = HydroModels.build_aggr_func(network)
        hru_types = [1, 1, 1, 2, 2, 2, 3, 3, 3]
    end

    input_arr = rand(1, 9, 20)
    timeidx = collect(1:20)
    config = (solver=ManualSolver(mutable=true), interp=LinearInterpolation)
    sol_2 = route(input_arr, pas, config, initstates=initstates, timeidx=timeidx)
    @test size(sol_2) == size(ones(2, 9, 20))
end

# @testset "test rapid route based on muskingum route flux" begin
#     @variables q q_routed

#     network = DiGraph(9)
#     add_edge!(network, 1, 2)
#     add_edge!(network, 2, 5)
#     add_edge!(network, 3, 5)
#     add_edge!(network, 4, 5)
#     add_edge!(network, 5, 8)
#     add_edge!(network, 6, 9)
#     add_edge!(network, 7, 8)
#     add_edge!(network, 8, 9)

#     ndtypes = [:ntype1, :ntype2, :ntype3]
#     params = ComponentVector(rapid_k=fill(2.0, length(ndtypes)), rapid_x=fill(0.2, length(ndtypes)))
#     ptyidx = [1, 1, 1, 2, 2, 2, 3, 3, 3]
#     pas = ComponentVector(params=params)
#     vroute = HydroModels.RapidRoute([q] => [q_routed], network=network)
#     input_arr = ones(1, 9, 20)
#     timeidx = collect(1:20)
#     sol_2 = vroute(input_arr, pas, timeidx=timeidx, ptyidx=ptyidx, solver=ManualSolver(mutable=true))
#     @test size(sol_2) == size(ones(1, 9, 20))
# end
