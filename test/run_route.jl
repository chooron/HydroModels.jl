@testset "test grid route based on discharge route flux" begin
    flwdir = [1 4 8; 1 4 4; 1 1 2]
    positions = [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]
    
    ndtypes = [:ntype1, :ntype2, :ntype3]
    hrunames = [:nid1, :nid2, :nid3, :nid4, :nid5, :nid6, :nid7, :nid8, :nid9]
    @variables q1 q1_routed s_river q_gen
    @parameters lag
    rflux = HydroModels.HydroFlux([q1, s_river] => [q_gen, q1_routed], [lag], exprs=[q1, s_river / (1 + lag) + q1])

    @test HydroModels.get_input_names(rflux) == [:q1, :s_river]
    @test HydroModels.get_output_names(rflux) == [:q_gen,:q1_routed]
    @test HydroModels.get_param_names(rflux) == [:lag]
    params = ComponentVector(lag=fill(0.2, length(ndtypes)))
    initstates = ComponentVector(s_river=fill(0.0, length(hrunames)))
    pas = ComponentVector(; params, initstates)
    route = HydroModels.GridRoute(rfunc=rflux, rstates=[s_river], flwdir=flwdir, positions=positions)
    @test HydroModels.get_input_names(route) == [:q1]
    @test HydroModels.get_output_names(route) == [:q_gen,:q1_routed]
    @test HydroModels.get_param_names(route) == [:lag]
    @test HydroModels.get_state_names(route) == [:s_river]
    
    ptyidx = [1, 1, 1, 2, 2, 2, 3, 3, 3]
    styidx = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    input_arr = ones(1, 9, 20)
    timeidx = collect(1:20)

    config = (solver=ManualSolver{true}(), interp=LinearInterpolation, ptyidx=ptyidx, styidx=styidx, timeidx=timeidx)
    output_arr = route(input_arr, pas, config=config)
    @test size(output_arr) == size(ones(3, 9, 20))
end

@testset "test vector route based on discharge route flux" begin
    @variables q1 q1_routed s_river q_gen
    @parameters lag
    rflux = HydroModels.HydroFlux([q1, s_river] => [q_gen, q1_routed], [lag], exprs=[q1, s_river / (1 + lag) + q1])
    
    network = DiGraph(9)
    add_edge!(network, 1, 2)
    add_edge!(network, 2, 5)
    add_edge!(network, 3, 5)
    add_edge!(network, 4, 5)
    add_edge!(network, 5, 8)
    add_edge!(network, 6, 9)
    add_edge!(network, 7, 8)
    add_edge!(network, 8, 9)
    
    ndtypes = [:ntype1, :ntype2, :ntype3]
    hrunames = [:nid1, :nid2, :nid3, :nid4, :nid5, :nid6, :nid7, :nid8, :nid9]
    params = ComponentVector(lag=fill(0.2, length(ndtypes)))
    initstates = ComponentVector(s_river=fill(0.0, length(hrunames)))
    pas = ComponentVector(; params, initstates)
    vroute = HydroModels.VectorRoute(rfunc=rflux, rstates=[s_river], network=network)

    ptyidx = [1, 1, 1, 2, 2, 2, 3, 3, 3]
    styidx = [1, 2, 3, 4, 5, 6, 7, 8, 9]

    input_arr = rand(1, 9, 20)
    timeidx = collect(1:20)
    sol_2 = vroute(input_arr, pas, config=(timeidx=timeidx, ptyidx=ptyidx, styidx=styidx, solver=HydroModels.ManualSolver{true}()))
    @test size(sol_2) == size(ones(3, 9, 20))
end

@testset "test rapid route based on muskingum route flux" begin
    @variables q q_routed

    network = DiGraph(9)
    add_edge!(network, 1, 2)
    add_edge!(network, 2, 5)
    add_edge!(network, 3, 5)
    add_edge!(network, 4, 5)
    add_edge!(network, 5, 8)
    add_edge!(network, 6, 9)
    add_edge!(network, 7, 8)
    add_edge!(network, 8, 9)

    ndtypes = [:ntype1, :ntype2, :ntype3]
    params = ComponentVector(rapid_k=fill(2.0, length(ndtypes)), rapid_x=fill(0.2, length(ndtypes)))
    ptyidx = [1, 1, 1, 2, 2, 2, 3, 3, 3]
    pas = ComponentVector(; params)
    vroute = HydroModels.RapidRoute([q]=>[q_routed], network=network)
    input_arr = ones(1, 9, 20)
    timeidx = collect(1:20)
    sol_2 = vroute(input_arr, pas, config=(timeidx=timeidx, ptyidx=ptyidx, solver=ManualSolver{true}()))
    @test size(sol_2) == size(ones(1, 9, 20))
end
