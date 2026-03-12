@testset "ChannelRoute constructor" begin
    single_calls = Ref(0)
    network_calls = Ref(0)

    simulate_single = function (route, inflow, params, delta_t)
        single_calls[] += 1
        @test route.algorithm == :test_channel
        @test route.name == :test_channel
        @test inflow == Float64[1.0, 2.0, 3.0]
        @test params == [2.5, 9.0]
        @test delta_t == 2.0
        return inflow .* params[1]
    end

    simulate_network = function (route, input, params, delta_t)
        network_calls[] += 1
        @test route.algorithm == :test_channel
        @test route.topo_order == [1, 2, 3]
        @test size(input) == (3, 4)
        @test params == [[10.0, 20.0, 10.0], [1.0, 2.0, 1.0]]
        @test delta_t == 3.0
        return input .* reshape(params[1], :, 1) .+ reshape(params[2], :, 1)
    end

    route = ChannelRoute(
        simulate_single,
        simulate_network;
        algorithm=:test_channel,
        inputs=[:qin],
        outputs=[:qout],
        params=[:lag, :bias],
        adjacency=[0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 1.0 0.0],
        max_lag=7,
        htypes=[1, 2, 1],
    )

    @test route isa ChannelRoute
    @test route.algorithm == :test_channel
    @test route.name == :test_channel
    @test route.max_lag == 7
    @test route.adjacency == [0.0 0.0 0.0; 1.0 0.0 0.0; 0.0 1.0 0.0]
    @test route.topo_order == [1, 2, 3]
    @test HydroModels.get_input_names(route) == [:qin]
    @test HydroModels.get_output_names(route) == [:qout]
    @test HydroModels.get_param_names(route) == [:lag, :bias]

    single_input = reshape(Float32[1, 2, 3], 1, :)
    single_params = ComponentVector(params=(lag=2.5, bias=9.0))
    single_output = route(single_input, single_params; delta_t=2.0)
    @test single_calls[] == 1
    @test network_calls[] == 0
    @test single_output == reshape(Float32[2.5, 5.0, 7.5], 1, :)

    network_input = zeros(Float32, 1, 3, 4)
    network_input[1, :, :] .= Float32[1 2 3 4; 5 6 7 8; 9 10 11 12]
    network_params = ComponentVector(params=(lag=[10.0, 20.0], bias=[1.0, 2.0]))
    network_output = route(network_input, network_params; delta_t=3.0)
    @test single_calls[] == 1
    @test network_calls[] == 1
    @test size(network_output) == (1, 3, 4)
    @test network_output[1, 1, :] == Float32[11, 21, 31, 41]
    @test network_output[1, 2, :] == Float32[102, 122, 142, 162]
    @test network_output[1, 3, :] == Float32[91, 101, 111, 121]
end

@testset "ChannelRoute validation" begin
    noop_single = (route, inflow, params, delta_t) -> Float64.(inflow)
    noop_network = (route, input, params, delta_t) -> Float64.(input)

    @test_throws ArgumentError ChannelRoute(
        noop_single,
        noop_network;
        algorithm=:bad_io,
        inputs=[:q1, :q2],
        outputs=[:qout],
    )

    @test_throws ArgumentError ChannelRoute(
        noop_single,
        noop_network;
        algorithm=:cyclic,
        inputs=[:qin],
        outputs=[:qout],
        adjacency=[0.0 1.0; 1.0 0.0],
    )

    route = ChannelRoute(
        noop_single,
        noop_network;
        algorithm=:validator,
        inputs=[:qin],
        outputs=[:qout],
        params=[:lag],
    )

    @test route.topo_order === nothing
    @test route.adjacency === nothing

    bad_single_input = zeros(Float64, 2, 4)
    @test_throws ArgumentError route(bad_single_input, ComponentVector(params=(lag=1.0,)))

    bad_network_input = zeros(Float64, 2, 3, 4)
    @test_throws ArgumentError route(bad_network_input, ComponentVector(params=(lag=[1.0, 2.0, 3.0],)))

    bad_scalar_params = ComponentVector(params=(lag=[1.0, 2.0],))
    @test_throws ArgumentError route(reshape(Float64[1, 2, 3], 1, :), bad_scalar_params)

    route_htypes = ChannelRoute(
        noop_single,
        noop_network;
        algorithm=:typed,
        inputs=[:qin],
        outputs=[:qout],
        params=[:lag],
        htypes=Int[],
    )
    @test route_htypes.htypes === nothing
end


@testset "ChannelRoute macro" begin
    macro_single = (route, inflow, params, delta_t) -> Float64.(inflow) .* params[1]
    macro_network = (route, input, params, delta_t) -> Float64.(input) .* reshape(params[1], :, 1)

    @variables qin qout
    @parameters lag bias

    named_route = @channelroute :macro_route begin
        algorithm = :macro_algorithm
        simulate_single = macro_single
        simulate_network = macro_network
        inputs = qin
        outputs = qout
        params = begin
            lag
            bias
        end
        adjacency = [0.0 0.0; 1.0 0.0]
        max_lag = 5
        htypes = [1, 2]
    end

    @test named_route isa ChannelRoute
    @test named_route.name == :macro_route
    @test named_route.algorithm == :macro_algorithm
    @test named_route.topo_order == [1, 2]
    @test named_route.max_lag == 5
    @test HydroModels.get_param_names(named_route) == [:lag, :bias]

    block_named_route = @channelroute begin
        name = :block_named_route
        algorithm = :block_algorithm
        simulate_single = macro_single
        simulate_network = macro_network
        inputs = begin
            qin
        end
        outputs = begin
            qout
        end
        params = lag
    end

    @test block_named_route.name == :block_named_route
    @test block_named_route.algorithm == :block_algorithm

    @test_throws Exception @eval @channelroute begin
        algorithm = :invalid_macro
        simulate_single = macro_single
        inputs = [:qin]
        outputs = [:qout]
    end
end

