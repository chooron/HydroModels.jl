@testset "d8 grid routing" begin
    d8_codes = [1, 2, 4, 8, 16, 32, 64, 128]
    d8_nn_pads = [(1, 1, 2, 0), (2, 0, 2, 0), (2, 0, 1, 1), (2, 0, 0, 2), (1, 1, 0, 2), (0, 2, 0, 2), (0, 2, 1, 1), (0, 2, 2, 0)]

    #* input dims: node_num * ts_len
    function grid_routing(input::AbstractVector, positions::AbstractVector, flwdir::AbstractMatrix)
        #* Convert input to sparse matrix
        input_arr = Array(sparse([pos[1] for pos in positions], [pos[2] for pos in positions], input, size(flwdir)[1], size(flwdir)[2]))
        #* Calculate weighted summation
        input_routed = sum(collect([pad_zeros(input_arr .* (flwdir .== code), arg) for (code, arg) in zip(d8_codes, d8_nn_pads)]))
        #* Clip input matrix border
        clip_arr = input_routed[2:size(input_arr)[1]+1, 2:size(input_arr)[2]+1]
        #* Convert input matrix to vector
        collect([clip_arr[pos[1], pos[2]] for pos in positions])
    end

    @testset "d8 grid routing (regular grid)" begin
        input = ones(9)
        positions = [(1, 1), (1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3)]
        flwdir = [1 4 8; 1 4 4; 1 1 2]
        result = grid_routing(input, positions, flwdir)
        target = [0.0, 1.0, 0.0, 0.0, 3.0, 0.0, 0.0, 2.0, 2.0]
        @test result == target
    end

    @testset "d8 grid routing (irregular grid)" begin
        input = ones(10)
        positions = [(1, 2), (1, 3), (2, 1), (2, 2), (2, 3), (2, 4), (3, 2), (3, 3), (4, 3), (4, 4)]
        flwdir = [0 4 4 0; 1 2 4 8; 0 1 4 0; 0 0 1 1]
        result = grid_routing(input, positions, flwdir)
        target = [0.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0, 4.0, 1.0, 1.0]
        @test result == target
    end
end
