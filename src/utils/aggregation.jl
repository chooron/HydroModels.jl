"""
    Build the aggregation function from a given graph
"""
function build_aggr_func(graph::DiGraph)
    adjacency = adjacency_matrix(graph)'
    aggr_func = (outflow) -> adjacency * outflow
    return aggr_func
end

"""
    Build the aggregation function from a given flow direction matrix and positions
"""
function build_aggr_func(flwdir::AbstractMatrix, positions::AbstractVector)
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
    #* build the outflow projection function
    aggr_func = (outflow) -> grid_routing(outflow, positions, flwdir)
    return aggr_func
end
