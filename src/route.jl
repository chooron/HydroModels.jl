"""
    WeightSumRoute <: AbstractRouteFlux

Represents a weighted cumulative sum routing structure for hydrological modeling.

# Fields
- `infos::NamedTuple`: A named tuple containing routing information with keys:
  - `name::Symbol`: The name of the routing component.
  - `input::Symbol`: The symbol representing the input variable.
  - `output::Symbol`: The symbol representing the output variable.
  - `param::Symbol`: The symbol representing the weight parameter.

# Constructor
    WeightSumRoute(
        name::Symbol;
        input::Num,
        output::Num,
        param::Num
    )

Constructs a WeightSumRoute instance.

# Arguments
- `name::Symbol`: The name of the routing component.
- `input::Num`: The input variable.
- `output::Num`: The output variable.
- `param::Num`: The weight parameter.

# Description
WeightSumRoute applies a weighted cumulative sum operation to the input,
where the weights are specified by the `param` parameter. This routing method
is useful for scenarios where the contribution of each input node needs to be
weighted differently in the cumulative output.
"""
struct DirectRoute <: AbstractDirectRoute
    rfunc::AbstractFlux
    "grid subarea information, km2"
    subareas::AbstractVector
    "Metadata: contains keys for input, output, param, state, and nn"
    meta::HydroMeta
    # "output identifier" 
    #=     
    This part doesn't need to add an output id, because the matrix cannot be directly modified,
    so this part is all one value, representing the output result at the output id
    =#
    # outid::Symbol

    function DirectRoute(;
        rfunc::AbstractFlux,
        subareas::AbstractVector,
        name::Union{Symbol,Nothing}=nothing
    )
        #* Extract all variable names of funcs and dfuncs
        input_names, output_names, state_names = get_var_names(rfunc)
        #* Extract all parameters names of funcs and dfuncs
        param_names = get_param_names(rfunc)
        push!(param_names, :route_weight)
        #* Extract all neuralnetwork names of the funcs
        nn_names = get_nn_names(rfunc)
        #* Setup the name information of the hydrobucket
        route_name = name === nothing ? Symbol(Symbol(reduce((x, y) -> Symbol(x, y), output_names)), :_weight_route) : name
        meta = HydroMeta(name=route_name, inputs=input_names, outputs=output_names, states=state_names, params=param_names, nns=nn_names)
        return new(
            rfunc,
            subareas,
            meta,
        )
    end
end

function (route::DirectRoute)(
    input::AbstractArray,
    pas::ComponentVector;
    config::NamedTuple=NamedTuple(),
    kwargs...
)
    ptypes = get(config, :ptypes, collect(keys(pas[:params])))
    #* 计算出每个节点的面积转换系数
    area_coefs = @. 24 * 3600 / (route.subareas * 1e6) * 1e3

    rfunc_output = route.rfunc(input, pas, ptypes=ptypes)
    weight_params = [pas[:params][ptype][:route_weight] for ptype in ptypes]
    weight_result = sum(rfunc_output[1, :, :] .* weight_params .* area_coefs, dims=1)
    # expand dims
    output_arr = repeat(weight_result, outer=(size(input)[1], 1))
    reshape(output_arr, 1, size(output_arr)...)
end

"""
    GridRoute(name::Symbol; rfunc::AbstractRouteFlux, flwdir::AbstractMatrix, positions::AbstractVector)

Represents a grid-based routing structure for hydrological modeling.

# Arguments
- `name::Symbol`: A symbol representing the name of the GridRoute instance.
- `rfunc::AbstractRouteFlux`: The routing function used for flow calculations.
- `flwdir::AbstractMatrix`: A matrix representing the flow direction for each grid cell.
- `positions::AbstractVector`: A vector of positions for each node in the grid.

# Fields
- `rfunc::AbstractRouteFlux`: The routing function used for flow calculations.
- `flwdir::AbstractMatrix`: A matrix representing the flow direction for each grid cell.
- `positions::AbstractVector`: A vector of positions for each node in the grid.
- `infos::NamedTuple`: Contains information about the GridRoute instance, including input, output, state, and parameter names.

# Description
GridRoute is a structure that represents a grid-based routing system in a hydrological model. 
It uses a specified routing function (`rfunc`) to calculate flow between grid cells based on 
the provided flow direction matrix (`flwdir`) and node positions (`positions`).

The `infos` field stores metadata about the GridRoute instance, including names of inputs, 
outputs, states, parameters, and neural networks (if applicable) derived from the routing function.

This structure is particularly useful for modeling water flow across a landscape represented as a grid, 
where each cell has a specific flow direction and contributes to downstream flow.
"""
struct HydroRoute <: AbstractHydroRoute
    "Routing function"
    rfunc::AbstractHydroFlux
    "Outflow projection function"
    projfunc::Function
    "grid subarea information, km2"
    subareas::AbstractVector
    "node index"
    nodeids::Vector{Symbol}
    "Metadata: contains keys for input, output, param, state, and nn"
    meta::HydroMeta

    function HydroRoute(;
        rfunc::AbstractHydroFlux,
        rstate::Num,
        projfunc::Function,
        subareas::Union{AbstractVector,Number},
        nodeids::Vector{Symbol},
        name::Union{Symbol,Nothing}=nothing,
    )
        #* Extract all variable names of funcs and dfuncs
        input_names, output_names = get_var_names(rfunc)
        state_name = Symbolics.tosymbol(rstate)
        input_names = setdiff(input_names, output_names)
        #* Extract all parameters names of funcs and dfuncs
        param_names = get_param_names(rfunc)
        #* Extract all neuralnetwork names of the funcs
        nn_names = get_nn_names(rfunc)
        #* Setup the name information of the hydrobucket
        route_name = name === nothing ? Symbol(Symbol(reduce((x, y) -> Symbol(x, y), input_names)), :_grid_route) : name
        meta = HydroMeta(route_name, input_names, output_names, param_names, [state_name], nn_names)
        #* fill the subareas vector if it's a single value
        subareas = subareas isa AbstractVector ? subareas : fill(subareas, length(nodeids))

        return new(
            rfunc,
            projfunc,
            subareas,
            nodeids,
            meta,
        )
    end
end

function GridRoute(;
    rfunc::AbstractHydroFlux,
    rstate::Num,
    flwdir::AbstractMatrix,
    positions::AbstractVector,
    subareas::Union{AbstractVector,Number},
    nodeids::Vector{Symbol},
    name::Union{Symbol,Nothing}=nothing,
)
    d8_codes = [1, 2, 4, 8, 16, 32, 64, 128]
    d8_nn_pads = [(1, 1, 2, 0), (2, 0, 2, 0), (2, 0, 1, 1), (2, 0, 0, 2), (1, 1, 0, 2), (0, 2, 0, 2), (0, 2, 1, 1), (0, 2, 2, 0),]

    """
    input dims: node_num * ts_len
    """
    function grid_routing(input::AbstractVector, positions::AbstractVector, flwdir::AbstractMatrix)
        #* 转换为input的稀疏矩阵
        input_arr = Array(sparse([pos[1] for pos in positions], [pos[2] for pos in positions], input, size(flwdir)[1], size(flwdir)[2]))
        #* 计算权重求和结果
        input_routed = sum(collect([pad_zeros(input_arr .* (flwdir .== code), arg) for (code, arg) in zip(d8_codes, d8_nn_pads)]))
        #* 裁剪输入矩阵边框
        clip_arr = input_routed[2:size(input_arr)[1]+1, 2:size(input_arr)[2]+1]
        #* 将输入矩阵转换为向量
        collect([clip_arr[pos[1], pos[2]] for pos in positions])
    end
    @assert length(positions) == length(nodeids) "The length of positions must be the same as the length of nodeids, but got positions: $(length(positions)) and nodeids: $(length(nodeids))"
    #* build the outflow projection function
    projfunc = (outflow) -> grid_routing(outflow, positions, flwdir)
    return HydroRoute(; rfunc, rstate, projfunc, subareas, nodeids, name)
end

function VectorRoute(;
    rfunc::AbstractHydroFlux,
    rstate::Num,
    network::DiGraph,
    subareas::Union{AbstractVector,Number},
    nodeids::Vector{Symbol},
    name::Union{Symbol,Nothing}=nothing,
)
    @assert length(nodeids) == nv(network) "The length of nodeids must be the same as the number of nodes, but got nodeids: $(length(nodeids)) and nodes: $(nv(network))"
    #* generate adjacency matrix from network
    adjacency = adjacency_matrix(network)'
    #* build the outflow projection function
    projfunc = (outflow) -> adjacency * outflow
    return HydroRoute(; rfunc, rstate, projfunc, subareas, nodeids, name)
end

function (route::HydroRoute)(
    input::Array,
    pas::ComponentVector;
    config::NamedTuple=NamedTuple(),
    kwargs...
)
    #* get the parameter types and state types
    ptypes = get(config, :ptypes, collect(keys(pas[:params])))
    stypes = get(config, :stypes, collect(keys(pas[:initstates])))
    #* get the interpolation type and solver type
    interp = get(config, :interp, LinearInterpolation)
    solver = get(config, :solver, ODESolver())
    #* get the time index
    timeidx = get(config, :timeidx, collect(1:size(input, 3)))

    #* check the length of subareas, ptypes, stypes
    @assert length(route.subareas) == size(input, 2) "The length of subareas must be the same as the number of nodes, but got subareas: $(length(route.subareas)) and input: $(size(input, 2))"
    @assert length(ptypes) == length(stypes) == size(input, 2) "The length of ptypes and stypes must be the same as the number of nodes, but got ptypes is $(length(ptypes)) and stypes is $(length(stypes))"
    @assert all(ptype in keys(pas[:params]) for ptype in ptypes) "Missing required parameters. Expected all of $(keys(pas[:params])), but got ptypes is $(length(ptypes))."
    @assert all(stype in keys(pas[:initstates]) for stype in stypes) "Missing required initial states. Expected all of $(keys(pas[:initstates])), but got stypes is $(length(stypes))."

    #* Extract the idx range of each variable in params,
    #* this extraction method is significantly more efficient than extracting by name
    #* Construct a function for the ode_func input variable. Because of the difference in t, the ode_func input is not fixed.
    #* The input format is the input variable plus the intermediate state, which is consistent with the input of ode_func
    if route.rfunc isa AbstractNeuralFlux
        nn_params_idx = [getaxes(pas[:nn])[1][nm].idx for nm in get_nn_names(route.rfunc)]
        param_func = (p) -> Ref([p[:nn][idx] for idx in nn_params_idx])
    else
        rflux_params_idx = [getaxes(pas[:params][ptypes[1]])[1][nm].idx for nm in get_param_names(route.rfunc)]
        param_func = (p) -> [p[:params][ptype][rflux_params_idx] for ptype in ptypes]
    end

    #* solve the problem
    sol_arr = solve_prob(route, input, pas, param_func, timeidx=timeidx, stypes=stypes, solver=solver, interp=interp)
    sol_arr_permuted = permutedims(sol_arr, (2, 1, 3))
    cont_arr = cat(input, sol_arr_permuted, dims=1)
    output_vec = [route.rfunc.func.(eachslice(cont_arr[:, :, i], dims=2), param_func(pas), timeidx[i]) for i in 1:size(input)[3]]
    out_arr = reduce(hcat, reduce.(vcat, output_vec))
    #  return route_states and q_out
    return cat(sol_arr_permuted, reshape(out_arr, 1, size(out_arr)...), dims=1)
end

function solve_prob(
    route::HydroRoute,
    input::Array,
    pas::ComponentVector,
    paramfunc::Function;
    timeidx::Vector{<:Number}=collect(1:size(input, 3)),
    stypes::Vector{Symbol}=collect(keys(pas[:initstates])),
    solver::AbstractSolver=ODESolver(),
    interp::Type{<:AbstractInterpolation}=LinearInterpolation,
)
    #* Interpolate the input data. Since ordinary differential calculation is required, the data input must be continuous,
    #* so an interpolation function can be constructed to apply to each time point.
    itp_funcs = interp.(eachslice(input[1, :, :], dims=1), Ref(timeidx), extrapolate=true)

    #* 准备初始状态
    init_states_vec = collect([collect(pas[:initstates][stype][get_state_names(route)]) for stype in stypes])
    #* prepare the initial states matrix (dims: state_num * node_num)
    init_states_mat = reduce(hcat, init_states_vec)'

    #* Construct a temporary function that couples multiple ode functions to construct the solution for all states under the bucket
    function route_ode!(du, u, p, t)
        # 提取单元产流
        q_gen = [itp_func(t) for itp_func in itp_funcs]
        # 计算单元出流,更新单元出流状态
        q_out_vec = route.rfunc.func.(eachslice(hcat(q_gen, u), dims=1), paramfunc(p), Ref(t))
        q_out = reduce(vcat, q_out_vec)
        q_in = route.projfunc(q_out)
        # 更新状态
        du[:] = q_in .+ q_gen .- q_out
    end

    #* Solve the problem using the solver wrapper
    sol = solver(route_ode!, pas, init_states_mat, timeidx, convert_to_array=true)
    # return route_states and q_out
    sol
end
