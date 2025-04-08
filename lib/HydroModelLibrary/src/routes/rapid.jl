"""
    VectorRoute <: AbstractVectorRoute

A structure representing a vector-based routing scheme for hydrological modeling.

# Fields
- `rfunc::AbstractVector{<:AbstractRouteFlux}`: Vector of routing flux functions for each node.
- `network::DiGraph`: A directed graph representing the routing network topology.
- `infos::NamedTuple`: Contains information about the VectorRoute instance, including input, output, state, and parameter names.

# Constructor
    VectorRoute(
        name::Symbol;
        rfunc::AbstractRouteFlux,
        network::DiGraph
    )

Constructs a `VectorRoute` object with the given name, routing flux function, and network structure.

# Arguments
- `name::Symbol`: A symbol representing the name of the routing scheme.
- `rfunc::AbstractRouteFlux`: The routing flux function to be applied at each node.
- `network::DiGraph`: A directed graph representing the routing network topology.

The constructor extracts variable names, parameter names, and neural network names from the provided
routing flux function to set up the internal information structure of the `VectorRoute` object.

Note: 来源于Rapid汇流模型
"""
struct RapidRoute <: AbstractRapidRoute
    "Routing adjacency matrix"
    adjacency::AbstractMatrix
    "node index"
    nodeids::Vector{Symbol}
    "Metadata: contains keys for input, output, param, state, and nn"
    meta::HydroMeta

    function RapidRoute(;
        input::Num, # 计算qout的函数
        output::Num,
        adjacency::AbstractMatrix,
        nodeids::Vector{Symbol},
        name::Union{Symbol,Nothing}=nothing,
    )
        input_name = Symbolics.tosymbol(input)
        output_name = Symbolics.tosymbol(output)
        #* Setup the name information of the hydrobucket
        route_name = name === nothing ? Symbol(Symbol(reduce((x, y) -> Symbol(x, y), input_name)), :_rapid_route) : name
        meta = HydroMeta(name=route_name, inputs=[input_name], outputs=[output_name], states=[], params=[:k, :x], nns=[])

        return new(
            adjacency,
            nodeids,
            meta,
        )
    end
end

function (route::RapidRoute)(
    input::Array,
    pas::ComponentVector;
    config::NamedTuple=NamedTuple(),
    kwargs...
)
    ptypes = get(config, :ptypes, collect(keys(pas[:params])))
    interp = get(config, :interp, LinearInterpolation)
    timeidx = get(config, :timeidx, collect(1:size(input, 3)))
    delta_t = get(config, :delta_t, 1.0)
    solver = DiscreteSolver(alg=FunctionMap{true}())

    @assert all(ptype in keys(pas[:params]) for ptype in ptypes) "Missing required parameters. Expected all of $(keys(pas[:params])), but got $(ptypes)."

    #* var num * node num * ts len
    #* 计算出每个node结果的插值函数
    itp_funcs = interp.(eachslice(input[1, :, :], dims=1), Ref(timeidx), extrapolate=true)

    #* prepare the parameters for the routing function
    #* 参数可能存在转换需求,其中包括A通常是固定值
    k_ps = [pas[:params][ptype][:k] for ptype in ptypes]
    x_ps = [pas[:params][ptype][:x] for ptype in ptypes]
    c0 = @. ((delta_t / k_ps) - (2 * x_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    c1 = @. ((delta_t / k_ps) + (2 * x_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    c2 = @. ((2 * (1 - x_ps)) - (delta_t / k_ps)) / ((2 * (1 - x_ps)) + (delta_t / k_ps))
    A = (p) -> Matrix(I, size(route.adjacency)...) .- diagm(p.c0) * route.adjacency

    function route_ode!(du, u, p, t)
        q_out_t1 = u
        q_gen = [itp_func(t) for itp_func in itp_funcs]
        #* Ax = b, x is the q_out(t+1)
        q_in_t1 = route.adjacency * q_out_t1
        rflux_b = @. p.c0 * q_gen + p.c1 * (q_in_t1 + q_gen) + p.c2 * q_out_t1
        #* solve the linear equation (simple solve by matrix inversion)
        du[:] = inv(A(p)) * (rflux_b .- A(p) * q_out_t1)
    end

    #* solve the ode
    sol_arr = solver(route_ode!, ComponentVector(c0=c0, c1=c1, c2=c2), zeros(size(input)[2]), timeidx, convert_to_array=true)
    return reshape(sol_arr, 1, size(sol_arr)...)
end
