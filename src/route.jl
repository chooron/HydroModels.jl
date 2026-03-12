"""
Route module - defines hydrological routing components, supporting both symbolic
and functional construction approaches, plus explicit IRF and sparse-kernel routing.
"""

"""
    ChannelRoute{SF, NF, AT, OT, HT, I} <: AbstractHydroRoute

Generic channel-routing component backed by algorithm-specific callables for
single-reach and network simulations.
"""
_channelroute_names(value::AbstractVector) = value
_channelroute_names(value) = [value]

struct ChannelRoute{SF,NF,AT,OT,HT,I} <: AbstractHydroRoute
    name::Symbol
    algorithm::Symbol
    simulate_single::SF
    simulate_network::NF
    adjacency::AT
    topo_order::OT
    max_lag::Int
    htypes::HT
    infos::I

    function ChannelRoute(
        simulate_single::SF,
        simulate_network::NF;
        algorithm::Symbol,
        inputs,
        outputs,
        params=Symbol[],
        adjacency::Union{Nothing,AbstractMatrix}=nothing,
        max_lag::Integer=32,
        htypes=nothing,
        name::Optional{Symbol}=nothing,
    ) where {SF,NF}
        input_names = _channelroute_names(inputs)
        output_names = _channelroute_names(outputs)
        param_names = _channelroute_names(params)

        length(input_names) == length(output_names) == 1 || throw(ArgumentError("ChannelRoute supports one input and one output."))

        htypes_norm = htypes isa Vector{Int} && isempty(htypes) ? nothing : htypes
        adjacency_mat = isnothing(adjacency) ? nothing : Matrix{Float64}(adjacency)
        topo_order = isnothing(adjacency_mat) ? nothing : _topological_order(adjacency_mat)
        infos = HydroInfos(
            inputs=tosymbol.(input_names),
            outputs=tosymbol.(output_names),
            params=!isempty(param_names) ? tosymbol.(param_names) : Symbol[],
        )
        route_name = isnothing(name) ? algorithm : name

        return new{SF,NF,typeof(adjacency_mat),typeof(topo_order),typeof(htypes_norm),typeof(infos)}(
            route_name,
            algorithm,
            simulate_single,
            simulate_network,
            adjacency_mat,
            topo_order,
            Int(max_lag),
            htypes_norm,
            infos,
        )
    end
end

macro channelroute(args...)
    name_arg = nothing
    expr = nothing

    if length(args) == 1
        expr = args[1]
    elseif length(args) == 2
        name_arg = args[1]
        expr = args[2]
    else
        error("@channelroute accepts either a begin...end block or `name begin ... end`.")
    end

    @assert Meta.isexpr(expr, :block) "Expected a begin...end block after channel route name"

    normalize_names(rhs) = Meta.isexpr(rhs, :block) ? Expr(:vect, filter(x -> !(x isa LineNumberNode), rhs.args)...) : rhs

    simulate_single_expr = nothing
    simulate_network_expr = nothing
    algorithm_expr = nothing
    inputs_expr = nothing
    outputs_expr = nothing
    params_expr = :(Symbol[])
    adjacency_expr = :(nothing)
    max_lag_expr = 32
    htypes_expr = :(nothing)
    block_name_expr = nothing

    for assign in filter(x -> !(x isa LineNumberNode), expr.args)
        @assert Meta.isexpr(assign, :(=)) "Expected assignments in the form `algorithm = ...`, `simulate_single = ...`, etc."
        lhs, rhs = assign.args

        if lhs == :simulate_single
            simulate_single_expr = rhs
        elseif lhs == :simulate_network
            simulate_network_expr = rhs
        elseif lhs == :algorithm
            algorithm_expr = rhs
        elseif lhs == :inputs
            inputs_expr = normalize_names(rhs)
        elseif lhs == :outputs
            outputs_expr = normalize_names(rhs)
        elseif lhs == :params
            params_expr = normalize_names(rhs)
        elseif lhs == :adjacency
            adjacency_expr = rhs
        elseif lhs == :max_lag
            max_lag_expr = rhs
        elseif lhs == :htypes
            htypes_expr = rhs
        elseif lhs == :name
            block_name_expr = rhs
        else
            error("Unknown assignment: $(lhs). Expected `algorithm`, `simulate_single`, `simulate_network`, `inputs`, `outputs`, `params`, `adjacency`, `max_lag`, `htypes`, or `name`.")
        end
    end

    isnothing(name_arg) || isnothing(block_name_expr) || error("Specify channel route name either as macro argument or as `name = ...`, not both.")
    name_expr = isnothing(block_name_expr) ? name_arg : block_name_expr

    err_msg = "`algorithm`, `simulate_single`, `simulate_network`, `inputs`, and `outputs` must all be specified"
    @assert !isnothing(simulate_single_expr) && !isnothing(simulate_network_expr) &&
            !isnothing(algorithm_expr) && !isnothing(inputs_expr) && !isnothing(outputs_expr) err_msg

    return esc(quote
        ChannelRoute(
            $simulate_single_expr,
            $simulate_network_expr;
            algorithm=$algorithm_expr,
            inputs=$inputs_expr,
            outputs=$outputs_expr,
            params=$params_expr,
            adjacency=$adjacency_expr,
            max_lag=$max_lag_expr,
            htypes=$htypes_expr,
            name=$(name_expr),
        )
    end)
end
function _topological_order(adjacency::AbstractMatrix)
    size(adjacency, 1) == size(adjacency, 2) || throw(ArgumentError("adjacency must be square."))
    nn = size(adjacency, 1)
    indegree = [count(!iszero, @view adjacency[row, :]) for row in 1:nn]
    pending = collect(findall(==(0), indegree))
    order = Int[]

    while !isempty(pending)
        node = popfirst!(pending)
        push!(order, node)
        for downstream in 1:nn
            if !iszero(adjacency[downstream, node])
                indegree[downstream] -= 1
                if indegree[downstream] == 0
                    push!(pending, downstream)
                end
            end
        end
    end

    length(order) == nn || throw(ArgumentError("adjacency must define an acyclic upstream-to-downstream network."))
    return order
end

function _scalar_channel_param(value, name::Symbol)
    if value isa Number
        return Float64(value)
    elseif value isa AbstractVector && length(value) == 1
        return Float64(first(value))
    end
    throw(ArgumentError("Parameter $name must be scalar for single-node routing."))
end

function _node_channel_param(value, name::Symbol, nn::Int, htypes)
    if !isnothing(htypes)
        if value isa Number
            return fill(Float64(value), length(htypes))
        end
        return Float64.(collect(value[htypes]))
    elseif value isa Number
        return fill(Float64(value), nn)
    elseif value isa AbstractVector && length(value) == nn
        return Float64.(collect(value))
    elseif value isa AbstractVector && length(value) == 1
        return fill(Float64(first(value)), nn)
    end
    throw(ArgumentError("Parameter $name must be scalar, length 1, or length $nn."))
end

function _extract_scalar_channel_params(route::ChannelRoute, params::ComponentVector)
    param_names = collect(get_param_names(route))
    return [_scalar_channel_param(params[:params][name], name) for name in param_names]
end

function _extract_node_channel_params(route::ChannelRoute, params::ComponentVector, nn::Int)
    param_names = collect(get_param_names(route))
    return [_node_channel_param(params[:params][name], name, nn, route.htypes) for name in param_names]
end

function (route::ChannelRoute)(
    input::AbstractArray{T,2},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
) where {T}
    size(input, 1) == 1 || throw(ArgumentError("ChannelRoute expects a single input variable."))
    params_cv = _as_componentvector(params)
    delta_t = Float64(get(kwargs, :delta_t, 1.0))
    routed = route.simulate_single(route, vec(input[1, :]), _extract_scalar_channel_params(route, params_cv), delta_t)
    return reshape(convert.(T, routed), 1, :)
end

function (route::ChannelRoute)(
    input::AbstractArray{T,3},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
) where {T}
    size(input, 1) == 1 || throw(ArgumentError("ChannelRoute expects a single input variable."))
    params_cv = _as_componentvector(params)
    delta_t = Float64(get(kwargs, :delta_t, 1.0))
    routed = route.simulate_network(route, input[1, :, :], _extract_node_channel_params(route, params_cv, size(input, 2)), delta_t)
    return reshape(convert.(T, routed), 1, size(routed, 1), size(routed, 2))
end
"""
    HydroRoute{FF, OF, AF, HT, I} <: AbstractHydroRoute

Represents a hydrological routing component that combines flux calculations with
an aggregation/distribution step.
"""
struct HydroRoute{FF,OF,AF,HT,I} <: AbstractHydroRoute
    name::Symbol
    flux_func::FF
    ode_func::OF
    aggr_func::AF
    htypes::HT
    infos::I

    function HydroRoute(;
        rfluxes::Vector{<:AbstractHydroFlux},
        dfluxes::Vector{<:AbstractStateFlux},
        htypes::Vector{Int},
        aggr_func::AF,
        name::Optional{Symbol}=nothing,
    ) where AF
        @assert !isempty(htypes) "htypes must be specified for HydroRoute"

        inputs, outputs, states = get_var_names(vcat(rfluxes, dfluxes))
        params = reduce(union, get_param_names.(vcat(rfluxes, dfluxes)); init=Symbol[])
        nns = reduce(union, get_nn_names.(rfluxes); init=Symbol[])

        infos = HydroInfos(
            inputs=inputs, states=states, outputs=outputs,
            params=params, nns=nns
        )
        route_name = isnothing(name) ? Symbol("##route#", hash(infos)) : name
        flux_func, ode_func = build_route_func(rfluxes, dfluxes, infos)

        return new{typeof(flux_func),typeof(ode_func),AF,typeof(htypes),typeof(infos)}(
            route_name, flux_func, ode_func, aggr_func, htypes, infos
        )
    end

    function HydroRoute(
        flux_func::Function,
        ode_func::Function,
        aggr_func::AF;
        name::Symbol,
        inputs::Vector{Symbol},
        outputs::Vector{Symbol},
        states::Vector{Symbol},
        params::Vector{Symbol}=Symbol[],
        htypes::Vector{Int}
    ) where AF
        @assert !isempty(htypes) "htypes must be specified for HydroRoute"

        infos = HydroInfos(
            inputs=inputs, states=states, outputs=outputs,
            params=params, nns=Symbol[]
        )

        return new{typeof(flux_func),typeof(ode_func),AF,typeof(htypes),typeof(infos)}(
            name, flux_func, ode_func, aggr_func, htypes, infos
        )
    end
end

macro hydroroute(args...)
    name = length(args) == 1 ? nothing : args[1]
    expr = length(args) == 1 ? args[1] : args[2]

    @assert Meta.isexpr(expr, :block) "Expected a begin...end block after route name"

    fluxes_expr, dfluxes_expr, aggr_func_expr, htypes_expr = nothing, nothing, nothing, nothing

    for assign in filter(x -> !(x isa LineNumberNode), expr.args)
        @assert Meta.isexpr(assign, :(=)) "Expected assignments in the form 'fluxes = begin...end', etc."
        lhs, rhs = assign.args

        if lhs == :fluxes
            @assert Meta.isexpr(rhs, :block) "Expected 'fluxes' to be defined in a begin...end block"
            fluxes_expr = Expr(:vect, filter(x -> !(x isa LineNumberNode), rhs.args)...)
        elseif lhs == :dfluxes
            @assert Meta.isexpr(rhs, :block) "Expected 'dfluxes' to be defined in a begin...end block"
            dfluxes_expr = Expr(:vect, filter(x -> !(x isa LineNumberNode), rhs.args)...)
        elseif lhs == :aggr_func
            aggr_func_expr = rhs
        elseif lhs == :htypes
            htypes_expr = rhs
        else
            error("Unknown assignment: $(lhs). Expected 'fluxes', 'dfluxes', 'htypes', or 'aggr_func'")
        end
    end

    err_msg = "'fluxes', 'dfluxes', 'htypes', and 'aggr_func' must all be specified"
    @assert !isnothing(fluxes_expr) && !isnothing(dfluxes_expr) &&
            !isnothing(aggr_func_expr) && !isnothing(htypes_expr) err_msg

    return esc(quote
        HydroRoute(
            rfluxes=$fluxes_expr,
            dfluxes=$dfluxes_expr,
            aggr_func=$aggr_func_expr,
            htypes=$htypes_expr,
            name=$(name)
        )
    end)
end

function build_aggr_func(graph::DiGraph)
    adjacency = adjacency_matrix(graph)'
    aggr_func = (outflow) -> adjacency * outflow
    return aggr_func
end

function build_aggr_func(flwdir::AbstractMatrix, positions::AbstractVector)
    d8_codes = [1, 2, 4, 8, 16, 32, 64, 128]
    d8_nn_pads = [
        (1, 1, 2, 0), (2, 0, 2, 0), (2, 0, 1, 1), (2, 0, 0, 2),
        (1, 1, 0, 2), (0, 2, 0, 2), (0, 2, 1, 1), (0, 2, 2, 0)
    ]

    cartesian_positions = [CartesianIndex(p...) for p in positions]

    function grid_routing(input::AbstractVector)
        rows, cols = size(flwdir)
        T = eltype(input)

        input_arr = zeros(T, rows, cols)
        input_arr[cartesian_positions] = input

        input_routed = sum(pad_zeros(input_arr .* (flwdir .== code), pad) for (code, pad) in zip(d8_codes, d8_nn_pads))

        clip_arr = input_routed[2:rows+1, 2:cols+1]
        return clip_arr[cartesian_positions]
    end

    return grid_routing
end

function (route::HydroRoute)(
    input::AbstractArray{T,3},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
) where {T}
    params = _as_componentvector(params)
    config_norm = normalize_config(config)
    solve_type = config_norm.solver
    interp_type = config_norm.interpolator

    input_dims, num_nodes, time_len = size(input)
    new_pas = expand_component_params(params, get_param_names(route), route.htypes)

    num_states = length(get_state_names(route))
    initstates = get(kwargs, :initstates, zeros(T, num_states * num_nodes))
    timeidx = isempty(config_norm.timeidx) ? collect(1:time_len) : config_norm.timeidx

    itpfunc = hydrointerp(interp_type, reshape(input, input_dims * num_nodes, time_len), timeidx)

    solved_states = hydrosolve(
        solve_type,
        (u, p, t) -> begin
            tmp_outflow, tmp_states = route.ode_func(
                reshape(itpfunc(t), input_dims, num_nodes),
                u,
                p
            )
            tmp_inflow_arr = route.aggr_func(tmp_outflow)
            tmp_states .+ tmp_inflow_arr
        end,
        new_pas, initstates, timeidx, config_norm
    )

    solved_states_reshape = reshape(solved_states, num_states, num_nodes, time_len)
    output = route.flux_func(input, solved_states_reshape, new_pas)
    cat(solved_states_reshape, stack(output, dims=1), dims=1)
end

function (route::HydroRoute)(
    input::AbstractArray{T,2},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
) where {T}
    error("HydroRoute only accepts 3D input (variables x nodes x time).\n" *
          "Routing is inherently spatial and requires multiple nodes.\n" *
          "Got input shape: $(size(input))")
end

struct RouteIRF{KF,NT}
    name::Symbol
    kernel_func::KF
    infos::NT

    function RouteIRF(
        params::AbstractVector,
        kernel_func::Function;
        name::Optional{Symbol}=nothing,
    )
        infos = HydroInfos(params=!isempty(params) ? tosymbol.(params) : Symbol[])
        irf_name = isnothing(name) ? Symbol("##route_irf#", hash(infos)) : name
        return new{typeof(kernel_func),typeof(infos)}(irf_name, kernel_func, infos)
    end
end

function (irf::RouteIRF)(params::AbstractVector; delta_t::Real=1.0, horizon::Integer=32)
    horizon > 0 || throw(ArgumentError("horizon must be positive, got $horizon"))
    kernel = irf.kernel_func(Float64.(collect(params)), Float64(delta_t), Int(horizon))
    return _pad_route_kernel(kernel, Int(horizon))
end

function build_irf_kernels(irf::RouteIRF, params::AbstractMatrix; delta_t::Real=1.0, horizon::Integer=32)
    kernels = Vector{Vector{Float64}}(undef, size(params, 1))
    for node in axes(params, 1)
        kernels[node] = irf(vec(@view(params[node, :])); delta_t=delta_t, horizon=horizon)
    end
    return kernels
end

function build_irf_kernels(irf::RouteIRF, params::AbstractVector{<:AbstractVector}; delta_t::Real=1.0, horizon::Integer=32)
    nn = isempty(params) ? 0 : length(first(params))
    kernels = Vector{Vector{Float64}}(undef, nn)
    for node in 1:nn
        kernels[node] = irf([param[node] for param in params]; delta_t=delta_t, horizon=horizon)
    end
    return kernels
end

struct SparseRouteKernel{IT,WT}
    indices::IT
    weights::WT
    n_out::Int
    n_in::Int
    horizon::Int

    function SparseRouteKernel(indices::AbstractMatrix{<:Integer}, weights::AbstractMatrix{<:Real}, n_out::Integer, n_in::Integer)
        size(indices, 2) == 2 || throw(ArgumentError("indices must have shape [nnz, 2]."))
        size(indices, 1) == size(weights, 1) || throw(ArgumentError("indices and weights must agree on nnz."))
        horizon = size(weights, 2)
        horizon > 0 || throw(ArgumentError("weights must contain at least one lag."))
        return new{Matrix{Int},Matrix{Float64}}(
            Int.(indices),
            Float64.(weights),
            Int(n_out),
            Int(n_in),
            horizon,
        )
    end
end

struct SparseRouteConvolution{KT,NT} <: AbstractHydroRoute
    name::Symbol
    kernel::KT
    infos::NT

    function SparseRouteConvolution(
        kernel::SparseRouteKernel;
        name::Optional{Symbol}=nothing,
        inputs::AbstractVector=Symbol[],
        outputs::AbstractVector=Symbol[],
        params::AbstractVector=Symbol[],
    )
        infos = HydroInfos(
            inputs=!isempty(inputs) ? tosymbol.(inputs) : Symbol[],
            outputs=!isempty(outputs) ? tosymbol.(outputs) : Symbol[],
            params=!isempty(params) ? tosymbol.(params) : Symbol[],
        )
        route_name = isnothing(name) ? Symbol("##sparse_route#", hash((kernel.n_out, kernel.n_in, kernel.horizon))) : name
        return new{typeof(kernel),typeof(infos)}(route_name, kernel, infos)
    end
end

@inline function _pad_route_kernel(kernel::AbstractVector, horizon::Int)
    padded = zeros(Float64, horizon)
    upper = min(length(kernel), horizon)
    @inbounds for idx in 1:upper
        padded[idx] = Float64(kernel[idx])
    end
    return padded
end

function _truncated_route_convolution(a::AbstractVector, b::AbstractVector, horizon::Int)
    out = zeros(Float64, horizon)
    upper_a = min(length(a), horizon)
    upper_b = min(length(b), horizon)
    @inbounds for ia in 1:upper_a
        aval = Float64(a[ia])
        iszero(aval) && continue
        max_ib = min(upper_b, horizon - ia + 1)
        for ib in 1:max_ib
            out[ia + ib - 1] += aval * Float64(b[ib])
        end
    end
    return out
end

function aggregate_route_kernel(
    adjacency::Union{Nothing,AbstractMatrix},
    node_kernels::AbstractVector{<:AbstractVector};
    topo_order::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    horizon::Integer,
)
    nn = length(node_kernels)
    horizon_int = Int(horizon)
    horizon_int > 0 || throw(ArgumentError("horizon must be positive, got $horizon"))

    if isnothing(adjacency)
        indices = hcat(collect(1:nn), collect(1:nn))
        weights = reduce(vcat, [_pad_route_kernel(node_kernels[node], horizon_int)' for node in 1:nn]; init=zeros(Float64, 0, horizon_int))
        return SparseRouteKernel(indices, weights, nn, nn)
    end

    size(adjacency, 1) == size(adjacency, 2) == nn || throw(ArgumentError("adjacency and node_kernels size mismatch."))
    order = isnothing(topo_order) ? collect(1:nn) : Int.(collect(topo_order))
    pairwise = [Dict{Int,Vector{Float64}}() for _ in 1:nn]

    for node in order
        self_kernel = _pad_route_kernel(node_kernels[node], horizon_int)
        pairwise[node][node] = copy(self_kernel)

        for upstream in 1:nn
            edge_weight = Float64(adjacency[node, upstream])
            iszero(edge_weight) && continue
            for (source, upstream_kernel) in pairwise[upstream]
                candidate = _truncated_route_convolution(self_kernel, upstream_kernel, horizon_int)
                edge_weight == 1.0 || (candidate .*= edge_weight)
                if haskey(pairwise[node], source)
                    pairwise[node][source] .+= candidate
                else
                    pairwise[node][source] = candidate
                end
            end
        end
    end

    nnz = sum(length(dict) for dict in pairwise)
    indices = Matrix{Int}(undef, nnz, 2)
    weights = Matrix{Float64}(undef, nnz, horizon_int)
    row = 1
    for dest in 1:nn
        for source in sort!(collect(keys(pairwise[dest])))
            indices[row, 1] = dest
            indices[row, 2] = source
            weights[row, :] .= pairwise[dest][source]
            row += 1
        end
    end

    return SparseRouteKernel(indices, weights, nn, nn)
end

function _causal_conv_add!(dest::AbstractVector, src::AbstractVector, kernel::AbstractVector)
    horizon = min(length(kernel), length(dest))
    @inbounds for tidx in eachindex(dest)
        upper = min(tidx, horizon)
        acc = zero(eltype(dest))
        for lag in 1:upper
            acc += kernel[lag] * src[tidx - lag + 1]
        end
        dest[tidx] += acc
    end
    return dest
end

function (conv::SparseRouteConvolution)(input::AbstractMatrix)
    kernel = conv.kernel
    size(input, 1) == kernel.n_in || throw(ArgumentError("Input node count $(size(input, 1)) does not match kernel.n_in $(kernel.n_in)."))
    T = promote_type(eltype(input), eltype(kernel.weights))
    out = zeros(T, kernel.n_out, size(input, 2))

    for row in axes(kernel.indices, 1)
        dest = kernel.indices[row, 1]
        source = kernel.indices[row, 2]
        _causal_conv_add!(view(out, dest, :), view(input, source, :), view(kernel.weights, row, :))
    end

    return out
end

function (conv::SparseRouteConvolution)(input::AbstractVector)
    out = conv(reshape(input, 1, :))
    return vec(out)
end

export HydroRoute,
       @hydroroute,
       build_aggr_func,
       RouteIRF,
       SparseRouteKernel,
       SparseRouteConvolution,
       build_irf_kernels,
       aggregate_route_kernel






