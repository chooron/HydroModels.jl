function CascadeRouteFlux(
    input::Num,
    output::Union{Num,Nothing}=nothing,
)
    @parameters k [description = "水库的平均滞留时间"]
    @parameters n [description = "水库的数目"]

    if isnothing(output)
        input_name = Symbolics.tosymbol(input, escape=false)
        output_name = Symbol(input_name, :_routed)
        output = first(@variables $output_name)
    end

    return RouteFlux(
        input,
        [k, n],
        Num[],
        routetype=:cascade,
        output=output
    )
end

function (flux::RouteFlux{:cascade})(input::AbstractMatrix, pas::ComponentVector; kwargs::NamedTuple=NamedTuple())
    n = Int(pas[:params].n)
    init_states = zeros(n)
    input_len = size(input)[2]
    input_itp = LinearInterpolation(input[1, :], collect(1:input_len))

    function nash_unithydro!(du, u, p, t)
        k = p
        du[1] = input_itp(t) - u[1] / k
        for i in 2:n
            du[i] = (u[i-1] - u[i]) / k
        end
    end

    prob = ODEProblem(nash_unithydro!, init_states, (1, length(input_vec)), (pas[:params].k,))
    sol = solve(prob, Tsit5())
    sol_vec = sol.u[:, end] .* pas[:params].k
    reshape(sol_vec, 1, input_len)
end

function get_rflux_initstates(::RouteFlux{:cascade}; input::AbstractMatrix, pas::ComponentVector, stypes::AbstractVector{Symbol}, ptypes::AbstractVector{Symbol})
    reduce(vcat, [zeros(eltype(pas), Int(pas[:params][ptype][:n])) for ptype in ptypes])
end

function get_rflux_func(::RouteFlux{:cascade}; pas::ComponentVector, ptypes::AbstractVector{Symbol})
    #* node_num * ts_len
    node_iuh_nums = [pas[:params][ptype][:n] for ptype in ptypes]
    start_idxes = Int[1; cumsum(node_iuh_nums)[1:end-1] .+ 1]
    end_idxes = Int.(cumsum(node_iuh_nums))
    iuh_states_idxes = [Int.(sidx:eidx) for (sidx, eidx) in zip(start_idxes, end_idxes)]

    function rflux_func(uh_states, q_in, q_gen, p)
        k_ps = [p[ptype][:k] for ptype in ptypes]
        dstates = zeros(eltype(uh_states), length(uh_states))
        dstates[start_idxes] = @.(q_in + q_gen - uh_states[start_idxes] / k_ps)
        for (i, (k, n)) in enumerate(zip(k_ps, node_iuh_nums))
            iuh_states_idx = iuh_states_idxes[i]
            for j in 2:Int(n)
                dstates[iuh_states_idx[j]] = (uh_states[iuh_states_idx[j-1]] - uh_states[iuh_states_idx[j]]) / k
            end
        end
        q_out = [k_ps[i] * uh_states[end_idxes[i]] for i in eachindex(k_ps)]
        q_out, dstates
    end

    return rflux_func
end
