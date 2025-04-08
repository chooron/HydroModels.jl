function MuskingumRouteFlux(
    input::Num,
    output::Union{Num,Nothing}=nothing,
)
    @parameters k, x
    @variables s_river

    if isnothing(output)
        input_name = Symbolics.tosymbol(input, escape=false)
        output_name = Symbol(input_name, :_routed)
        output = first(@variables $output_name)
    end

    return RouteFlux(
        input,
        [k, x],
        [s_river],
        routetype=:muskingum,
        output=output
    )
end

function (flux::RouteFlux{:muskingum})(input::Matrix, pas::ComponentVector; kwargs...)
    timeidx = get(kwargs, :timeidx, collect(1:size(input)[2]))
    input_itp = LinearInterpolation(input[1, :], timeidx)
    params = pas[:params]

    function msk_prob!(du, u, p, t)
        s_river = u[1]
        q_in = input_itp(t)
        k, x = p
        q_out = (s_river - k * x * q_in) / (k * (1 - x))
        du[1] = q_in - q_out
    end

    init_states = [params.k * input[1, 1]]
    prob = ODEProblem(msk_prob!, init_states, (timeidx[1], timeidx[end]), params)
    sol = solve(prob, Rosenbrock23(), saveat=timeidx)

    s_river_vec = Array(sol)
    q_out_vec = @.((s_river_vec - params.k * params.x * input) / (params.k * (1 - params.x)))
    q_out_vec
end

function get_rflux_initstates(::RouteFlux{:muskingum}; input::AbstractMatrix, pas::ComponentVector, stypes::AbstractVector{Symbol}, ptypes::AbstractVector{Symbol})
    [pas[:params][ptype][:k] for ptype in ptypes] .* input[:, 1]
end

function get_rflux_func(::RouteFlux{:muskingum}; pas::ComponentVector, ptypes::AbstractVector{Symbol})

    function rflux_func(s_river, q_in, q_gen, p)
        k_ps = [p[ptype][:k] for ptype in ptypes]
        x_ps = [p[ptype][:x] for ptype in ptypes]

        q_rf = @.((s_river - k_ps * x_ps * q_in) / (k_ps * (1 - x_ps)))
        d_state = q_in .- q_rf
        q_out = q_rf .+ q_gen
        q_out, d_state
    end

    return rflux_func
end
