@variables exchange

ExchangeFlux(::Val{:1}; input::NamedTuple, params::NamedTuple, output=exchange) = begin
    @assert haskey(input, :S) "ExchangeFlux{:1}: input must contain :S"
    @assert haskey(params, :p1) "ExchangeFlux{:1}: params must contain :p1"
    @assert haskey(params, :p2) "ExchangeFlux{:1}: params must contain :p2"
    @assert haskey(params, :p3) "ExchangeFlux{:1}: params must contain :p3"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[(params.p1 * abs(input.S) + params.p2 * (1 - exp(-params.p3 * abs(input.S)))) * sign(input.S)]
    )
end

ExchangeFlux(::Val{:2}; input::NamedTuple, params::NamedTuple, output=exchange) = begin
    @assert haskey(input, :S1) "ExchangeFlux{:2}: input must contain :S1"
    @assert haskey(input, :S2) "ExchangeFlux{:2}: input must contain :S2"
    @assert haskey(params, :S1max) "ExchangeFlux{:2}: params must contain :S1max"
    @assert haskey(params, :S2max) "ExchangeFlux{:2}: params must contain :S2max"
    @assert haskey(params, :p1) "ExchangeFlux{:2}: params must contain :p1"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * (input.S1 / params.S1max - input.S2 / params.S2max)]
    )
end

ExchangeFlux(::Val{:3}; input::NamedTuple, params::NamedTuple, output=exchange) = begin
    @assert haskey(input, :S) "ExchangeFlux{:3}: input must contain :S"
    @assert haskey(params, :p1) "ExchangeFlux{:3}: params must contain :p1"
    @assert haskey(params, :p2) "ExchangeFlux{:3}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * (input.S - params.p2)]
    )
end

export ExchangeFlux