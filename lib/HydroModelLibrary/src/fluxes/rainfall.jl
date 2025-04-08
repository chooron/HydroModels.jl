@variables rainfall

RainfallFlux(::Val{:1}; input::NamedTuple, params::NamedTuple, output=rainfall) = begin
    @assert haskey(input, :P) "RainfallFlux{:1}: input must contain :P"
    @assert haskey(input, :T) "RainfallFlux{:1}: input must contain :T"
    @assert haskey(params, :thres) "RainfallFlux{:1}: params must contain :thres"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[ifelse(input.T > params.thres, input.P, 0)]
    )
end

RainfallFlux(::Val{:2}; input::NamedTuple, params::NamedTuple, output=rainfall) = begin
    @assert haskey(input, :P) "RainfallFlux{:2}: input must contain :P"
    @assert haskey(input, :T) "RainfallFlux{:2}: input must contain :T"
    @assert haskey(params, :p1) "RainfallFlux{:2}: params must contain :p1"
    @assert haskey(params, :p2) "RainfallFlux{:2}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[ifelse(input.T ≤ params.p1 - 1 / 2 * params.p2, 0,
            ifelse(input.T ≤ params.p1 + 1 / 2 * params.p2, input.P * (params.p1 + 1 / 2 * params.p2 - input.T) / params.p2, input.P))]
    )
end

export RainfallFlux