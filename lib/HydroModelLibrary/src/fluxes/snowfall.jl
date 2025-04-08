@variables snowfall

SnowfallFlux(::Val{:1}; input::NamedTuple, params::NamedTuple, output=snowfall) = begin
    @assert haskey(input, :P) "SnowfallFlux{:1}: input must contain :P"
    @assert haskey(input, :T) "SnowfallFlux{:1}: input must contain :T"
    @assert haskey(params, :thres) "SnowfallFlux{:1}: params must contain :thres"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[ifelse(input.T > params.thres, 0, input.P)]
    )
end

SnowfallFlux(::Val{:2}; input::NamedTuple, params::NamedTuple, output=snowfall) = begin
    @assert haskey(input, :P) "SnowfallFlux{:2}: input must contain :P"
    @assert haskey(input, :T) "SnowfallFlux{:2}: input must contain :T"
    @assert haskey(params, :p1) "SnowfallFlux{:2}: params must contain :p1"
    @assert haskey(params, :p2) "SnowfallFlux{:2}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[ifelse(input.T ≤ params.p1 - 1 / 2 * params.p2, input.P,
            ifelse(input.T ≤ params.p1 + 1 / 2 * params.p2, input.P * (params.p1 + 1 / 2 * params.p2 - input.T) / params.p2, 0))]
    )
end

export SnowfallFlux