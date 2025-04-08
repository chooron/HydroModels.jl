@variables melt

MeltFlux(::Val{:1}; input::NamedTuple, params::NamedTuple, output=melt) = begin
    @assert haskey(input, :S) "MeltFlux{:1}: input must contain :S"
    @assert haskey(input, :t) "MeltFlux{:1}: input must contain :t"
    @assert haskey(params, :p1) "MeltFlux{:1}: params must contain :p1"
    @assert haskey(params, :t0) "MeltFlux{:1}: params must contain :t0"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * max(0, input.t - params.t0) * input.S]
    )
end

MeltFlux(::Val{:2}; input::NamedTuple, params::NamedTuple, output=melt) = begin
    @assert haskey(input, :S) "MeltFlux{:2}: input must contain :S"
    @assert haskey(input, :t) "MeltFlux{:2}: input must contain :t"
    @assert haskey(params, :p1) "MeltFlux{:2}: params must contain :p1"
    @assert haskey(params, :t0) "MeltFlux{:2}: params must contain :t0"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[min(input.S, params.p1 * max(0, input.t - params.t0))]
    )
end

MeltFlux(::Val{:3}; input::NamedTuple, params::NamedTuple, output=melt) = begin
    @assert haskey(input, :S) "MeltFlux{:3}: input must contain :S"
    @assert haskey(input, :t) "MeltFlux{:3}: input must contain :t"
    @assert haskey(params, :p1) "MeltFlux{:3}: params must contain :p1"
    @assert haskey(params, :t0) "MeltFlux{:3}: params must contain :t0"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * max(0, input.t - params.t0)]
    )
end

MeltFlux(::Val{:4}; input::NamedTuple, params::NamedTuple, output=melt) = begin
    @assert haskey(input, :S) "MeltFlux{:4}: input must contain :S"
    @assert haskey(input, :t) "MeltFlux{:4}: input must contain :t"
    @assert haskey(params, :p1) "MeltFlux{:4}: params must contain :p1"
    @assert haskey(params, :t0) "MeltFlux{:4}: params must contain :t0"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[min(input.S, params.p1 * input.t)]
    )
end

MeltFlux(::Val{:5}; input::NamedTuple, params::NamedTuple, output=melt) = begin
    @assert haskey(input, :S) "MeltFlux{:5}: input must contain :S"
    @assert haskey(input, :t) "MeltFlux{:5}: input must contain :t"
    @assert haskey(params, :p1) "MeltFlux{:5}: params must contain :p1"
    @assert haskey(params, :p2) "MeltFlux{:5}: params must contain :p2"
    @assert haskey(params, :t0) "MeltFlux{:5}: params must contain :t0"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[min(input.S, params.p1 * max(0, input.t - params.t0)^params.p2)]
    )
end

MeltFlux(::Val{:6}; input::NamedTuple, params::NamedTuple, output=melt) = begin
    @assert haskey(input, :S) "MeltFlux{:6}: input must contain :S"
    @assert haskey(input, :t) "MeltFlux{:6}: input must contain :t"
    @assert haskey(params, :p1) "MeltFlux{:6}: params must contain :p1"
    @assert haskey(params, :t0) "MeltFlux{:6}: params must contain :t0"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[min(input.S, params.p1 * (1 - exp(-max(0, input.t - params.t0))))]
    )
end

export MeltFlux