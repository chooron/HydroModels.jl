@variables refreeze

RefreezeFlux(::Val{:1}; input::NamedTuple, params::NamedTuple, output=refreeze) = begin
    @assert haskey(input, :S) "RefreezeFlux{:1}: input must contain :S"
    @assert haskey(input, :t) "RefreezeFlux{:1}: input must contain :t"
    @assert haskey(params, :p1) "RefreezeFlux{:1}: params must contain :p1"
    @assert haskey(params, :t0) "RefreezeFlux{:1}: params must contain :t0"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * min(0, input.t - params.t0) * input.S]
    )
end

RefreezeFlux(::Val{:2}; input::NamedTuple, params::NamedTuple, output=refreeze) = begin
    @assert haskey(input, :S) "RefreezeFlux{:2}: input must contain :S"
    @assert haskey(input, :t) "RefreezeFlux{:2}: input must contain :t"
    @assert haskey(params, :p1) "RefreezeFlux{:2}: params must contain :p1"
    @assert haskey(params, :t0) "RefreezeFlux{:2}: params must contain :t0"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[min(input.S, params.p1 * min(0, input.t - params.t0))]
    )
end

RefreezeFlux(::Val{:3}; input::NamedTuple, params::NamedTuple, output=refreeze) = begin
    @assert haskey(input, :S) "RefreezeFlux{:3}: input must contain :S"
    @assert haskey(input, :t) "RefreezeFlux{:3}: input must contain :t"
    @assert haskey(params, :p1) "RefreezeFlux{:3}: params must contain :p1"
    @assert haskey(params, :t0) "RefreezeFlux{:3}: params must contain :t0"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * min(0, input.t - params.t0)]
    )
end

RefreezeFlux(::Val{:4}; input::NamedTuple, params::NamedTuple, output=refreeze) = begin
    @assert haskey(input, :S) "RefreezeFlux{:4}: input must contain :S"
    @assert haskey(input, :t) "RefreezeFlux{:4}: input must contain :t"
    @assert haskey(params, :p1) "RefreezeFlux{:4}: params must contain :p1"
    @assert haskey(params, :p2) "RefreezeFlux{:4}: params must contain :p2"
    @assert haskey(params, :t0) "RefreezeFlux{:4}: params must contain :t0"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[min(input.S, params.p1 * min(0, input.t - params.t0)^params.p2)]
    )
end

RefreezeFlux(::Val{:5}; input::NamedTuple, params::NamedTuple, output=refreeze) = begin
    @assert haskey(input, :S) "RefreezeFlux{:5}: input must contain :S"
    @assert haskey(input, :t) "RefreezeFlux{:5}: input must contain :t"
    @assert haskey(params, :p1) "RefreezeFlux{:5}: params must contain :p1"
    @assert haskey(params, :t0) "RefreezeFlux{:5}: params must contain :t0"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[min(input.S, params.p1 * (1 - exp(min(0, input.t - params.t0))))]
    )
end

export RefreezeFlux