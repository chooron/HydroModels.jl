@variables recharge

RechargeFlux(::Val{:1}; input::NamedTuple, params::NamedTuple, output=recharge) = begin
    @assert haskey(input, :S) "RechargeFlux{:1}: input must contain :S"
    @assert haskey(params, :p1) "RechargeFlux{:1}: params must contain :p1"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * input.S]
    )
end

RechargeFlux(::Val{:2}; input::NamedTuple, params::NamedTuple, output=recharge) = begin
    @assert haskey(input, :S) "RechargeFlux{:2}: input must contain :S"
    @assert haskey(params, :p1) "RechargeFlux{:2}: params must contain :p1"
    @assert haskey(params, :p2) "RechargeFlux{:2}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * input.S^params.p2]
    )
end

RechargeFlux(::Val{:3}; input::NamedTuple, params::NamedTuple, output=recharge) = begin
    @assert haskey(input, :S) "RechargeFlux{:3}: input must contain :S"
    @assert haskey(params, :p1) "RechargeFlux{:3}: params must contain :p1"
    @assert haskey(params, :Smax) "RechargeFlux{:3}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * input.S / params.Smax]
    )
end

RechargeFlux(::Val{:4}; input::NamedTuple, params::NamedTuple, output=recharge) = begin
    @assert haskey(input, :S) "RechargeFlux{:4}: input must contain :S"
    @assert haskey(params, :p1) "RechargeFlux{:4}: params must contain :p1"
    @assert haskey(params, :p2) "RechargeFlux{:4}: params must contain :p2"
    @assert haskey(params, :Smax) "RechargeFlux{:4}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * (input.S / params.Smax)^params.p2]
    )
end

RechargeFlux(::Val{:5}; input::NamedTuple, params::NamedTuple, output=recharge) = begin
    @assert haskey(input, :S) "RechargeFlux{:5}: input must contain :S"
    @assert haskey(params, :p1) "RechargeFlux{:5}: params must contain :p1"
    @assert haskey(params, :Smax) "RechargeFlux{:5}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * max(0, input.S - params.Smax)]
    )
end

RechargeFlux(::Val{:6}; input::NamedTuple, params::NamedTuple, output=recharge) = begin
    @assert haskey(input, :S) "RechargeFlux{:6}: input must contain :S"
    @assert haskey(params, :p1) "RechargeFlux{:6}: params must contain :p1"
    @assert haskey(params, :p2) "RechargeFlux{:6}: params must contain :p2"
    @assert haskey(params, :Smax) "RechargeFlux{:6}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * max(0, input.S - params.Smax)^params.p2]
    )
end

export RechargeFlux