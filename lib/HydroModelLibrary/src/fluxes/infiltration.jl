@variables infiltration

InfiltrationFlux(::Val{:1}; input::NamedTuple, params::NamedTuple, output=infiltration) = begin
    @assert haskey(input, :S) "InfiltrationFlux{:1}: input must contain :S"
    @assert haskey(params, :p1) "InfiltrationFlux{:1}: params must contain :p1"
    @assert haskey(params, :p2) "InfiltrationFlux{:1}: params must contain :p2"
    @assert haskey(params, :Smax) "InfiltrationFlux{:1}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * exp(-params.p2 * input.S / params.Smax)]
    )
end

InfiltrationFlux(::Val{:2}; input::NamedTuple, params::NamedTuple, output=infiltration) = begin
    @assert haskey(input, :S) "InfiltrationFlux{:2}: input must contain :S"
    @assert haskey(input, :used) "InfiltrationFlux{:2}: input must contain :used"
    @assert haskey(params, :p1) "InfiltrationFlux{:2}: params must contain :p1"
    @assert haskey(params, :p2) "InfiltrationFlux{:2}: params must contain :p2"
    @assert haskey(params, :Smax) "InfiltrationFlux{:2}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * exp(-params.p2 * input.S / params.Smax) - input.used]
    )
end

InfiltrationFlux(::Val{:3}; input::NamedTuple, params::NamedTuple, output=infiltration) = begin
    @assert haskey(input, :S) "InfiltrationFlux{:3}: input must contain :S"
    @assert haskey(input, :in) "InfiltrationFlux{:3}: input must contain :in"
    @assert haskey(params, :Smax) "InfiltrationFlux{:3}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps, exprs=[ifelse(input.S > params.Smax, input.in, 0)])
end

InfiltrationFlux(::Val{:4}; input::NamedTuple, params::NamedTuple, output=infiltration) = begin
    @assert haskey(params, :p1) "InfiltrationFlux{:4}: params must contain :p1"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps, exprs=[params.p1])
end

InfiltrationFlux(::Val{:5}; input::NamedTuple, params::NamedTuple, output=infiltration) = begin
    @assert haskey(input, :S1) "InfiltrationFlux{:5}: input must contain :S1"
    @assert haskey(input, :S2) "InfiltrationFlux{:5}: input must contain :S2"
    @assert haskey(params, :S1max) "InfiltrationFlux{:5}: params must contain :S1max"
    @assert haskey(params, :S2max) "InfiltrationFlux{:5}: params must contain :S2max"
    @assert haskey(params, :p1) "InfiltrationFlux{:5}: params must contain :p1"
    @assert haskey(params, :p2) "InfiltrationFlux{:5}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * (1 - input.S1 / params.S1max) * abs(input.S2 / params.S2max)^(-params.p2)]
    )
end

InfiltrationFlux(::Val{:6}; input::NamedTuple, params::NamedTuple, output=infiltration) = begin
    @assert haskey(input, :S) "InfiltrationFlux{:6}: input must contain :S"
    @assert haskey(input, :in) "InfiltrationFlux{:6}: input must contain :in"
    @assert haskey(params, :Smax) "InfiltrationFlux{:6}: params must contain :Smax"
    @assert haskey(params, :p1) "InfiltrationFlux{:6}: params must contain :p1"
    @assert haskey(params, :p2) "InfiltrationFlux{:6}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * abs(input.S / params.Smax)^params.p2 * input.in]
    )
end

InfiltrationFlux(::Val{:7}; input::NamedTuple, params::NamedTuple, output=infiltration) = begin
    @assert haskey(input, :S) "InfiltrationFlux{:7}: input must contain :S"
    @assert haskey(params, :Smax) "InfiltrationFlux{:7}: params must contain :Smax"
    @assert haskey(params, :p1) "InfiltrationFlux{:7}: params must contain :p1"
    @assert haskey(params, :p2) "InfiltrationFlux{:7}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[ifelse(input.S ≥ params.Smax, params.p1 * exp(-params.p2 * input.S / params.Smax), 0)]
    )
end

export InfiltrationFlux