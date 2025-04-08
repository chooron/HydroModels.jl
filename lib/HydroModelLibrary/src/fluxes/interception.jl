@variables interception

InterceptionFlux(::Val{:1}; input::NamedTuple, params::NamedTuple, output=interception) = begin
    @assert haskey(input, :S) "InterceptionFlux{:1}: input must contain :S"
    @assert haskey(input, :pet) "InterceptionFlux{:1}: input must contain :pet"
    @assert haskey(params, :Smax) "InterceptionFlux{:1}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[min(input.S, input.pet)]
    )
end

InterceptionFlux(::Val{:2}; input::NamedTuple, params::NamedTuple, output=interception) = begin
    @assert haskey(input, :S) "InterceptionFlux{:2}: input must contain :S"
    @assert haskey(input, :pet) "InterceptionFlux{:2}: input must contain :pet"
    @assert haskey(params, :Smax) "InterceptionFlux{:2}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[min(input.S, input.pet * input.S / params.Smax)]
    )
end

InterceptionFlux(::Val{:3}; input::NamedTuple, params::NamedTuple, output=interception) = begin
    @assert haskey(input, :S) "InterceptionFlux{:3}: input must contain :S"
    @assert haskey(input, :pet) "InterceptionFlux{:3}: input must contain :pet"
    @assert haskey(params, :Smax) "InterceptionFlux{:3}: params must contain :Smax"
    @assert haskey(params, :p1) "InterceptionFlux{:3}: params must contain :p1"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[min(input.S, input.pet * (input.S / params.Smax)^params.p1)]
    )
end

InterceptionFlux(::Val{:4}; input::NamedTuple, params::NamedTuple, output=interception) = begin
    @assert haskey(input, :S) "InterceptionFlux{:4}: input must contain :S"
    @assert haskey(input, :pet) "InterceptionFlux{:4}: input must contain :pet"
    @assert haskey(params, :Smax) "InterceptionFlux{:4}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[min(input.S, input.pet * (1 - exp(-input.S / params.Smax)))]
    )
end

InterceptionFlux(::Val{:5}; input::NamedTuple, params::NamedTuple, output=interception) = begin
    @assert haskey(input, :S) "InterceptionFlux{:5}: input must contain :S"
    @assert haskey(input, :pet) "InterceptionFlux{:5}: input must contain :pet"
    @assert haskey(params, :Smax) "InterceptionFlux{:5}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[min(input.S, input.pet * (input.S / params.Smax))]
    )
end

InterceptionFlux(::Val{:6}; input::NamedTuple, params::NamedTuple, output=interception) = begin
    @assert haskey(input, :S) "InterceptionFlux{:6}: input must contain :S"
    @assert haskey(input, :pet) "InterceptionFlux{:6}: input must contain :pet"
    @assert haskey(params, :Smax) "InterceptionFlux{:6}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[min(input.S, input.pet * (input.S / params.Smax)^2)]
    )
end

InterceptionFlux(::Val{:7}; input::NamedTuple, params::NamedTuple, output=interception) = begin
    @assert haskey(input, :S) "InterceptionFlux{:7}: input must contain :S"
    @assert haskey(input, :pet) "InterceptionFlux{:7}: input must contain :pet"
    @assert haskey(params, :Smax) "InterceptionFlux{:7}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[min(input.S, input.pet * (input.S / params.Smax)^0.5)]
    )
end

export InterceptionFlux