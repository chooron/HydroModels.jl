@variables interflow

InterflowFlux(::Val{:1}; input::NamedTuple, params::NamedTuple, output=interflow) = begin
    @assert haskey(input, :S) "InterflowFlux{:1}: input must contain :S"
    @assert haskey(params, :p1) "InterflowFlux{:1}: params must contain :p1"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * input.S]
    )
end

InterflowFlux(::Val{:2}; input::NamedTuple, params::NamedTuple, output=interflow) = begin
    @assert haskey(input, :S) "InterflowFlux{:2}: input must contain :S"
    @assert haskey(params, :p1) "InterflowFlux{:2}: params must contain :p1"
    @assert haskey(params, :p2) "InterflowFlux{:2}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * input.S^params.p2]
    )
end

InterflowFlux(::Val{:3}; input::NamedTuple, params::NamedTuple, output=interflow) = begin
    @assert haskey(input, :S) "InterflowFlux{:3}: input must contain :S"
    @assert haskey(params, :p1) "InterflowFlux{:3}: params must contain :p1"
    @assert haskey(params, :p2) "InterflowFlux{:3}: params must contain :p2"
    @assert haskey(params, :Smax) "InterflowFlux{:3}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * (input.S / params.Smax)^params.p2]
    )
end

InterflowFlux(::Val{:4}; input::NamedTuple, params::NamedTuple, output=interflow) = begin
    @assert haskey(input, :S) "InterflowFlux{:4}: input must contain :S"
    @assert haskey(params, :p1) "InterflowFlux{:4}: params must contain :p1"
    @assert haskey(params, :p2) "InterflowFlux{:4}: params must contain :p2"
    @assert haskey(params, :Smax) "InterflowFlux{:4}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * input.S * (input.S / params.Smax)^params.p2]
    )
end

InterflowFlux(::Val{:5}; input::NamedTuple, params::NamedTuple, output=interflow) = begin
    @assert haskey(input, :S) "InterflowFlux{:5}: input must contain :S"
    @assert haskey(params, :p1) "InterflowFlux{:5}: params must contain :p1"
    @assert haskey(params, :Smax) "InterflowFlux{:5}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * max(0, input.S - params.Smax)]
    )
end

InterflowFlux(::Val{:6}; input::NamedTuple, params::NamedTuple, output=interflow) = begin
    @assert haskey(input, :S) "InterflowFlux{:6}: input must contain :S"
    @assert haskey(params, :p1) "InterflowFlux{:6}: params must contain :p1"
    @assert haskey(params, :p2) "InterflowFlux{:6}: params must contain :p2"
    @assert haskey(params, :Smax) "InterflowFlux{:6}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * max(0, input.S - params.Smax)^params.p2]
    )
end

InterflowFlux(::Val{:7}; input::NamedTuple, params::NamedTuple, output=interflow) = begin
    @assert haskey(input, :S) "InterflowFlux{:7}: input must contain :S"
    @assert haskey(params, :p1) "InterflowFlux{:7}: params must contain :p1"
    @assert haskey(params, :p2) "InterflowFlux{:7}: params must contain :p2"
    @assert haskey(params, :Smax) "InterflowFlux{:7}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * input.S * (max(0, input.S - params.Smax) / params.Smax)^params.p2]
    )
end

InterflowFlux(::Val{:8}; input::NamedTuple, params::NamedTuple, output=interflow) = begin
    @assert haskey(input, :S) "InterflowFlux{:8}: input must contain :S"
    @assert haskey(params, :p1) "InterflowFlux{:8}: params must contain :p1"
    @assert haskey(params, :Smax) "InterflowFlux{:8}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * input.S * (1 - exp(-max(0, input.S - params.Smax) / params.Smax))]
    )
end

export InterflowFlux