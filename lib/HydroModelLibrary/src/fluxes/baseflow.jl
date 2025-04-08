@variables baseflow

BaseflowFlux(::Val{1}; input::NamedTuple, params::NamedTuple, output=baseflow) = begin
    @assert haskey(input, :S) "BaseflowFlux{1}: input must contain :S"
    @assert haskey(params, :p1) "BaseflowFlux{1}: params must contain :p1"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps, exprs=[params.p1 * input.S])
end

BaseflowFlux(::Val{2}; input::NamedTuple, params::NamedTuple, output=baseflow) = begin
    @assert haskey(input, :S) "BaseflowFlux{2}: input must contain :S"
    @assert haskey(params, :p1) "BaseflowFlux{2}: params must contain :p1"
    @assert haskey(params, :p2) "BaseflowFlux{2}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps, exprs=[abs(1 / params.p1 * input.S)^(1 / params.p2)])
end

BaseflowFlux(::Val{3}; input::NamedTuple, params::NamedTuple, output=baseflow) = begin
    @assert haskey(input, :S) "BaseflowFlux{3}: input must contain :S"
    @assert haskey(params, :Smax) "BaseflowFlux{3}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps, exprs=[abs(params.Smax)^(-4) / (4) * (abs(input.S)^5)])
end

BaseflowFlux(::Val{4}; input::NamedTuple, params::NamedTuple, output=baseflow) = begin
    @assert haskey(input, :S) "BaseflowFlux{4}: input must contain :S"
    @assert haskey(params, :p1) "BaseflowFlux{4}: params must contain :p1"
    @assert haskey(params, :p2) "BaseflowFlux{4}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps, exprs=[params.p1 * exp(-params.p2 * input.S)])
end

BaseflowFlux(::Val{5}; input::NamedTuple, params::NamedTuple, output=baseflow) = begin
    @assert haskey(input, :S) "BaseflowFlux{5}: input must contain :S"
    @assert haskey(params, :p1) "BaseflowFlux{5}: params must contain :p1"
    @assert haskey(params, :p2) "BaseflowFlux{5}: params must contain :p2"
    @assert haskey(params, :Smax) "BaseflowFlux{5}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps, exprs=[params.p1 * (abs(input.S / params.Smax)^params.p2)])
end

BaseflowFlux(::Val{6}; input::NamedTuple, params::NamedTuple, output=baseflow) = begin
    @assert haskey(input, :S) "BaseflowFlux{6}: input must contain :S"
    @assert haskey(params, :p1) "BaseflowFlux{6}: params must contain :p1"
    @assert haskey(params, :p2) "BaseflowFlux{6}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps, exprs=[ifelse(input.S > params.p2, params.p1 * abs(input.S)^params.p2, 0)])
end

BaseflowFlux(::Val{7}; input::NamedTuple, params::NamedTuple, output=baseflow) = begin
    @assert haskey(input, :S) "BaseflowFlux{7}: input must contain :S"
    @assert haskey(params, :p1) "BaseflowFlux{7}: params must contain :p1"
    @assert haskey(params, :p2) "BaseflowFlux{7}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps, exprs=[params.p1 * abs(input.S)^params.p2])
end

BaseflowFlux(::Val{8}; input::NamedTuple, params::NamedTuple, output=baseflow) = begin
    @assert haskey(input, :S) "BaseflowFlux{8}: input must contain :S"
    @assert haskey(params, :p1) "BaseflowFlux{8}: params must contain :p1"
    @assert haskey(params, :p2) "BaseflowFlux{8}: params must contain :p2"
    @assert haskey(params, :Smax) "BaseflowFlux{8}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps, exprs=[params.p1 * (exp(-params.p2 * input.S / params.Smax) - 1)])
end

BaseflowFlux(::Val{9}; input::NamedTuple, params::NamedTuple, output=baseflow) = begin
    @assert haskey(input, :S) "BaseflowFlux{9}: input must contain :S"
    @assert haskey(params, :p1) "BaseflowFlux{9}: params must contain :p1"
    @assert haskey(params, :p2) "BaseflowFlux{9}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps, exprs=[ifelse(input.S > params.p2, params.p1 * (input.S - params.p2), 0)])
end

export BaseflowFlux