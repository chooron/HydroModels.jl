@variables saturation

SaturationFlux(::Val{:1}; input::NamedTuple, params::NamedTuple, output=saturation) = begin
    @assert haskey(input, :S) "SaturationFlux{:1}: input must contain :S"
    @assert haskey(params, :Smax) "SaturationFlux{:1}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[ifelse(input.S > params.Smax, input.in, 0)]
    )
end

SaturationFlux(::Val{:2}; input::NamedTuple, params::NamedTuple, output=saturation) = begin
    @assert haskey(input, :in) "SaturationFlux{:2}: input must contain :in"
    @assert haskey(params, :p1) "SaturationFlux{:2}: params must contain :p1"
    @assert haskey(params, :Smax) "SaturationFlux{:2}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[(1 - min(1, (1 - input.S / params.Smax)^params.p1)) * input.in]
    )
end

SaturationFlux(::Val{:3}; input::NamedTuple, params::NamedTuple, output=saturation) = begin
    @assert haskey(input, :S) "SaturationFlux{:3}: input must contain :S"
    @assert haskey(input, :in) "SaturationFlux{:3}: input must contain :in"
    @assert haskey(params, :p1) "SaturationFlux{:3}: params must contain :p1"
    @assert haskey(params, :Smax) "SaturationFlux{:3}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[(1 - min(1, 1 / (1 + exp((input.S / params.Smax + 0.5) / params.p1)))) * input.in]
    )
end

SaturationFlux(::Val{:4}; input::NamedTuple, params::NamedTuple, output=saturation) = begin
    @assert haskey(input, :S) "SaturationFlux{:4}: input must contain :S"
    @assert haskey(input, :in) "SaturationFlux{:4}: input must contain :in"
    @assert haskey(params, :p1) "SaturationFlux{:4}: params must contain :p1"
    @assert haskey(params, :Smax) "SaturationFlux{:4}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[(1 - min(1, input.S / params.Smax)^2) * input.in]
    )
end

SaturationFlux(::Val{:5}; input::NamedTuple, params::NamedTuple, output=saturation) = begin
    @assert haskey(input, :S) "SaturationFlux{:5}: input must contain :S"
    @assert haskey(input, :in) "SaturationFlux{:5}: input must contain :in"
    @assert haskey(params, :p1) "SaturationFlux{:5}: params must contain :p1"
    @assert haskey(params, :p2) "SaturationFlux{:5}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[(1 - min(1, (input.S / params.p1)^params.p2)) * input.in]
    )
end

SaturationFlux(::Val{:6}; input::NamedTuple, params::NamedTuple, output=saturation) = begin
    @assert haskey(input, :S) "SaturationFlux{:6}: input must contain :S"
    @assert haskey(input, :in) "SaturationFlux{:6}: input must contain :in"
    @assert haskey(params, :p1) "SaturationFlux{:6}: params must contain :p1"
    @assert haskey(params, :Smax) "SaturationFlux{:6}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * input.S / params.Smax * input.in]
    )
end

#TODO Saturation excess from a store with different degrees of saturation (gamma function variant) 

SaturationFlux(::Val{:8}; input::NamedTuple, params::NamedTuple, output=saturation) = begin
    @assert haskey(input, :S) "SaturationFlux{:8}: input must contain :S"
    @assert haskey(input, :in) "SaturationFlux{:8}: input must contain :in"
    @assert haskey(params, :p1) "SaturationFlux{:8}: params must contain :p1"
    @assert haskey(params, :p2) "SaturationFlux{:8}: params must contain :p2"
    @assert haskey(params, :Smax) "SaturationFlux{:8}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[(params.p1 + max(0, params.p2 - params.p1) * input.S / params.Smax) * input.in]
    )
end

SaturationFlux(::Val{:9}; input::NamedTuple, params::NamedTuple, output=saturation) = begin
    @assert haskey(input, :S) "SaturationFlux{:9}: input must contain :S"
    @assert haskey(input, :in) "SaturationFlux{:9}: input must contain :in"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[ifelse(input.S == 0, input.in, 0)]
    )
end

SaturationFlux(::Val{:10}; input::NamedTuple, params::NamedTuple, output=saturation) = begin
    @assert haskey(input, :S) "SaturationFlux{:10}: input must contain :S"
    @assert haskey(input, :in) "SaturationFlux{:10}: input must contain :in"
    @assert haskey(params, :p1) "SaturationFlux{:10}: params must contain :p1"
    @assert haskey(params, :p2) "SaturationFlux{:10}: params must contain :p2"
    @assert haskey(params, :p3) "SaturationFlux{:10}: params must contain :p3"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[input.in * min(params.p1, params.p2 + params.p2 * exp(params.p3 * input.S))]
    )
end

SaturationFlux(::Val{:11}; input::NamedTuple, params::NamedTuple, output=saturation) = begin
    @assert haskey(input, :S) "SaturationFlux{:11}: input must contain :S"
    @assert haskey(input, :in) "SaturationFlux{:11}: input must contain :in"
    @assert haskey(params, :p1) "SaturationFlux{:11}: params must contain :p1"
    @assert haskey(params, :p2) "SaturationFlux{:11}: params must contain :p2"
    @assert haskey(params, :Smax) "SaturationFlux{:11}: params must contain :Smax"
    @assert haskey(params, :Smin) "SaturationFlux{:11}: params must contain :Smin"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * clamp((input.S - params.Smin) / (params.Smax - params.Smin), 0, 1)^params.p2 * input.in]
    )
end

SaturationFlux(::Val{:12}; input::NamedTuple, params::NamedTuple, output=saturation) = begin
    @assert haskey(input, :in) "SaturationFlux{:12}: input must contain :in"
    @assert haskey(params, :p1) "SaturationFlux{:12}: params must contain :p1"
    @assert haskey(params, :p2) "SaturationFlux{:12}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[(params.p1 - params.p2) / (1 - params.p2) * input.in]
    )
end

# TODO Saturation excess flow from a store with different degrees of saturation (normal distribution variant) 

SaturationFlux(::Val{:14}; input::NamedTuple, params::NamedTuple, output=saturation) = begin
    @assert haskey(input, :S) "SaturationFlux{:14}: input must contain :S"
    @assert haskey(input, :in) "SaturationFlux{:14}: input must contain :in"
    @assert haskey(params, :p1) "SaturationFlux{:14}: params must contain :p1"
    @assert haskey(params, :p2) "SaturationFlux{:14}: params must contain :p2"
    @assert haskey(params, :p3) "SaturationFlux{:14}: params must contain :p3"
    @assert haskey(params, :Smax) "SaturationFlux{:14}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[ifelse(input.S / params.Smax ≤ 0.5 - params.p1,
            (0.5 - params.p1)^(1 - params.p2) * max(0, input.S / params.Smax)^params.p3,
            1 - (0.5 - params.p1)^(1 - params.p2) * max(0, 1 - input.S / params.Smax)^params.p3) * input.in]
    )
end

export SaturationFlux