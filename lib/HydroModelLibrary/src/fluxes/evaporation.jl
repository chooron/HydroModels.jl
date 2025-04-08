@variables evaporation

EvaporationFlux(::Val{1}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S) "EvaporationFlux{1}: input must contain :S"
    @assert haskey(input, :pet) "EvaporationFlux{1}: input must contain :pet"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps, exprs=[ifelse(input.S > 0, input.pet, 0)])
end

EvaporationFlux(::Val{2}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S) "EvaporationFlux{2}: input must contain :S"
    @assert haskey(params, :p1) "EvaporationFlux{2}: params must contain :p1"
    @assert haskey(params, :Smax) "EvaporationFlux{2}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps, exprs=[params.p1 * input.S / params.Smax])
end

EvaporationFlux(::Val{3}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S) "EvaporationFlux{3}: input must contain :S"
    @assert haskey(input, :pet) "EvaporationFlux{3}: input must contain :pet"
    @assert haskey(params, :p1) "EvaporationFlux{3}: params must contain :p1"
    @assert haskey(params, :Smax) "EvaporationFlux{3}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[ifelse(input.S < params.p1 * params.Smax, input.pet * input.S / params.p1 / params.Smax, input.pet)])
end

EvaporationFlux(::Val{4}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S) "EvaporationFlux{4}: input must contain :S"
    @assert haskey(input, :pet) "EvaporationFlux{4}: input must contain :pet"
    @assert haskey(params, :p1) "EvaporationFlux{4}: params must contain :p1"
    @assert haskey(params, :p2) "EvaporationFlux{4}: params must contain :p2"
    @assert haskey(params, :Smax) "EvaporationFlux{4}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[input.pet * max(0, params.p1 * (input.S - params.p2 * params.Smax) / (params.Smax - params.p2 * params.Smax))])
end

EvaporationFlux(::Val{5}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S) "EvaporationFlux{5}: input must contain :S"
    @assert haskey(input, :pet) "EvaporationFlux{5}: input must contain :pet"
    @assert haskey(params, :p1) "EvaporationFlux{5}: params must contain :p1"
    @assert haskey(params, :Smax) "EvaporationFlux{5}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[(1 - params.p1) * input.S / params.Smax * input.pet])
end

EvaporationFlux(::Val{6}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S) "EvaporationFlux{6}: input must contain :S"
    @assert haskey(input, :pet) "EvaporationFlux{6}: input must contain :pet"
    @assert haskey(params, :p1) "EvaporationFlux{6}: params must contain :p1"
    @assert haskey(params, :p2) "EvaporationFlux{6}: params must contain :p2"
    @assert haskey(params, :Smax) "EvaporationFlux{6}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[min(1.0, input.S / params.p2 * params.Smax) * params.p1 * input.pet])
end

EvaporationFlux(::Val{7}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S) "EvaporationFlux{7}: input must contain :S"
    @assert haskey(input, :pet) "EvaporationFlux{7}: input must contain :pet"
    @assert haskey(params, :Smax) "EvaporationFlux{7}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[input.S / params.Smax * input.pet])
end

EvaporationFlux(::Val{8}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S1) "EvaporationFlux{8}: input must contain :S1"
    @assert haskey(input, :S2) "EvaporationFlux{8}: input must contain :S2"
    @assert haskey(input, :pet) "EvaporationFlux{8}: input must contain :pet"
    @assert haskey(params, :Smax) "EvaporationFlux{8}: params must contain :Smax"
    @assert haskey(params, :p1) "EvaporationFlux{8}: params must contain :p1"
    @assert haskey(params, :p2) "EvaporationFlux{8}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[input.S1 / (input.S1 + input.S2) * params.p1 * input.pet * min(1, input.S1 / params.p2)])
end

EvaporationFlux(::Val{9}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S1) "EvaporationFlux{9}: input must contain :S1"
    @assert haskey(input, :S2) "EvaporationFlux{9}: input must contain :S2"
    @assert haskey(input, :pet) "EvaporationFlux{9}: input must contain :pet"
    @assert haskey(params, :Smax) "EvaporationFlux{9}: params must contain :Smax"
    @assert haskey(params, :p1) "EvaporationFlux{9}: params must contain :p1"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[input.S1 / (input.S1 + input.S2) * (1 - params.p1) * (input.S1 / (params.Smax - input.S2))])
end

EvaporationFlux(::Val{10}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S) "EvaporationFlux{10}: input must contain :S"
    @assert haskey(input, :pet) "EvaporationFlux{10}: input must contain :pet"
    @assert haskey(params, :Smax) "EvaporationFlux{10}: params must contain :Smax"
    @assert haskey(params, :p1) "EvaporationFlux{10}: params must contain :p1"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[input.S / params.Smax * params.p1 * input.pet])
end

EvaporationFlux(::Val{11}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S) "EvaporationFlux{11}: input must contain :S"
    @assert haskey(input, :pet) "EvaporationFlux{11}: input must contain :pet"
    @assert haskey(params, :Smax) "EvaporationFlux{11}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[(2 * input.S / params.Smax - (input.S / params.Smax)^2) * input.pet])
end

EvaporationFlux(::Val{12}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S) "EvaporationFlux{12}: input must contain :S"
    @assert haskey(input, :pet) "EvaporationFlux{12}: input must contain :pet"
    @assert haskey(params, :p1) "EvaporationFlux{12}: params must contain :p1"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[min(1, exp(2 * (1 - input.S / params.p1))) * input.pet])
end

EvaporationFlux(::Val{13}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :pet) "EvaporationFlux{13}: input must contain :pet"
    @assert haskey(params, :p1) "EvaporationFlux{13}: params must contain :p1"
    @assert haskey(params, :p2) "EvaporationFlux{13}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[(params.p1^params.p2) * input.pet])
end

EvaporationFlux(::Val{14}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S) "EvaporationFlux{14}: input must contain :S"
    @assert haskey(input, :pet) "EvaporationFlux{14}: input must contain :pet"
    @assert haskey(params, :p1) "EvaporationFlux{14}: params must contain :p1"
    @assert haskey(params, :p2) "EvaporationFlux{14}: params must contain :p2"
    @assert haskey(params, :Smin) "EvaporationFlux{14}: params must contain :Smin"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[ifelse(input.S > params.Smin, 0.0, (params.p1^params.p2) * input.pet)])
end

EvaporationFlux(::Val{15}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S1) "EvaporationFlux{15}: input must contain :S1"
    @assert haskey(input, :S2) "EvaporationFlux{15}: input must contain :S2"
    @assert haskey(input, :pet) "EvaporationFlux{15}: input must contain :pet"
    @assert haskey(params, :p1) "EvaporationFlux{15}: params must contain :p1"
    @assert haskey(params, :Smax) "EvaporationFlux{15}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[ifelse(input.S2 < params.p1, input.S1 / params.Smax * input.pet, 0.0)])
end

EvaporationFlux(::Val{16}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S) "EvaporationFlux{16}: input must contain :S"
    @assert haskey(input, :pet) "EvaporationFlux{16}: input must contain :pet"
    @assert haskey(params, :p1) "EvaporationFlux{16}: params must contain :p1"
    @assert haskey(params, :p2) "EvaporationFlux{16}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[ifelse(input.S < params.p2, params.p1 * input.pet, 0.0)])
end

EvaporationFlux(::Val{17}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S) "EvaporationFlux{17}: input must contain :S"
    @assert haskey(input, :pet) "EvaporationFlux{17}: input must contain :pet"
    @assert haskey(params, :p1) "EvaporationFlux{17}: params must contain :p1"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[1 / (1 + exp(-params.p1 * input.S)) * input.pet])
end

EvaporationFlux(::Val{18}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S) "EvaporationFlux{18}: input must contain :S"
    @assert haskey(input, :pet) "EvaporationFlux{18}: input must contain :pet"
    @assert haskey(params, :p1) "EvaporationFlux{18}: params must contain :p1"
    @assert haskey(params, :p2) "EvaporationFlux{18}: params must contain :p2"
    @assert haskey(params, :p3) "EvaporationFlux{18}: params must contain :p3"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * exp(-params.p2 * input.S / params.p3) * input.pet])
end

EvaporationFlux(::Val{19}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S) "EvaporationFlux{19}: input must contain :S"
    @assert haskey(input, :pet) "EvaporationFlux{19}: input must contain :pet"
    @assert haskey(params, :p1) "EvaporationFlux{19}: params must contain :p1"
    @assert haskey(params, :p2) "EvaporationFlux{19}: params must contain :p2"
    @assert haskey(params, :Smax) "EvaporationFlux{19}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * (input.S / params.Smax)^params.p2 * input.pet])
end

EvaporationFlux(::Val{20}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S) "EvaporationFlux{20}: input must contain :S"
    @assert haskey(input, :pet) "EvaporationFlux{20}: input must contain :pet"
    @assert haskey(params, :p1) "EvaporationFlux{20}: params must contain :p1"
    @assert haskey(params, :p2) "EvaporationFlux{20}: params must contain :p2"
    @assert haskey(params, :Smax) "EvaporationFlux{20}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[ifelse(input.S < params.p2 * params.Smax, params.p1 * input.S / (params.Smax * params.p2), input.pet)])
end

EvaporationFlux(::Val{21}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S) "EvaporationFlux{21}: input must contain :S"
    @assert haskey(input, :pet) "EvaporationFlux{21}: input must contain :pet"
    @assert haskey(params, :p1) "EvaporationFlux{21}: params must contain :p1"
    @assert haskey(params, :p2) "EvaporationFlux{21}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[ifelse(input.S < params.p1, input.pet,
            ifelse(input.S < params.p2 * params.p1, input.S / params.p1 * input.pet, params.p2 * input.pet))
        ])
end

EvaporationFlux(::Val{22}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S) "EvaporationFlux{22}: input must contain :S"
    @assert haskey(input, :pet) "EvaporationFlux{22}: input must contain :pet"
    @assert haskey(params, :p1) "EvaporationFlux{22}: params must contain :p1"
    @assert haskey(params, :p2) "EvaporationFlux{22}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[ifelse(input.S < params.p1, input.pet,
            ifelse(input.S < params.p2 * params.p1, (input.S - params.p1) / (params.p1 - params.p2) * input.pet, 0))
        ])
end

EvaporationFlux(::Val{23}; input::NamedTuple, params::NamedTuple, output=evaporation) = begin
    @assert haskey(input, :S) "EvaporationFlux{23}: input must contain :S"
    @assert haskey(input, :pet) "EvaporationFlux{23}: input must contain :pet"
    @assert haskey(params, :p1) "EvaporationFlux{23}: params must contain :p1"
    @assert haskey(params, :p2) "EvaporationFlux{23}: params must contain :p2"
    @assert haskey(params, :Smax) "EvaporationFlux{23}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[min(input.S / params.p2 * params.Smax, 1.0) * params.p1 * input.pet + (1 + params.p1) * input.S / params.Smax * input.pet])
end

export EvaporationFlux