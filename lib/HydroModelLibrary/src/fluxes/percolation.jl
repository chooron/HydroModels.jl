@variables percolation

PercolationFlux(::Val{:1}; input::NamedTuple, params::NamedTuple, output=percolation) = begin
    @assert haskey(input, :S) "PercolationFlux{:1}: input must contain :S"
    @assert haskey(params, :p1) "PercolationFlux{:1}: params must contain :p1"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[ifelse(input.S > 0, params.p1, 0)]
    )
end

PercolationFlux(::Val{:2}; input::NamedTuple, params::NamedTuple, output=percolation) = begin
    @assert haskey(input, :S) "PercolationFlux{:2}: input must contain :S"
    @assert haskey(params, :p1) "PercolationFlux{:2}: params must contain :p1"
    @assert haskey(params, :Smax) "PercolationFlux{:2}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * input.S / params.Smax]
    )
end

PercolationFlux(::Val{:3}; input::NamedTuple, params::NamedTuple, output=percolation) = begin
    @assert haskey(input, :S) "PercolationFlux{:3}: input must contain :S"
    @assert haskey(params, :p1) "PercolationFlux{:3}: params must contain :p1"
    @assert haskey(params, :Smax) "PercolationFlux{:3}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.Smax^(-4) / 4 * (4 / 9)^(-4) * input.S^5]
    )
end

# TODO Demand-based percolation scaled by available moisture 

PercolationFlux(::Val{:5}; input::NamedTuple, params::NamedTuple, output=percolation) = begin
    @assert haskey(input, :S) "PercolationFlux{:5}: input must contain :S"
    @assert haskey(params, :p1) "PercolationFlux{:5}: params must contain :p1"
    @assert haskey(params, :p2) "PercolationFlux{:5}: params must contain :p2"
    @assert haskey(params, :Smax) "PercolationFlux{:5}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[params.p1 * (input.S / params.Smax)^params.p2]
    )
end

PercolationFlux(::Val{:6}; input::NamedTuple, params::NamedTuple, output=percolation) = begin
    @assert haskey(input, :S) "PercolationFlux{:6}: input must contain :S"
    @assert haskey(params, :p1) "PercolationFlux{:6}: params must contain :p1"
    @assert haskey(params, :p2) "PercolationFlux{:6}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps,
        exprs=[min(1, max(0, input.S) / params.p2) * params.p1]
    )
end

export PercolationFlux