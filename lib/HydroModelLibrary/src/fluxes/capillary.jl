"""
    CapillaryFlux

A collection of capillary flux functions for hydrological modeling.

# Available Models

## Scaled Linear (`:scaled`)
```math
q = p_1\\left(1 - \\frac{S}{S_{max}}\\right)
```
Linear decrease of capillary flux with relative soil moisture.
- `p1`: Maximum capillary flux rate [L/T]
- `Smax`: Maximum soil moisture storage [L]

## Constant (`:constant`)
```math
q = \\begin{cases}
p_1 & \\text{if } S > 0 \\\\
0 & \\text{otherwise}
\\end{cases}
```
Constant capillary flux when soil moisture is available.
- `p1`: Constant capillary flux rate [L/T]

## Threshold (`:threshold`)
```math
q = \\begin{cases}
p_1\\left(1 - \\frac{S}{p_2}\\right) & \\text{if } S < p_2 \\\\
0 & \\text{otherwise}
\\end{cases}
```
Linear decrease with threshold cutoff.
- `p1`: Maximum capillary flux rate [L/T]
- `p2`: Threshold soil moisture storage [L]
"""

@variables capillary

CapillaryFlux(::Val{:scaled}; input::NamedTuple, params::NamedTuple, output=capillary) = begin
    @assert haskey(input, :S) "CapillaryFlux{:scaled}: input must contain :S"
    @assert haskey(params, :p1) "CapillaryFlux{:scaled}: params must contain :p1"
    @assert haskey(params, :Smax) "CapillaryFlux{:scaled}: params must contain :Smax"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps, exprs=[params.p1 * (1 - input.S / params.Smax)])
end

CapillaryFlux(::Val{:constant}; input::NamedTuple, params::NamedTuple, output=capillary) = begin
    @assert haskey(input, :S) "CapillaryFlux{:constant}: input must contain :S"
    @assert haskey(params, :p1) "CapillaryFlux{:constant}: params must contain :p1"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps, exprs=[ifelse(input.S > 0, params.p1, 0)])
end

CapillaryFlux(::Val{:threshold}; input::NamedTuple, params::NamedTuple, output=capillary) = begin
    @assert haskey(input, :S) "CapillaryFlux{:threshold}: input must contain :S"
    @assert haskey(params, :p1) "CapillaryFlux{:threshold}: params must contain :p1"
    @assert haskey(params, :p2) "CapillaryFlux{:threshold}: params must contain :p2"
    var, ps = split_vars_and_params(input, params)
    HydroFlux(var => [output], ps, exprs=[ifelse(input.S > params.p2, 0, params.p1 * (1 - input.S / params.p2))])
end

export CapillaryFlux