function MinMaxNormFlux(
    fluxes::Pair{Vector{Num},Vector{Num}},
    params::Vector{Vector{Num}}
)
    flux_exprs = [(input - param[2]) / (param[1] - param[2]) for (input, param) in zip(fluxes[1], params)]
    params_vec = reduce(vcat, params)

    HydroFlux(
        fluxes,
        params_vec,
        exprs=flux_exprs,
    )
end

function StdMeanNormFlux(
    fluxes::Pair{Vector{Num},Vector{Num}},
    params::Vector{Vector{Num}}
)
    flux_exprs = [(input - param[1]) / param[2] for (input, param) in zip(fluxes[1], params)]
    params_vec = reduce(vcat, params)

    HydroFlux(
        fluxes,
        params_vec,
        exprs=flux_exprs,
    )
end

function TranparentFlux(
    flux_names::Pair{Vector{Symbol},Vector{Symbol}},
)
    old_names, new_names = flux_names[1], flux_names[2]
    inputs = [@variables $nm = 0.0 for nm in old_names]
    outputs = [@variables $nm = 0.0 for nm in new_names]
    flux_exprs = [var for var in inputs]

    HydroFlux(inputs => outputs, exprs=flux_exprs)
end

function RenameFlux(
    fluxes::Pair{Vector{Num},Vector{Num}},
)
    flux_exprs = [var for var in fluxes[1]]
    HydroFlux(fluxes, Num[], exprs=flux_exprs)
end

export MinMaxNormFlux, StdMeanNormFlux, TranparentFlux, RenameFlux