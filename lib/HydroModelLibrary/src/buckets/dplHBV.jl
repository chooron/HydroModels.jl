@reexport module dPLHBV

using ..HydroModels
using ..Lux

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

@variables soilwater snowpack meltwater suz slz
@variables prcp pet temp
@variables rainfall snowfall melt refreeze infil excess recharge evap q0 q1 q2 q perc
@parameters TT CFMAX CFR CWH LP FC PPERC UZL k0 k1 k2 kp
@variables BETA GAMMA

function LSTMCompact(in_dims, hidden_dims, out_dims)
    lstm_cell = LSTMCell(in_dims => hidden_dims)
    classifier = Dense(hidden_dims => out_dims, sigmoid)
    return @compact(; lstm_cell, classifier) do x::AbstractArray{T,2} where {T}
        x = reshape(x, size(x)..., 1)
        x_init, x_rest = Iterators.peel(LuxOps.eachslice(x, Val(2)))
        y, carry = lstm_cell(x_init)
        output = [vec(classifier(y))]
        for x in x_rest
            y, carry = lstm_cell((x, carry))
            output = vcat(output, [vec(classifier(y))])
        end
        @return reduce(hcat, output)
    end
end

#* parameters estimate by NN
params_nn = LSTMCompact(3, 10, 2)
nn_wrapper = HydroModels.NeuralWrapper([prcp, temp, pet] => [BETA, GAMMA], params_nn, :pnn)

#* snowfall and rainfall split flux
split_flux = HydroFlux([prcp, temp] => [snowfall, rainfall], [TT], exprs=[step_func(TT - temp) * prcp, step_func(temp - TT) * prcp])
#* snowpack bucket
snow_funcs = [
    HydroFlux([temp, snowpack] => [melt], [TT, CFMAX], exprs=[min(snowpack, max(0.0, temp - TT) * CFMAX)]),
    HydroFlux([temp, meltwater] => [refreeze], [TT, CFR, CFMAX], exprs=[min(max((TT - temp), 0.0) * CFR * CFMAX, meltwater)]),
    HydroFlux([meltwater] => [infil], [CWH], exprs=[max(0.0, meltwater - snowpack * CWH)])
]
snow_dfuncs = [StateFlux([snowfall, refreeze] => [melt], snowpack), StateFlux([melt] => [refreeze, infil], meltwater)]
snow_bucket = HydroBucket(name=:hbv_snow, funcs=snow_funcs, dfuncs=snow_dfuncs)

#* soilwater bucket
soil_funcs = [
    HydroFlux([infil, rainfall, BETA] => [recharge], [FC], exprs=[(rainfall + infil) * clamp(max(0.0, soilwater / FC)^(BETA * 5 + 1), 0, 1)]),
    HydroFlux([soilwater] => [excess], [FC], exprs=[max(soilwater - FC, 0.0)]),
    HydroFlux([soilwater, pet, GAMMA] => [evap], [LP, FC], exprs=[clamp(max(0.0, soilwater / (LP * FC))^(GAMMA + 1), 0, 1) * pet]),
]
soil_dfuncs = [StateFlux([rainfall, infil] => [recharge, excess, evap], soilwater)]
soil_bucket = HydroBucket(name=:hbv_soil, funcs=soil_funcs, dfuncs=soil_dfuncs)

#* up and low zone bucket
zone_funcs = [
    HydroFlux([suz] => [perc], [PPERC], exprs=[suz * PPERC]),
    HydroFlux([suz] => [q0], [UZL, k0], exprs=[max(0.0, suz - UZL) * k0]),
    HydroFlux([suz] => [q1], [k1], exprs=[suz * k1]),
    HydroFlux([slz] => [q2], [k2], exprs=[slz * k2]),
    HydroFlux([q0, q1, q2] => [q], exprs=[q0 + q1 + q2]),
]
zone_dfuncs = [StateFlux([recharge, excess] => [perc, q0, q1], suz), StateFlux([perc] => [q2], slz),]
zone_bucket = HydroBucket(name=:hbv_zone, funcs=zone_funcs, dfuncs=zone_dfuncs)
dpl_hbv_model = HydroModel(name=:dpl_hbv, components=[nn_wrapper, split_flux, snow_bucket, soil_bucket, zone_bucket])

export nn_wrapper, split_flux, snow_bucket, soil_bucket, zone_bucket, dpl_hbv_model
end