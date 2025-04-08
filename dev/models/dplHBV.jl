using Lux

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

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

@variables soilwater snowpack meltwater suz slz
@variables prcp pet temp
@variables rainfall snowfall melt refreeze infil excess recharge evap q0 q1 q2 q perc
@parameters TT CFMAX CFR CWH LP FC PPERC UZL k0 k1 k2 kp
@variables BETA GAMMA
#* parameters estimate by NN
params_nn = LSTMCompact(3, 10, 2)
nn_wrapper = HydroModels.NeuralWrapper([prcp, temp, pet] => [BETA, GAMMA], params_nn, name=:pnn)

#* snowfall and rainfall split flux
split_flux = HydroFlux([prcp, temp] => [snowfall, rainfall], [TT],
    exprs=[step_func(TT - temp) * prcp, step_func(temp - TT) * prcp])
#* snowpack bucket
snow_fluxes = [
    HydroFlux([temp, snowpack] => [melt], [TT, CFMAX], exprs=[min(snowpack, max(0.0, temp - TT) * CFMAX)]),
    HydroFlux([temp, meltwater] => [refreeze], [TT, CFR, CFMAX], exprs=[min(max((TT - temp), 0.0) * CFR * CFMAX, meltwater)]),
    HydroFlux([meltwater] => [infil], [CWH], exprs=[max(0.0, meltwater - snowpack * CWH)])
]
snow_dfluxes = [StateFlux([snowfall, refreeze] => [melt], snowpack), StateFlux([melt] => [refreeze, infil], meltwater)]
snow_bucket = HydroBucket(name=:hbv_snow, fluxes=snow_fluxes, dfluxes=snow_dfluxes)

#* soilwater bucket
soil_fluxes = [
    HydroFlux([infil, rainfall, BETA] => [recharge], [FC], exprs=[(rainfall + infil) * clamp(max(0.0, soilwater / FC)^(BETA * 5 + 1), 0, 1)]),
    HydroFlux([soilwater] => [excess], [FC], exprs=[max(soilwater - FC, 0.0)]),
    HydroFlux([soilwater, pet, GAMMA] => [evap], [LP, FC], exprs=[clamp(max(0.0, soilwater / (LP * FC))^(GAMMA + 1), 0, 1) * pet]),
]
soil_dfluxes = [StateFlux([rainfall, infil] => [recharge, excess, evap], soilwater)]
soil_bucket = HydroBucket(name=:hbv_soil, fluxes=soil_fluxes, dfluxes=soil_dfluxes)

#* up and low zone bucket
zone_fluxes = [
    HydroFlux([suz] => [perc], [PPERC], exprs=[suz * PPERC]),
    HydroFlux([suz] => [q0], [UZL, k0], exprs=[max(0.0, suz - UZL) * k0]),
    HydroFlux([suz] => [q1], [k1], exprs=[suz * k1]),
    HydroFlux([slz] => [q2], [k2], exprs=[slz * k2]),
    HydroFlux([q0, q1, q2] => [q], exprs=[q0 + q1 + q2]),
]
zone_dfluxes = [StateFlux([recharge, excess] => [perc, q0, q1], suz), StateFlux([perc] => [q2], slz),]
zone_bucket = HydroBucket(name=:hbv_zone, fluxes=zone_fluxes, dfluxes=zone_dfluxes)
dpl_hbv_model = HydroModel(name=:dpl_hbv, components=[nn_wrapper, split_flux, snow_bucket, soil_bucket, zone_bucket])

