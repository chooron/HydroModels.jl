# Neural Network Integration in Hydrological Modeling

A distinctive feature of the HydroModels.jl framework is its seamless integration of neural networks into hydrological model computations through the `NeuralFlux` type. This integration enables hybrid modeling approaches that combine process-based equations with data-driven components.

## NeuralFlux Architecture

The `NeuralFlux` type is a subtype of `AbstractNeuralFlux` in the HydroModels.jl framework, maintaining interface consistency with `HydroFlux` for seamless integration with `HydroBucket` components:

```julia
AbstractNeuralFlux <: AbstractFlux
```

Constructing a `NeuralFlux` requires three essential components:

- `chain`: A neural network model constructed using Lux.jl
- `inputs`: Input variables defined using ModelingToolkit.jl's symbolic system
- `outputs`: Output variables defined using ModelingToolkit.jl's symbolic system

The `@neuralflux` macro provides a concise syntax for defining neural network fluxes:

```julia
@neuralflux [output1, output2, ...] ~ chain([input1, input2, ...])
```

## Implementation Methodology

### Standalone Usage

`NeuralFlux` can be used independently to perform computations, functionally equivalent to direct Lux.jl neural network invocation but with the advantage of accepting `ComponentVector` parameters, enabling integration with `Element` and `Model` components:

```julia
using HydroModels, Lux

@variables i1 i2 i3 o1 o2
chain_nm = :testnn
chain = Lux.Chain(
    Lux.Dense(3 => 16, tanh),
    Lux.Dense(16 => 2, leakyrelu),
    name=chain_nm
)
neuralflux = @neuralflux [o1, o2] ~ chain([i1, i2, i3])
input = rand(3, 100)
nnps = ComponentVector(Lux.initialparameters(StableRNG(42), chain)) |> Vector
output = neuralflux(input, ComponentVector(nns=(chain_nm=nnps,)))
```

During computation, `NeuralFlux` accepts a `ComponentVector` containing neural network parameters under the `nns` namespace. The `chain_nm` corresponds to the name assigned during chain definition, with each neural network's parameters represented as a vector. During computation, automatic conversion occurs through the recorded `ComponentVector` axes.

## Application Paradigms

When used with `HydroBucket`, `NeuralFlux` functions as a flux component alongside `HydroFlux` and `StateFlux`. Two primary application paradigms have emerged: process substitution and parameter estimation.

### 1. Process Substitution

Replacing conventional hydrological formulations with neural networks represents a fundamental integration approach. The general formulation can be expressed as:

```math
flux = \text{NeuralNetwork}(\text{variables}; \theta)
```

This substitution can be implemented at varying levels of complexity:

#### Post-Processing Integration

Neural networks can integrate all hydrological fluxes to predict streamflow (functioning as a post-processor, e.g., XAJ-LSTM):

```julia
@neuralflux flow ~ lstm_model([soilwater, infiltration, prcp, ...])
```

This implementation constructs an LSTM model that predicts streamflow using intermediate states (soil moisture, infiltration, etc.) calculated by a conventional hydrological model.

#### Differential Equation Integration

A more complex integration involves neural networks participating in ordinary differential equation (ODE) calculations, where neural network inputs include state variables and outputs contribute to state variable updates. The M50 model exemplifies this approach:

```julia
@parameters snowpack_std snowpack_mean
@parameters soilwater_std soilwater_mean
@parameters prcp_std prcp_mean
@parameters temp_std temp_mean

# Input variables
@variables prcp temp lday
# State variables
@variables snowpack soilwater
# Process variables
@variables pet rainfall snowfall melt
# Neural network variables
@variables log_evap_div_lday log_flow flow
@variables norm_snw norm_slw norm_temp norm_prcp

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

# Neural network definitions
ep_nn = Lux.Chain(
    Lux.Dense(3 => 16, tanh),
    Lux.Dense(16 => 16, leakyrelu),
    Lux.Dense(16 => 1, leakyrelu),
    name=:epnn
)

q_nn = Lux.Chain(
    Lux.Dense(2 => 16, tanh),
    Lux.Dense(16 => 16, leakyrelu),
    Lux.Dense(16 => 1, leakyrelu),
    name=:qnn
)

# Soil water component
soil_bucket = @hydrobucket :m50_soil begin
    fluxes = begin
        @hydroflux norm_snw ~ (snowpack - snowpack_mean) / snowpack_std
        @hydroflux norm_slw ~ (soilwater - soilwater_mean) / soilwater_std
        @hydroflux norm_prcp ~ (prcp - prcp_mean) / prcp_std
        @hydroflux norm_temp ~ (temp - temp_mean) / temp_std
        @neuralflux log_evap_div_lday ~ ep_nn([norm_snw, norm_slw, norm_temp])
        @neuralflux log_flow ~ q_nn([norm_slw, norm_prcp])
    end
    dfluxes = begin
        @stateflux soilwater ~ rainfall + melt - step_func(soilwater) * lday * exp(log_evap_div_lday) - step_func(soilwater) * exp(log_flow)
    end
end
```

This implementation demonstrates the M50 model replacing evaporation and runoff calculations in the ExpHydro model. Two neural networks are defined using Lux.jl for predicting evaporation and runoff. The `HydroBucket` first defines normalization calculations for snowpack, soil water, precipitation, and temperature, then constructs two `NeuralFlux` components based on these neural networks. In the `dfluxes` section, `log_evap_div_lday` and `log_flow` participate in soil state updates, indicating neural network integration into the hydrological model solution (this substitution approach incurs significantly higher computational costs than post-processing integration).

### 2. Parameter Estimation

Using neural networks to estimate hydrological model parameters represents another significant application. The dPL-HBV model exemplifies this approach, using an LSTM model to predict sensitive HBV model parameters (BETA and GAMMA) based on forcing data and watershed static attributes, enabling dynamic parameter estimation to enhance effective rainfall and evaporation calculations:

```julia
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
params_nn_flux = NeuralFlux([prcp, temp, pet] => [BETA, GAMMA], params_nn, name=:pnn)

#* snowfall and rainfall split flux
split_flux = @hydroflux begin
    snowfall ~ step_func(TT - temp) * prcp
    rainfall ~ step_func(temp - TT) * prcp
end

snow_bucket = @hydrobucket :hbv_snow begin
    fluxes = begin
        @hydroflux melt ~ min(snowpack, max(0.0, temp - TT) * CFMAX)
        @hydroflux refreeze ~ min(max((TT - temp), 0.0) * CFR * CFMAX, meltwater)
        @hydroflux infil ~ max(0.0, meltwater - snowpack * CWH)
    end
    dfluxes = begin
        @stateflux snowpack ~ snowfall + refreeze - melt
        @stateflux meltwater ~ melt - refreeze - infil
    end
end

soil_bucket = @hydrobucket :hbv_soil begin
    fluxes = begin
        @hydroflux recharge ~ (rainfall + infil) * clamp((max(soilwater / FC, 0.0))^(BETA * 5 + 1), 0, 1)
        @hydroflux excess ~ max(soilwater - FC, 0.0)
        @hydroflux evap ~ clamp(soilwater / (LP * FC), 0, 1) * pet
    end
    dfluxes = begin
        @stateflux soilwater ~ (rainfall + infil) - (recharge + excess + evap)
    end
end

zone_bucket = @hydrobucket :hbv_zone begin
    fluxes = begin
        @hydroflux perc ~ suz * PPERC
        @hydroflux q0 ~ max(0.0, suz - UZL) * k0
        @hydroflux q1 ~ suz * k1
        @hydroflux q2 ~ slz * k2
        @hydroflux q ~ q0 + q1 + q2
    end
    dfluxes = begin
        @stateflux suz ~ recharge + excess - (perc + q0 + q1)
        @stateflux slz ~ perc - q2
    end
end

dpl_hbv_model = @hydromodel :dpl_hbv begin
    params_nn_flux
    split_flux
    snow_bucket
    soil_bucket
    zone_bucket
end
```

This model incorporates the LSTM-based parameter estimation flux `params_nn_flux` while maintaining consistency with the standard HBV model structure. In `params_nn_flux`, a custom LSTM model predicts GAMMA and BETA parameters using forcing data (precipitation, temperature, daylight duration). These parameters function as dynamic, time-varying dimensionless fluxes, necessitating their definition as `variables` rather than `parameters`. The LSTM model's output activation function is typically sigmoid, requiring denormalization to map parameters to appropriate ranges:

```julia
(rainfall + infil) * clamp((max(soilwater / FC, 0.0))^(BETA * 5 + 1), 0, 1)
```

This expression maps BETA to a range of 1-6.

This parameter estimation approach leverages LSTM to capture cumulative effects of forcing data for predicting GAMMA and BETA parameters. While this approach decouples neural network computation from hydrological model calculation, enhancing computational efficiency, the limited interpretability of LSTM models may contradict the original intent of neural network integration (as discussed in relevant literature).

An alternative integration approach directly uses state variables like soil moisture as neural network inputs to represent cumulative effects:

```julia
ps_nn = Lux.Chain(
    Lux.Dense(4 => 16, tanh),
    Lux.Dense(16 => 16, leakyrelu),
    Lux.Dense(16 => 2, leakyrelu),
    name=:psnn
)

soil_bucket = @hydrobucket :hbv_soil begin
    fluxes = begin
        @neuralflux [BETA, GAMMA] ~ ps_nn([soilwater, prcp, pet, temp])
        @hydroflux recharge ~ (rainfall + infil) * clamp((max(soilwater / FC, 0.0))^(BETA * 5 + 1), 0, 1)
        @hydroflux excess ~ max(soilwater - FC, 0.0)
        @hydroflux evap ~ clamp(soilwater / (LP * FC), 0, 1) * pet
    end
    dfluxes = begin
        @stateflux soilwater ~ (rainfall + infil) - (recharge + excess + evap)
    end
end
```

This implementation constructs a feedforward neural network that predicts BETA and GAMMA based on soil water, precipitation, potential evapotranspiration, and temperature, though at a significant computational cost.

## Summary

### Advancing Hybrid Modeling

The integration of neural networks into hydrological models through HydroModels.jl's `NeuralFlux` component represents a significant advancement in hybrid modeling approaches. This framework enables seamless combination of process-based equations with data-driven components, creating models that leverage the strengths of both paradigms.

### Integration Paradigms

The framework offers two primary integration paradigms:

- **Process Substitution**: Replacing conventional hydrological formulations with neural networks, ranging from simple post-processing to complex differential equation integration
- **Parameter Estimation**: Using neural networks to dynamically adapt model parameters based on environmental conditions

### Technical Implementation

The type-stable implementation and consistent interfaces ensure that neural network integration maintains computational efficiency while preserving the scientific rigor of traditional hydrological modeling approaches. The framework's design allows for seamless integration within the broader HydroModels.jl ecosystem.
