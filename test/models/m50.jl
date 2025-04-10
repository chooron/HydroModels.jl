step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

# Model parameters
# Physical process parameters
@parameters Tmin Tmax Df Smax f Qmax
# Normalization parameters
@parameters snowpack_std snowpack_mean
@parameters soilwater_std soilwater_mean
@parameters prcp_std prcp_mean
@parameters temp_std temp_mean

# Model variables
# Input variables
@variables prcp temp lday
# State variables
@variables snowpack soilwater
# Process variables
@variables pet rainfall snowfall melt
# Neural network variables
@variables log_evap_div_lday log_flow flow
@variables norm_snw norm_slw norm_temp norm_prcp

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

# Snow component
snow_bucket = @hydrobucket :m50_snow begin
    fluxes = begin
        @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
        @hydroflux begin
            snowfall ~ step_func(Tmin - temp) * prcp
            rainfall ~ step_func(temp - Tmin) * prcp
        end
        @hydroflux melt ~ step_func(temp - Tmax) * min(snowpack, Df * (temp - Tmax))
    end
    dfluxes = begin
        @stateflux snowpack ~ snowfall - melt
    end
end

# Neural network fluxes
ep_nn_flux = NeuralFlux([norm_snw, norm_slw, norm_temp] => [log_evap_div_lday], ep_nn)
q_nn_flux = NeuralFlux([norm_slw, norm_prcp] => [log_flow], q_nn)

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

# Flow conversion
flow_conversion = @hydroflux flow ~ exp(log_flow)

# Complete model
m50_model = @hydromodel :m50 begin
    snow_bucket
    soil_bucket
    flow_conversion
end

export m50_model, snow_bucket, soil_bucket
