step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5
# define variables and parameters
@variables temp lday pet prcp snowfall rainfall snowpack melt
@variables soilwater pet evap baseflow surfaceflow flow rainfall
@parameters Tmin Tmax Df Smax Qmax f

# define model components
bucket_1 = @hydrobucket :surface begin
    fluxes = begin
        @hydroflux begin
            snowfall ~ step_func(Tmin - temp) * prcp
            rainfall ~ step_func(temp - Tmin) * prcp
        end
        @hydroflux melt ~ step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))
        @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
    end
    dfluxes = begin
        @stateflux snowpack ~ snowfall - melt
    end
end

bucket_2 = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux evap ~ step_func(soilwater) * pet * min(1.0, soilwater / Smax)
        @hydroflux baseflow ~ step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))
        @hydroflux surfaceflow ~ max(0.0, soilwater - Smax)
        @hydroflux flow ~ baseflow + surfaceflow
    end
    dfluxes = begin
        @stateflux soilwater ~ (rainfall + melt) - (evap + flow)
    end
end

exphydro_model = @hydromodel :exphydro begin
    bucket_1
    bucket_2
end

export exphydro_model, bucket_1, bucket_2
