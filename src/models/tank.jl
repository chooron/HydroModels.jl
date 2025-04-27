module tank
using ..HydroModels
using ..HydroModels: step_func

# Model variables
@variables S1 [description = "current storage in the upper one", unit = "mm"]
@variables P [description = "precipitation", unit = "mm/d"]
@variables E1 [description = "evaporation", unit = "mm/d"]
@variables E12 [description = "drainage", unit = "mm/d"]
@variables Y1 [description = "surface runoff", unit = "mm/d"]
@variables Y2 [description = "surface runoff", unit = "mm/d"]
@variables Ep [description = "potential rate", unit = "mm/d"]

@variables S2 [description = "current storage in the intermediate zone", unit = "mm"]
@variables F12 [description = "drainage", unit = "mm/d"]
@variables E2 [description = "evaporation", unit = "mm/d"]
@variables F23 [description = "drainage", unit = "mm/d"]
@variables Y3 [description = "intermediate discharge", unit = "mm/d"]

@variables S3 [description = "current storage in the sub-base zone", unit = "mm"]
@variables E3 [description = "evaporation", unit = "mm/d"]
@variables F34 [description = "drainage", unit = "mm/d"]
@variables Y4 [description = "ub-base runoff", unit = "mm/d"]

@variables S4 [description = "current storage in the base layer", unit = "mm"]
@variables E4 [description = "evaporation", unit = "mm/d"]
@variables Y5 [description = "baseflow", unit = "mm/d"]

@variables Qt [description = "total runoff", unit = "mm/d"]

# Model parameters
@parameters a0 [description = "Time parameter for drainage 1>2", bounds = (0, 1), unit = "1/d"]
@parameters b0 [description = "Time parameter for drainage 2>3", bounds = (0, 1), unit = "1/d"]
@parameters c0 [description = "Time parameter for drainage 3>4", bounds = (0, 1), unit = "1/d"]
@parameters a1 [description = "Time parameter for surface runoff 1", bounds = (0, 1), unit = "1/d"]
@parameters fa [description = "Fraction of a1 that is a2", bounds = (0, 1), unit = "-"]
@parameters fb [description = "Fraction of a2 that is b1", bounds = (0, 1), unit = "-"]
@parameters fc [description = "Fraction of b1 that is c1", bounds = (0, 1), unit = "-"]
@parameters fd [description = "Fraction of c1 that is d1", bounds = (0, 1), unit = "-"]
@parameters st [description = "Maximum soil depth (sum of runoff thresholds) affected the results!", bounds = (1, 2000), unit = "mm"]
@parameters f2 [description = "Fraction of st that consitutes threshold t2 ", bounds = (0.01, 0.99), unit = "-"]
@parameters f1 [description = "Fraction of st-t2 that is added to t2 to find threshold 1", bounds = (0.01, 0.99), unit = "-"]
@parameters f3 [description = "Fraction of st-t1-t2 that consitutes threshold 3", bounds = (0.01, 0.99), unit = "-"]

t2 = f2 * st         # Threshold before surface runoff 2 starts [mm]
t1 = t2 + f1 * (st - t2) # Threshold before surface runoff 1 starts [mm]
t3 = f3 * (st - t1)    # Threshold before intermediate runoff starts [mm]
t4 = st - t1 - t3      # Threshold before sub-base runoff starts [mm]
a2 = fa * a1         # Time parameter for surface runoff 2 [d-1]
b1 = fb * a2         # Time parameter for intermediate runoff 1 [d-1]
c1 = fc * b1         # Time parameter for sub-base runoff 1 [d-1]
d1 = fd * c1         # Time parameter for base runoff 1 [d-1]

bucket1 = @hydrobucket :bucket1 begin
    fluxes = begin
        @hydroflux E1 ~ min(max(0.0, S1), Ep)
        @hydroflux F12 ~ a0 * max(0.0, S1)
        @hydroflux Y2 ~ a2 * max(0.0, S1 - t2)
        @hydroflux Y1 ~ a1 * max(0.0, S1 - t1)
    end
    dfluxes = begin
        @stateflux S1 ~ P - E1 - F12 - Y2 - Y1
    end
end

bucket2 = @hydrobucket :bucket2 begin
    fluxes = begin
        @hydroflux E2 ~ min(max(0.0, S2), Ep - E1)
        @hydroflux F23 ~ b0 * max(0.0, S2)
        @hydroflux Y3 ~ b1 * max(0.0, S2 - t3)
    end
    dfluxes = begin
        @stateflux S2 ~ F12 - E2 - F23 - Y3
    end
end

bucket3 = @hydrobucket :bucket3 begin
    fluxes = begin
        @hydroflux E3 ~ min(max(0.0, S3), Ep - E1 - E2)
        @hydroflux F34 ~ c0 * max(0.0, S3)
        @hydroflux Y4 ~ c1 * max(0.0, S3 - t4)
    end
    dfluxes = begin
        @stateflux S3 ~ F23 - E3 - F34 - Y4
    end
end

bucket4 = @hydrobucket :bucket4 begin
    fluxes = begin
        @hydroflux E4 ~ min(max(0.0, S4), Ep - E1 - E2 - E3)
        @hydroflux Y5 ~ d1 * max(0.0, S4)
    end
    dfluxes = begin
        @stateflux S4 ~ F34 - E4 - Y5
    end
end

flux1 = @hydroflux Qt ~ Y1 + Y2 + Y3 + Y4 + Y5

model = @hydromodel :tank begin
    bucket1
    bucket2
    bucket3
    bucket4
    flux1
end

end