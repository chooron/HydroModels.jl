module GammaUnitHydro
using ..HydroModels
using ..HydroModels: AbstractRoute

function gamma_uh()
    uh = HydroModels.@unithydro :maxbas_uh begin
        uh_func = begin
            10 => (1 / t) * (beta2 * t)^alpha2 / gamma(alpha2) * exp(-beta2 * t)
        end
        uh_vars = [flow => flow_lag]
        configs = (solvetype=:SPARSE, )
    end
    return uh
end

end