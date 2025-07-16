module UnitHydro

function CONVOL_GAMMA(;
)
    uh1 = HydroModels.@unithydro :maxbas_uh begin
        uh_func = begin
            lag => (t / lag)^2.5
        end
        uh_vars = [q1]
        configs = (solvetype=:SPARSE, suffix=:_lag)
    end

end

end