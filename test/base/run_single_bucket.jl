Df_v, Tmax_v, Tmin_v = 2.674548848, 0.175739196, -2.092959084
params = ComponentVector(Df=Df_v, Tmax=Tmax_v, Tmin=Tmin_v)
init_states = ComponentVector(snowpack=0.0)

ts = collect(1:10)
df = DataFrame(CSV.File("../data/exphydro/01013500.csv"))
input_ntp = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])

input = Matrix(reduce(hcat, collect(input_ntp[[:temp, :lday, :prcp]]))')
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

@variables temp lday prcp pet snowfall rainfall melt snowpack
@parameters Tmin Tmax Df


@testset "test single hydro element (with state)" begin

    snow_single_ele = @hydrobucket :snow begin
        fluxes = begin
            @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
            @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
            @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
            @hydroflux melt ~ step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))
        end
        dfluxes = begin
            @stateflux snowpack ~ snowfall - melt
        end
    end

    @testset "test hydro element info" begin
        @test Set(HydroModels.get_input_names(snow_single_ele)) == Set([:temp, :lday, :prcp])
        @test Set(HydroModels.get_param_names(snow_single_ele)) == Set([:Tmin, :Tmax, :Df])
        @test Set(HydroModels.get_output_names(snow_single_ele)) == Set([:pet, :snowfall, :rainfall, :melt])
        @test Set(HydroModels.get_state_names(snow_single_ele)) == Set([:snowpack])
    end

    @testset "test run with single node input" begin
        pas = ComponentVector(params=ComponentVector(Df=Df_v, Tmax=Tmax_v, Tmin=Tmin_v))
        result = snow_single_ele(input, pas; initstates=init_states, timeidx=ts, solver=ManualSolver(mutable=true))
        ele_state_and_output_names = vcat(HydroModels.get_state_names(snow_single_ele), HydroModels.get_output_names(snow_single_ele))
        result = NamedTuple{Tuple(ele_state_and_output_names)}(eachslice(result, dims=1))
    end
end

@testset "test single hydro element (without state)" begin
    @variables temp lday prcp pet snowfall rainfall melt snowpack
    @parameters Tmin Tmax Df

    snow_single_ele = @hydrobucket :snow begin
        fluxes = begin
            @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
            @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
            @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
        end
    end

    @testset "test hydro element info" begin
        @test Set(HydroModels.get_input_names(snow_single_ele)) == Set((:temp, :lday, :prcp))
        @test Set(HydroModels.get_param_names(snow_single_ele)) == Set((:Tmin,))
        @test Set(HydroModels.get_output_names(snow_single_ele)) == Set((:pet, :snowfall, :rainfall))
        @test Set(HydroModels.get_state_names(snow_single_ele)) == Set()
    end

    pas = ComponentVector(params=(Df=Df_v, Tmax=Tmax_v, Tmin=Tmin_v))
    @testset "test run with single node input" begin
        result = snow_single_ele(input, pas; initstates=init_states, timeidx=ts, solver=ManualSolver(mutable=true))
        ele_state_and_output_names = vcat(HydroModels.get_state_names(snow_single_ele), HydroModels.get_output_names(snow_single_ele))
        result = NamedTuple{Tuple(ele_state_and_output_names)}(eachslice(result, dims=1))
    end
end