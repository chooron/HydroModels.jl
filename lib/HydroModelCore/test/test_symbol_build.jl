using Test
using HydroModelCore
using Symbolics

struct DummyHydroFlux <: HydroModelCore.AbstractHydroFlux
    name::Symbol
    infos::HydroInfos
    exprs::Vector{Num}
end

struct DummyStateFlux <: HydroModelCore.AbstractStateFlux
    name::Symbol
    infos::HydroInfos
    exprs::Vector{Num}
end

struct DummyNeuralFlux <: HydroModelCore.AbstractNeuralFlux
    name::Symbol
    infos::HydroInfos
    chain
end

function _bench_ns_per_call(f, args...; warmup::Int=20, iters::Int=1000)
    for _ in 1:warmup
        f(args...)
    end
    t0 = time_ns()
    for _ in 1:iters
        f(args...)
    end
    return (time_ns() - t0) / iters
end

@inline _approx_equal(a, b) = a ≈ b
@inline _approx_equal(a::AbstractArray, b::AbstractArray) = all(isapprox.(a, b))

@testset "SymbolicUtils Builders" begin
    @variables temp prcp storage k

    flux_infos = HydroInfos(
        inputs=[:temp, :prcp],
        outputs=[:q],
        params=[:k],
        states=Symbol[],
        nns=Symbol[]
    )
    flux_expr = k * (temp + sin(prcp))
    pas = (params=(k=0.5,),)

    @testset "Scalar Flux" begin
        flux_func = build_symbolic_flux_func([flux_expr], flux_infos)
        result = flux_func([1.0, 2.0], pas)
        @test length(result) == 1
        @test result[1] ≈ 0.5 * (1.0 + sin(2.0))

        repeated_expr = (temp + prcp)^2 + k * (temp + prcp)
        f_nocse = build_symbolic_flux_func([repeated_expr], flux_infos; enable_cse=false)
        f_cse = build_symbolic_flux_func([repeated_expr], flux_infos; enable_cse=true)

        x = [1.2, 3.4]
        y_nocse = f_nocse(x, pas)
        y_cse = f_cse(x, pas)
        @test _approx_equal(y_nocse, y_cse)
    end

    @testset "Vector Flux" begin
        flux_func = build_symbolic_flux_func([flux_expr], flux_infos; dims=1)
        inputs = [1.0 2.0 3.0; 4.0 5.0 6.0]
        result = flux_func(inputs, pas)
        expected = @. 0.5 * (inputs[1, :] + sin(inputs[2, :]))
        @test length(result) == 1
        @test result[1] ≈ expected
    end

    @testset "Matrix Flux" begin
        flux_func = build_symbolic_flux_func([flux_expr], flux_infos; dims=2)
        inputs = Array{Float64}(undef, 2, 2, 2)
        inputs[1, :, :] = [1.0 2.0; 3.0 4.0]
        inputs[2, :, :] = [5.0 6.0; 7.0 8.0]

        result = flux_func(inputs, pas)
        expected = @. 0.5 * (inputs[1, :, :] + sin(inputs[2, :, :]))
        @test length(result) == 1
        @test result[1] ≈ expected
    end

    @testset "Bucket Builder" begin
        bucket_infos = HydroInfos(
            inputs=[:temp],
            outputs=[:q],
            states=[:storage],
            params=[:k],
            nns=Symbol[]
        )
        flux = DummyHydroFlux(
            :surface,
            HydroInfos(inputs=[:temp], outputs=[:q], states=[:storage], params=[:k], nns=Symbol[]),
            Num[k * (temp + storage)]
        )
        dflux = DummyStateFlux(
            :storage_balance,
            HydroInfos(inputs=[:temp], outputs=Symbol[], states=[:storage], params=[:k], nns=Symbol[]),
            Num[temp - storage]
        )

        flux_func, diff_func = build_symbolic_bucket_func([flux], [dflux], bucket_infos, false)
        flux_result = flux_func(reshape([1.0, 2.0, 3.0], 1, :), reshape([4.0, 5.0, 6.0], 1, :), pas)
        diff_result = diff_func([1.0], [4.0], pas)

        @test flux_result[1] ≈ @. 0.5 * ([1.0, 2.0, 3.0] + [4.0, 5.0, 6.0])
        @test diff_result ≈ [-3.0]
    end

    @testset "Neural Flux Rejection" begin
        neural_infos = HydroInfos(inputs=[:temp], outputs=[:q], states=[:storage], params=[:k], nns=[:nn])
        neural_flux = DummyNeuralFlux(:nn_flux, neural_infos, nothing)
        dflux = DummyStateFlux(:storage_balance, neural_infos, Num[temp - storage])

        @test_throws ArgumentError build_symbolic_bucket_func([neural_flux], [dflux], neural_infos, false)
    end

    @testset "Old vs New Consistency" begin
        @variables temp2 storage2 k2

        bucket_infos = HydroInfos(
            inputs=[:temp2],
            outputs=[:q2],
            states=[:storage2],
            params=[:k2],
            nns=Symbol[]
        )

        flux = DummyHydroFlux(
            :surface,
            HydroInfos(inputs=[:temp2], outputs=[:q2], states=[:storage2], params=[:k2], nns=Symbol[]),
            Num[(temp2 + storage2)^2 + k2 * (temp2 + storage2)]
        )
        dflux = DummyStateFlux(
            :storage_balance,
            HydroInfos(inputs=[:temp2], outputs=Symbol[], states=[:storage2], params=[:k2], nns=Symbol[]),
            Num[temp2 - storage2]
        )

        pas2 = (params=(k2=0.3,),)

        # Vector case (multiply = false)
        old_flux_v, old_diff_v = build_bucket_func([flux], [dflux], bucket_infos, false)
        new_flux_v, new_diff_v = build_symbolic_bucket_func([flux], [dflux], bucket_infos, false; enable_cse=true)

        in_v = reshape([1.0, 2.0, 3.0], 1, :)
        st_v = reshape([4.0, 5.0, 6.0], 1, :)
        old_v = old_flux_v(in_v, st_v, pas2)
        new_v = new_flux_v(in_v, st_v, pas2)
        @test _approx_equal(old_v[1], new_v[1])

        old_dv = old_diff_v([1.5], [0.2], pas2)
        new_dv = new_diff_v([1.5], [0.2], pas2)
        @test _approx_equal(old_dv, new_dv)

        # Matrix case (multiply = true)
        old_flux_m, old_diff_m = build_bucket_func([flux], [dflux], bucket_infos, true)
        new_flux_m, new_diff_m = build_symbolic_bucket_func([flux], [dflux], bucket_infos, true; enable_cse=true)

        in_m = Array{Float64}(undef, 1, 2, 2)
        st_m = Array{Float64}(undef, 1, 2, 2)
        in_m[1, :, :] = [1.0 2.0; 3.0 4.0]
        st_m[1, :, :] = [2.0 3.0; 4.0 5.0]

        old_m = old_flux_m(in_m, st_m, pas2)
        new_m = new_flux_m(in_m, st_m, pas2)
        @test _approx_equal(old_m[1], new_m[1])

        in_dm = reshape([1.0, 2.0, 3.0, 4.0], 1, :)
        st_dm = reshape([0.5, 1.0, 1.5, 2.0], 1, :)
        old_dm = old_diff_m(in_dm, st_dm, pas2)
        new_dm = new_diff_m(in_dm, st_dm, pas2)
        @test _approx_equal(old_dm, new_dm)
    end

    @testset "Old vs New Timing By Dimension" begin
        @variables ta pa sa ka

        flux_infos_cmp = HydroInfos(
            inputs=[:ta, :pa],
            outputs=[:qa],
            params=[:ka],
            states=Symbol[],
            nns=Symbol[]
        )

        expr_cmp = (ta + pa)^2 + ka * (ta + pa) + sin(ta)
        pas_cmp = (params=(ka=0.2,),)

        # Scalar (Dim0)
        old_scalar = build_flux_func([expr_cmp], flux_infos_cmp)
        new_scalar = build_symbolic_flux_func([expr_cmp], flux_infos_cmp; dims=0, enable_cse=true)
        scalar_inputs = [1.0, 2.0]

        t_old_scalar = _bench_ns_per_call(old_scalar, scalar_inputs, pas_cmp; warmup=20, iters=3000)
        t_new_scalar = _bench_ns_per_call(new_scalar, scalar_inputs, pas_cmp; warmup=20, iters=3000)

        # Vector (Dim1 via bucket multiply=false)
        infos_v = HydroInfos(
            inputs=[:ta],
            outputs=[:qa],
            states=[:sa],
            params=[:ka],
            nns=Symbol[]
        )
        flux_v = DummyHydroFlux(
            :fv,
            HydroInfos(inputs=[:ta], outputs=[:qa], states=[:sa], params=[:ka], nns=Symbol[]),
            Num[(ta + sa)^2 + ka * (ta + sa) + sin(ta)]
        )
        dflux_v = DummyStateFlux(
            :dfv,
            HydroInfos(inputs=[:ta], outputs=Symbol[], states=[:sa], params=[:ka], nns=Symbol[]),
            Num[ta - sa]
        )

        old_vec, _ = build_bucket_func([flux_v], [dflux_v], infos_v, false)
        new_vec, _ = build_symbolic_bucket_func([flux_v], [dflux_v], infos_v, false; enable_cse=true)
        in_vec = reshape(collect(1.0:1.0:64.0), 1, :)
        st_vec = reshape(collect(2.0:1.0:65.0), 1, :)

        t_old_vec = _bench_ns_per_call(old_vec, in_vec, st_vec, pas_cmp; warmup=20, iters=1200)
        t_new_vec = _bench_ns_per_call(new_vec, in_vec, st_vec, pas_cmp; warmup=20, iters=1200)

        # Matrix (Dim2 via bucket multiply=true)
        old_mat, _ = build_bucket_func([flux_v], [dflux_v], infos_v, true)
        new_mat, _ = build_symbolic_bucket_func([flux_v], [dflux_v], infos_v, true; enable_cse=true)
        in_mat = Array{Float64}(undef, 1, 32, 32)
        st_mat = Array{Float64}(undef, 1, 32, 32)
        in_mat[1, :, :] .= reshape(collect(1.0:1024.0), 32, 32)
        st_mat[1, :, :] .= reshape(collect(2.0:1025.0), 32, 32)

        t_old_mat = _bench_ns_per_call(old_mat, in_mat, st_mat, pas_cmp; warmup=10, iters=200)
        t_new_mat = _bench_ns_per_call(new_mat, in_mat, st_mat, pas_cmp; warmup=10, iters=200)

        @test all(x -> x > 0, (t_old_scalar, t_new_scalar, t_old_vec, t_new_vec, t_old_mat, t_new_mat))

        println("Timing(ns/call) old vs new: Dim0=$(round(t_old_scalar, digits=2)) vs $(round(t_new_scalar, digits=2)); " *
                "Dim1=$(round(t_old_vec, digits=2)) vs $(round(t_new_vec, digits=2)); " *
                "Dim2=$(round(t_old_mat, digits=2)) vs $(round(t_new_mat, digits=2))")
    end
end
