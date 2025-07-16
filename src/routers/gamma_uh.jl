using Integrals
using Distributions
using SpecialFunctions
using DSP
using Zygote
using SparseArrays
using BenchmarkTools
β = 0.5
a = 1.0
gamma_func(t, β=3.0, a=2.0) = (1 / t) * (β * t)^a / gamma(a) * exp(-β * t)
uh_vector = gamma_func.(1:1:100, Ref(β), Ref(a))
uh_vector = uh_vector ./ sum(uh_vector)  # Normalize the UH vector
uh_vector_filter = filter(x -> x > maximum(uh_vector) * 0.001, uh_vector)
@btime runoff_vector = conv(uh_vector_filter, rand(1000))

function sparse_compute(input_vec, uh_weight)
    #* the weight of the unit hydrograph is normalized by the sum of the weights
    uh_result = [-(i - 1) => uh_wi .* input_vec ./ sum(uh_weight) for (i, uh_wi) in enumerate(uh_weight)]
    #* sum the matrix
    sum(spdiagm(uh_result...), dims=2)[1:end-length(uh_weight)+1]
end
@btime sparse_compute(rand(1000), uh_vector_filter)