function r2(simulated_array::AbstractVector{T}, observed_array::AbstractVector{T}; kwargs...)::T where {T}
    a = observed_array .- mean(observed_array)
    b = simulated_array .- mean(simulated_array)
    sum(a .* b)^2 / (sum(a .^ 2) * sum(b .^ 2))
end