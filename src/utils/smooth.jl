step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5
export step_func

smoothlogistic_func(S, Smax, r=0.01, e=5.0) = 1 / (1 + exp((S - r * e * Smax) / (r * Smax)))
export smoothlogistic_func