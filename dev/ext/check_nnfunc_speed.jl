using Lux
using ComponentArrays
using Symbolics
using SymbolicUtils
using SymbolicUtils.Code
using StableRNGs
using ModelingToolkit
using BenchmarkTools
using Zygote

q_nn = Lux.Chain(
    Lux.Dense(2 => 16, tanh),
    Lux.Dense(16 => 16, leakyrelu),
    Lux.Dense(16 => 16, leakyrelu),
    Lux.Dense(16 => 1, leakyrelu),
    name=:qnn
)
q_nn_ps, q_nn_st = Lux.setup(StableRNGs.LehmerRNG(1234), q_nn)
q_NN_stateful = Lux.StatefulLuxLayer{true}(q_nn, nothing, q_nn_st) # 不能用symbolic表示
q_nn_params_ca = ComponentVector(first(Lux.setup(StableRNGs.LehmerRNG(1234), q_nn)))
q_nn_states_ca = ComponentVector((Lux.setup(StableRNGs.LehmerRNG(1234), q_nn))[2])
q_nn_params_vec = Vector(q_nn_params_ca)
q_nn_states_vec = Vector(q_nn_states_ca)
q_axes = getaxes(q_nn_params_ca)
q_nn_func1 = (x, q) -> LuxCore.stateless_apply(q_nn, x, ComponentVector(q, q_axes))
q_nn_func2 = (x, q) -> LuxCore.stateless_apply(q_nn, x, q)
q_nn_func3 = (x, q) -> LuxCore.apply(q_NN_stateful, x, q)

loss1(q) = sum(q_nn_func1([1,2], q))
loss2(q) = sum(q_nn_func2([1,2], q))
@btime Zygote.gradient(loss1, q_nn_params_vec)
@btime Zygote.gradient(loss2, q_nn_params_ca)