using Lux
using ComponentArrays
using Symbolics
using Symbolics: unwrap
using StableRNGs
using ModelingToolkit
using SymbolicUtils.Code

@variables a b
@variables inputs[1:3] outputs[1:2]

model = Lux.Chain(
    Lux.Dense(3, 10, identity),
    Lux.Dense(10, 2, identity),
)
#* Initialize parameters
init_params = ComponentVector(Lux.initialparameters(StableRNG(42), model))
init_states = Lux.initialstates(StableRNG(42), model)
params_axes = getaxes(init_params)

#* Define parameter variables 
chain_params = first(@parameters ps[1:length(init_params)] [guess = Vector(init_params), description = "Neural network parameters"])
lazy_params = Symbolics.array_term((x, axes) -> ComponentVector(x, axes), chain_params, params_axes, size=size(chain_params))

chain_name = :nn
#* Define neural network input/output variables
nn_input_name, nn_output_name = Symbol(chain_name, :_input), Symbol(chain_name, :_output)
nn_input = first(@variables $(nn_input_name)[1:length(inputs)] [description = "$chain_name Neural network input"])
nn_output = first(@variables $(nn_output_name)[1:length(outputs)] [description = "$chain_name Neural network output"])

flux_expr = unwrap(LuxCore.stateless_apply(model, nn_input, lazy_params)) |> toexpr
nn_input = ones(3,10)
ps = Vector(init_params)
eval(flux_expr)
