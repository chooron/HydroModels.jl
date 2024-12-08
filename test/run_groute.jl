
include("../src/HydroModels.jl")
using ModelingToolkit
using ComponentArrays

@variables q1 q1_routed s_river
@parameters lag
flwdir = [1 4 8; 1 4 4; 1 1 2]
positions = [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]

ndtypes = [:ntype1, :ntype2, :ntype3]
rflux = HydroModels.HydroFlux([q1, s_river] => [q1_routed], [lag], exprs=[s_river / (1 + lag) + q1])
params = ComponentVector(NamedTuple{Tuple(ndtypes)}([(lag=0.2,) for _ in eachindex(ndtypes)]))
initstates = ComponentVector(NamedTuple{Tuple(ndtypes)}([(s_river=0.1,) for _ in eachindex(ndtypes)]))
pas = ComponentVector(; params, initstates)
nodeids = [:nid1, :nid2, :nid3, :nid4, :nid5, :nid6, :nid7, :nid8, :nid9]
route = HydroModels.GridRoute(rfunc=rflux, rstate=s_river, flwdir=flwdir, positions=positions, nodeids=nodeids)

input_arr = ones(1, 9, 20)
timeidx = collect(1:20)
node_types = [:ntype1, :ntype2, :ntype3, :ntype2, :ntype1, :ntype2, :ntype3, :ntype1, :ntype3]
config = (solver=HydroModels.ODESolver(saveat=timeidx), interp=LinearInterpolation, ptypes=node_types, stypes=node_types, timeidx=timeidx)
@btime output_arr = route(input_arr, pas, config=config);