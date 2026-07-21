# HydroModelsSimpleChainsExt

This optional extension activates after `using SimpleChains`. It integrates
the `SimpleChain`/`TurboDense` feed-forward subset with HydroModels' neural
components and factories.

SimpleChains is a CPU forward-only backend in HydroModels. Direct Mooncake
end-to-end differentiation is unsupported: the tested SimpleChains path uses
LLVM intrinsic kernels that Mooncake cannot translate.
