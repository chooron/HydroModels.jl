# 比较superflexpy, MARRMoT 和 HydroModels.jl三种库对于自定义模型的编程风格

本内容对应了论文的4.1节的内容,我们在这里希望能够讨论三种软件库对于模型计算公式的自定义的编程风格的差异. 我们以GR4J模型的实习过程为例,对模型计算公式的自定义的编程风格进行比较.

## [Superflexpy](https://github.com/dalmo1991/superflexPy/tree/master)

Superflexpy是基于SUPERFLEX思想下的python代码实现, 针对于水文模型中一些计算功能的实现, 该框架一般而言是直接使用python语法表述模型的计算公式, 如下所示.

```python
@staticmethod
def _flux_function_python(S, S0, ind, P, x1, alpha, beta, ni, PET, dt):
    if ind is None:
        return (
            [
                P * (1 - (S / x1) ** alpha), # Ps
                -PET * (2 * (S / x1) - (S / x1) ** alpha), # Evaporation
                -((x1 ** (1 - beta)) / ((beta - 1) * dt)) * (ni ** (beta - 1)) * (S**beta), # Perc
            ],
            0.0,
            S0 + P * (1 - (S / x1) ** alpha) * dt,
        )
    else:
        return (
            [
                P[ind] * (1 - (S / x1[ind]) ** alpha[ind]), # Ps
                -PET[ind] * (2 * (S / x1[ind]) - (S / x1[ind]) ** alpha[ind]), # Evaporation
                -((x1[ind] ** (1 - beta[ind])) / ((beta[ind] - 1) * dt[ind]))
                * (ni[ind] ** (beta[ind] - 1))
                * (S ** beta[ind]), # Perc
            ],
            0.0,
            S0 + P[ind] * (1 - (S / x1[ind]) ** alpha[ind]) * dt[ind],
            [
                -(P[ind] * alpha[ind] / x1[ind]) * ((S / x1[ind]) ** (alpha[ind] - 1)),
                -(PET[ind] / x1[ind]) * (2 - alpha[ind] * ((S / x1[ind]) ** (alpha[ind] - 1))),
                -beta[ind]
                * ((x1[ind] ** (1 - beta[ind])) / ((beta[ind] - 1) * dt[ind]))
                * (ni[ind] ** (beta[ind] - 1))
                * (S ** (beta[ind] - 1)),
            ],
        )

@staticmethod
@nb.jit(
    "Tuple((UniTuple(f8, 3), f8, f8, UniTuple(f8, 3)))"
    "(optional(f8), f8, i4, f8[:], f8[:], f8[:], f8[:], f8[:], f8[:], f8[:])",
    nopython=True,
)
def _flux_function_numba(S, S0, ind, P, x1, alpha, beta, ni, PET, dt):
    return (
        (
            P[ind] * (1 - (S / x1[ind]) ** alpha[ind]),  # Ps
            -PET[ind] * (2 * (S / x1[ind]) - (S / x1[ind]) ** alpha[ind]),  # Evaporation
            -((x1[ind] ** (1 - beta[ind])) / ((beta[ind] - 1) * dt[ind]))
            * (ni[ind] ** (beta[ind] - 1))
            * (S ** beta[ind]),  # Perc
        ),
        0.0,
        S0 + P[ind] * (1 - (S / x1[ind]) ** alpha[ind]) * dt[ind],
        (
            -(P[ind] * alpha[ind] / x1[ind]) * ((S / x1[ind]) ** (alpha[ind] - 1)),
            -(PET[ind] / x1[ind]) * (2 - alpha[ind] * ((S / x1[ind]) ** (alpha[ind] - 1))),
            -beta[ind]
            * ((x1[ind] ** (1 - beta[ind])) / ((beta[ind] - 1) * dt[ind]))
            * (ni[ind] ** (beta[ind] - 1))
            * (S ** (beta[ind] - 1)),
        ),
    )
```

可以看出, superflexpy框架在实现GR4J的production模块计算时[origin code](https://github.com/dalmo1991/superflexPy/blob/master/superflexpy/implementation/elements/gr4j.py),不仅需要分别构建numpy和numba的计算代码,在函数中需要给出模型中涉及的水文通量`Ps`,`Evaporation`,`Perc`,同时还要给出状态变量的平衡方程,一个代码中需要耦合多个计算过程的实现代码,不仅导致代码的复用的灵活性较差(仅能在模块层面进行复用),模型编写的代码量和重复度也较高.

## [MARRMoT](https://github.com/wknoben/MARRMoT)

MARRMoT是一个由MATLAB语言开发的水文模型框架,该框架根据不同水文模型会分别构建模型类型,储存模型涉及的计算公式和状态变化方程,GR4J模型在该框架下的实现代码(简化了变量命名)如下所示.

```matlab
function [dS, fluxes] = model_fun(obj, S)
    % definitions
    ...
    % fluxes functions
    flux_pn   = max(P-Ep,0);
    flux_en   = max(Ep-P,0);
    flux_ef   = P - flux_pn;
    flux_ps   = saturation_4(S1,x1,flux_pn);
    flux_es   = evap_11(S1,x1,flux_en);
    flux_perc = percolation_3(S1,x1);
    flux_q9   = route(.9.*(flux_pn - flux_ps + flux_perc), uh_q9);
    flux_q1   = route(.1.*(flux_pn - flux_ps + flux_perc), uh_q1);
    flux_fr   = recharge_2(3.5,S2,x3,x2);
    flux_fq   = flux_fr;
    flux_qr   = baseflow_3(S2,x3);
    flux_qt   = flux_qr + max(flux_q1 + flux_fq,0);
    % this flux is not included in original MARRMoT,
    % but it is useful to calculate the water balance
    flux_ex = flux_fr + max(flux_q1 + flux_fq,0) - flux_q1;      
    % stores ODEs
    dS1 = flux_ps - flux_es - flux_perc;
    dS2 = flux_q9 + flux_fr - flux_qr;
    % outputs
    dS = [dS1 dS2];
    fluxes = [flux_pn,   flux_en, flux_ef, flux_ps, flux_es,...
                flux_perc, flux_q9, flux_q1, flux_fr, flux_fq,...
                flux_qr,   flux_qt, flux_ex];
end
```

可以看出,MARRMoT框架将每一个Flux计算公式由函数进行表达,因此使用者可以直接调取函数进行计算(可以根据提供的文档进行调用),这种方式相比superflexpy的确简化了很多代码编写需求. 然而MARRMoT对于模块的重复率支持性较差,不能够像superflexpy那样进行复用. 同时将变量赋值,通量计算和状态变化方程以一个函数包括,计算内容的解耦性较差.

## HydroModels.jl

```julia
#* define the production store
prod_funcs = [
    HydroFlux([prcp, ep] => [pn, en], exprs=[prcp - min(prcp, ep), ep - min(prcp, ep)]),
    HydroFlux([pn, soilwater] => [ps], [x1], exprs=[max(0.0, pn * (1 - (soilwater / x1)^2))]),
    HydroFlux([en, soilwater] => [es], [x1], exprs=[en * (2 * soilwater / x1 - (soilwater / x1)^2)]),
    HydroFlux([soilwater] => [perc], [x1], exprs=[((x1)^(-4)) / 4 * ((4 / 9)^(4)) * (soilwater^5)]),
    HydroFlux([pn, ps, perc] => [pr], exprs=[pn - ps + perc]),
    HydroFlux([pr] => [slowflow, fastflow], exprs=[0.9 * pr, 0.1 * pr]),
]
prod_dfuncs = [StateFlux([ps] => [es, perc], soilwater)]
prod_ele = HydroBucket(name=:gr4j_prod, funcs=prod_funcs, dfuncs=prod_dfuncs)

#* uh function
uh_flux_1 = HydroModels.UnitHydrograph([slowflow] => [slowflow_routed], x4, uhfunc=HydroModels.UHFunction(:UH_1_HALF), solvetype=:SPARSE)
uh_flux_2 = HydroModels.UnitHydrograph([fastflow] => [fastflow_routed], x4, uhfunc=HydroModels.UHFunction(:UH_2_FULL), solvetype=:SPARSE)

#* define the routing store
rst_funcs = [
    HydroFlux([routingstore] => [exch], [x2, x3], exprs=[x2 * abs(routingstore / x3)^3.5]),
    HydroFlux([routingstore, slowflow, exch] => [routedflow], [x3], exprs=[x3^(-4) / 4 * (routingstore + slowflow + exch)^5]),
    HydroFlux([routedflow, fastflow_routed, exch] => [flow], exprs=[routedflow + max(fastflow_routed + exch, 0.0)]),
]
rst_dfuncs = [StateFlux([exch, slowflow] => [routedflow], routingstore)]
rst_ele = HydroBucket(name=:gr4j_rst, funcs=rst_funcs, dfuncs=rst_dfuncs)

gr4j_model = HydroModel(name=:gr4j, components=[prod_ele, uh_flux_1, uh_flux_2, rst_ele])
```

回到HydroModels.jl框架中对于GR4J模型的实现, 在该代码中我们依次对production store, unit hydrograph, routing store进行了定义, 在这个过程中,在变量对应的前提下,模型能够实现可复用的能力. 在这个搭建过程中, 包括HydroFlux, StateFlux, UnitHydrograph, HydroModel, HydroBucket等类的实例化,这些类型的实例化方式不是像前两者一样使用函数的方式构建的, 而是通过符号编程的方式对计算公式进行表达, 在模型构建完成后是通过输入数据和参数进行模型计算,可以看出这种设计思路是将模型结构与公式同参数和数据完全解耦,有着极高的可复用性.