# Embedding Neural Networks in ExpHydro Model: Implementation of [M50](https://hess.copernicus.org/articles/26/5085/2022/) Model

[marv-in](https://github.com/marv-in)在论文中提出了M50和M100两个模型,代码储存于[HydroNODE](https://github.com/marv-in/HydroNODE)中,开创了神经网络嵌入水文模型的先河(重点在于实现神经网络替换了ExpHydro模型的计算公式,并参与到平衡方程的求解中),这个模型也是对该仓库的启蒙,我们发现当前模型编写仍存在一定难度,因此我们创建这个仓库的首要目的就是简化神经网络嵌入水文模型这个过程,下面我们将对比源代码与本仓库在M50和M100模型上的实现差异,为使用者展现出一种新的模型搭建方式.

## HydroNODE中M50的实现方法

首先我们展示了HydroNODE中M50和M100模型的实现:

```julia
# define the ODE problem
function NeuralODE_M50_core!(dS,S,p,t)

    Tmin, Tmax, Df = (p_bucket_precal...,)

    Lday = itp_Lday(t)
    P    = itp_P(t)
    T    = itp_T(t)

    g_ET = ann_ET([norm_S0(S[1]), norm_S1(S[2]), norm_T(T)],p[:p1]) #p[idcs_params_ann_ET])
    g_Q = ann_Q([norm_S1(S[2]), norm_P(P)],p[:p2]) #p[idcs_params_ann_Q])

    melting = M(S[1], T, Df, Tmax)
    dS[1] = Ps(P, T, Tmin) - melting
    dS[2] = Pr(P, T, Tmin) + melting - step_fct(S[2])*Lday*exp(g_ET[1])- step_fct(S[2])*exp(g_Q[1])

end
# build and solve ODE problem
prob = ODEProblem(NeuralODE_M50_core!, S_init, Float64.((t_out[1], maximum(t_out))), p)
sol = solve(prob, BS3(), dt=1.0, saveat=t_out, reltol=1e-3, abstol=1e-3, sensealg=BacksolveAdjoint(autojacvec=ZygoteVJP()))
# calculate Qout
P_interp = norm_P.(itp_P.(t_out))
S1_ = norm_S1.(sol[2,:])
Qout_ =  exp.(ann_Q(permutedims([S1_ P_interp]),p[:p2])[1,:])
```

在HydroNODE的代码中可以看出,模型构建分为以下三个步骤:

- ODE函数的定义
- ODE问题的构建与求解
- 以及根据求解的中间状态计算Qout

这个代码的构建风格遵循了最标准的DifferentialEquations.jl库的风格,然而这个代码中需要手动搭建模型的计算公式(`M`,`Ps`,`Pr`)并与状态方程表达式组合在一起(`dS`),模型的耦合性极高,同时,为获取`Qout`等中间计算通量,还需要额外调用对应的计算公式,代码编写的繁琐性相对较高.

## HydroModels.jl中M50的实现方法

首先需要定义M50模型中设计的参数和变量:

```julia
using HydroModels

#! parameters in the Exp-Hydro model
@parameters Tmin Tmax Df Smax f Qmax
#! parameters in normalize flux
@parameters snowpack_std snowpack_mean
@parameters soilwater_std soilwater_mean
@parameters prcp_std prcp_mean
@parameters temp_std temp_mean

#! hydrological flux in the Exp-Hydro model
@variables prcp temp lday pet rainfall snowfall
@variables snowpack soilwater lday pet
@variables melt log_evap_div_lday log_flow flow
@variables norm_snw norm_slw norm_temp norm_prcp

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5
```


相比较普通的概念式水文模型,M50模型使用神经网络模型替换了水文通量的计算公式(`ET`和`Q`),这个神经网络相比普通的计算公式有着明显的差异,因此在HydroModels.jl中采用了`NeuralFlux`来表示神经网络, 用于ET和Q预测的NeuralFlux如下表示:

```julia
# define the ET NN and Q NN
ep_nn = Lux.Chain(
    Lux.Dense(3 => 16, tanh),
    Lux.Dense(16 => 16, leakyrelu),
    Lux.Dense(16 => 1, leakyrelu),
    name=:epnn
)
ep_nn_params = Vector(ComponentVector(first(Lux.setup(StableRNGs.LehmerRNG(1234), ep_nn))))
q_nn = Lux.Chain(
    Lux.Dense(2 => 16, tanh),
    Lux.Dense(16 => 16, leakyrelu),
    Lux.Dense(16 => 1, leakyrelu),
    name=:qnn
)
q_nn_params = Vector(ComponentVector(first(Lux.setup(StableRNGs.LehmerRNG(1234), q_nn))))


ep_nn_flux = NeuralFlux([norm_snw, norm_slw, norm_temp] => [log_evap_div_lday], ep_nn)
q_nn_flux = NeuralFlux([norm_slw, norm_prcp] => [log_flow], q_nn)
```

构建M50模型相对特殊的一点在于需要针对嵌入的神经网络使用NeuralFlux模型进行包装,为神经网络赋予输入输出的信息,并在构建过程中将神经网络转换为表达式.

然后我们可以构建M50模型的其余部分:


```julia
# define the snow pack reservoir
snow_fluxes = [
    HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * min(snowpack, Df * (temp - Tmax))]),
]
snow_dfluxes = [StateFlux([snowfall] => [melt], snowpack)]
snow_ele = HydroBucket(name=:exphydro_snow, fluxes=snow_fluxes, dfluxes=snow_dfluxes)

# define the soil water reservoir
soil_fluxes = [
    #* normalize
    HydroFlux([snowpack, soilwater, prcp, temp] => [norm_snw, norm_slw, norm_prcp, norm_temp],
        [snowpack_mean, soilwater_mean, prcp_mean, temp_mean, snowpack_std, soilwater_std, prcp_std, temp_std],
        exprs=[(var - mean) / std for (var, mean, std) in zip([snowpack, soilwater, prcp, temp],
            [snowpack_mean, soilwater_mean, prcp_mean, temp_mean],
            [snowpack_std, soilwater_std, prcp_std, temp_std]
        )]),
    ep_nn_flux,
    q_nn_flux
]
state_expr = rainfall + melt - step_func(soilwater) * lday * exp(log_evap_div_lday) - step_func(soilwater) * exp(log_flow)
soil_dfluxes = [StateFlux([soilwater, rainfall, melt, lday, log_evap_div_lday, log_flow], soilwater, expr=state_expr)]
soil_ele = HydroBucket(name=:m50_soil, fluxes=soil_fluxes, dfluxes=soil_dfluxes)
convert_flux = HydroFlux([log_flow] => [flow], exprs=[exp(log_flow)])
# define the Exp-Hydro model
m50_model = HydroModel(name=:m50, components=[snow_ele, soil_ele, convert_flux]);
```

上述代码中我们对M50模型进行了一个完整的表现,首先包括ExpHydro中Snowpack Bucket的实现,其次就是对于神经网络嵌入后Soilwater Bucket的实现. 在Soilwater Bucket模块实现中,可以通过定义的NeuralFlux(`ep_nn_flux`和`q_nn_flux`)来表达嵌入的神经网络,并与其他`HydroFlux`进行组合,结合`StateFlux`构建Soilwater Bucket.
