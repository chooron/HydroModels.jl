# 构建ExpHydro模型

## ExpHydro模型介绍

本节将基于 HydroModels.jl 构建一个结构较为简单的水文模型-ExpHydro 模型,以此来开启概念式模型搭建的教程。
首先对 ExpHydro 模型进行介绍,该模型由 Snowpack Bucket 和 Soilwater Bucket 两个计算模块组成,其计算公式如下：

```math
\begin{aligned}
& \text{Snowpack Bucket:} \\
& pet = 29.8 \cdot lday \cdot 24 \cdot 0.611 \cdot \frac{\exp(17.3 \cdot temp)}{temp + 237.3} \cdot \frac{1}{temp + 273.2} && (1) \\
& snowfall = H(T_{min} - temp) \cdot prcp && (2) \\
& rainfall = H(temp - T_{min}) \cdot prcp && (3) \\
& melt = H(temp - T_{max}) \cdot H(snowpack) \cdot \min(snowpack, D_f \cdot (temp - T_{max})) && (4) \\
& \frac{d(snowpack)}{dt} = snowfall - melt && (5) \\
\\
& \text{Soilwater Bucket:} \\
& evap = H(soilwater) \cdot pet \cdot \min(1.0, \frac{soilwater}{S_{max}}) && (6) \\
& baseflow = H(soilwater) \cdot Q_{max} \cdot \exp(-f \cdot \max(0.0, S_{max} - soilwater)) && (7) \\
& surfaceflow = \max(0.0, soilwater - S_{max}) && (8) \\
& flow = baseflow + surfaceflow && (9) \\
& \frac{d(soilwater)}{dt} = rainfall + melt - evap - flow && (10)
\end{aligned}
```

其中:
- $H(x)$ 表示Heaviside阶跃函数，当 $x > 0$ 时为1，否则为0
- $T_{min}, T_{max}, D_f, S_{max}, Q_{max}, f$ 为模型参数
- $temp, lday, prcp$ 为输入变量
- $snowpack, soilwater$ 为状态变量
- 其余变量为中间计算变量

## 完整的模型构建过程

```julia
# import packages
using HydroModels

# define variables and parameters
@variables temp lday pet prcp 
@variables snowfall rainfall melt evap baseflow surfaceflow flow
@variables snowpack soilwater
@parameters Tmin Tmax Df Smax Qmax f

# define step function
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

# define snowpack bucket
fluxes_1 = [
    HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
]
dfluxes_1 = [StateFlux([snowfall] => [melt], snowpack),]
snowpack_bucket = HydroBucket(name=:surface, fluxes=fluxes_1, dfluxes=dfluxes_1)

# define soilwater bucket
fluxes_2 = [
    HydroFlux([soilwater, pet] => [evap], [Smax], exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
    HydroFlux([soilwater] => [baseflow], [Smax, Qmax, f], exprs=[step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))]),
    HydroFlux([soilwater] => [surfaceflow], [Smax], exprs=[max(0.0, soilwater - Smax)]),
    HydroFlux([baseflow, surfaceflow] => [flow], exprs=[baseflow + surfaceflow]),
]
dfluxes_2 = [StateFlux([rainfall, melt] => [evap, flow], soilwater)]
soilwater_bucket = HydroBucket(name=:soil, fluxes=fluxes_2, dfluxes=dfluxes_2)

# define the Exp-Hydro model
exphydro_model = HydroModel(name=:exphydro, components=[snowpack_bucket, soilwater_bucket])
```

## 分段剖析

接下来我们将分块对模型构建过程进行介绍.

首先需要引入HydroModels.jl的依赖:

```julia
using HydroModels
```

通过导入HydroModels模块,可以直接访问HydroFlux,StateFlux,HydroBucket,HydroModel等模块定义的类型,同时还该包导出了ModelingToolkit.jl包中@variables,@parameters等宏,可以用于水文模型中变量和参数的定义和操作, 如下所示:

```julia
@variables temp lday prcp 
@variables snowfall rainfall melt evap baseflow surfaceflow flow pet
@variables snowpack soilwater
@parameters Tmin Tmax Df Smax Qmax f
```

在这段代码中,需要将ExpHydro模型中所有变量进行定义,包括输入变量温度(`temp`),日照时长(`lday`)和降水(`prcp`),同时给出了ExpHydro模型中所有变量的中间计算变量,如`pet`, `snowfall`, `rainfall`, `melt`, `evap`, `baseflow`, `surfaceflow`, `flow`等,以及ExpHydro模型的状态变量`snowpack`和`soilwater`,`以及ExpHydro模型的参数Tmin`, `Tmax`, `Df`, `Smax`, `Qmax`和`f`.

`HydroFlux`的定义需要根据计算公式确定模型的输入输出变量和模型参数,例如在模型的雨雪划分计算公式,该公式的输入变量为`prcp`和`temp`,输出变量为`snowfall`和`rainfall`,模型参数为`Tmin`, 公式转译为`HydroFlux`的结果如下所示:

```julia
# 定义平滑函数
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5
# 定义雨雪划分函数
split_flux = HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp])
```

`StateFlux`的定义是根据状态变量的平衡方程, 一般来说平衡方程是状态变量的变化率等于输入输出通量的差值得到,例如`snowpack`的平衡方程是`snowfall`和`melt`的差值, 公式转译为`StateFlux`的结果如下所示:

```julia
snowpack_dflux = StateFlux([snowfall] => [melt], snowpack)
```

接着根据模型计算公式通过`HydroFlux`和`StateFlux`进行定义, 并将其整合为`HydroBucket`,得到Snowpack Bucket和Soilwater Bucket如下所示:

```julia
# define snowpack bucket
fluxes_1 = [
    HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
]
dfluxes_1 = [StateFlux([snowfall] => [melt], snowpack),]
snowpack_bucket = HydroBucket(name=:surface, fluxes=fluxes_1, dfluxes=dfluxes_1)

# define soilwater bucket
fluxes_2 = [
    HydroFlux([soilwater, pet] => [evap], [Smax], exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
    HydroFlux([soilwater] => [baseflow], [Smax, Qmax, f], exprs=[step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))]),
    HydroFlux([soilwater] => [surfaceflow], [Smax], exprs=[max(0.0, soilwater - Smax)]),
    HydroFlux([baseflow, surfaceflow] => [flow], exprs=[baseflow + surfaceflow]),
]
dfluxes_2 = [StateFlux([rainfall, melt] => [evap, flow], soilwater)]
soilwater_bucket = HydroBucket(name=:soil, fluxes=fluxes_2, dfluxes=dfluxes_2)
```

`HydroBucket`的构建由`HydroFlux`和`StateFlux`组成, 其中`HydroFlux`用于定义模型的计算公式,`StateFlux`用于定义状态变量的平衡方程.

最后将`HydroBucket`组合成`HydroModel`,作为components输入到模型中,得到ExpHydro模型如下所示:

```julia
# define the Exp-Hydro model
exphydro_model = HydroModel(name=:exphydro, components=[snowpack_bucket, soilwater_bucket])
```