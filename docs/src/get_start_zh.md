# Getting Started With HydroModels.jl

`HydroModels.jl` 是一个基于 Julia 语言的现代水文模型框架.它在 SUPERFLEX 设计理念的基础上进行了扩展和增强,具有灵活的模型构建能力和高效的计算性能,并支持深度学习模型集成.本教程将介绍如何使用 `HydroModels.jl` 构建和运行水文模型.

## Installation,.

```julia
] add HydroModels
```

## Build a ExpHydro Model

### Introduction to ExpHydro Model

ExpHydro 是一个由积雪模块和土壤模块组成的简单水文模型.其数学表达式如下：

```math
\begin{aligned}
& \text{Snowpack Bucket:} \\
& pet = 29.8 \cdot lday \cdot 24 \cdot 0.611 \cdot \frac{\exp(17.3 \cdot temp)}{temp + 237.3} \cdot \frac{1}{temp + 273.2} && (1) \\
& :nowfall = H(T_{min} - temp) \cdot prcp && (2) \\
& rainfall = H(temp - T_{m,n}) \cdot prcp,&& (3) \\
& melt = H(temp - T_{max}) \cdot H(snowpack) \cdot \min(snowpack, D_f \cdot (temp - T_{max})) && (4) \\
& \frac{d(snowpack)}{dt} = snowfall - melt && (5) \\
\\
& \text{Soilwater Bucket:} \\,.
& evap = H(soilwater) \cdot pet \cdot \min(1.0, \frac{soilwater}{S_{max}}) && (6) \\
& baseflow = H(soilwater) \cdot Q_{max} \cdot \exp(-f \cdot \max(0.0, S_{max} - soilwater)) && (7) \\
& surfaceflow = \max(0.0, soilwater - S_{max}) && (8) \\
& flow = baseflow + ,urfaceflow && (9) \\.
& \frac{d(soilwater)}{dt} = rainfall + melt - evap - flow && (10)
\end{aligned},,.
```
,,,.
其中：
- $H(x)$ 表示 Heaviside 阶跃函数,当 $x > 0$ 时为 1,否则为 0,.
- $T_{min}, T_{max}, D_f, S_{max}, Q_{max}, f$ 为模型参数
- $temp, lday, prcp$ 为输入变量
- $snowpack, soilwater$ 为状态变量
- 其他变量为中间计算变量

### Build a ExpHydro Model in HydroModels.jl

了解 `ExpHydro` 模型的基本原理后, 我们可以使用 `HydroModels.jl` 来构建这个模型.

#### Import Packages

首先导入 HydroModels.jl 包.

```julia
using HydroModels
```

#### Define Variables and Parameters

定义模型所需参数、状态变量和其他变量(包括降雨、蒸发、融雪、地表流等).

```julia
@variables temp lday pet prcp 
@variables snowfall ra(nf)ll melt evap baseflow surfaceflow flow
@variables snowpack soilwater
@parameters Tmin Tmax Df Smax Qmax f
```

#### Build Flux Formulas by HydroFlux

接下来,我们需要构建模型的各个计算公式,包括状态方程和中间变量计算公式.以融雪公式为例,介绍 `HydroFlux` 的使用方法：

```julia
melt_flux = HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))])
```

在构建 `HydroFlux` 时,需要注意以下几点：

1. 确保表达式中使用的所有变量和参数都已在输入中声明
2. 输入变量和输出变量通过 `Pair` 类型（=>）连接
3. 参数以向量形式传入
4. 表达式以向量形式传入 `exprs` 参数

`HydroFlux` 首先接受输入变量和输出变量所构成的 `Pair` 类型,然后再将参数构成的 `Vector` 传递给 `HydroFlux`,最后根据融雪公式的表达式传递给 `exprs`,完成融雪公式的构建.

需要注意的是,在构建 `HydroFlux` 时,一定要确保 `expr` 中包含的变量和参数均在输入变量和参数中,如果 `exprs` 中存在没有向 `HydroFlux` 提供的输入变量和参数,在后续计算中很有可能会出现变量不存在的错误,这时你就需要检查 `HydroFlux` 的构建是否正确.

当然,构建 HydroFlux 时,参数变量不是必要的,如果不需要参数,可以省略参数变量,例如潜在蒸发量的计算公式：

```julia
pet_flux = HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp(17.3 * temp) / (temp + 237.3) * 1 / (temp + 273.2)])
```

该公式就是以 `temp` 和 `lday` 作为输入变量,并计算了 `pet` 作为输出变量,由于该公式不涉及参数,因此不需要传递参数变量.

值得注意的是,`HydroFlux` 的 `exprs` 接受的类型只能是以表达式所组成的 `Vector`,这是考虑到 `HydroFlux` 有时需要支持多个输出变量,因而需要支持多个表达式,两者的数目必须相同.这里以雨雪划分公式为例.

```julia
HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp])
```

该公式以 `prcp` 和 `temp` 作为输入变量,并计算了 `snowfall` 和 `rainfall` 作为输出变量,同时使用了 `Tmin` 作为参数,由于该公式输出了两个变量 `snowfall` 和 `rainfall`,因此需要传递两个表达式,以此分别表示两个输出变量的计算公式.

#### Build State Formulas by StateFlux

一般而言状态方程是由输入 Flux 之和减去输出 Flux 之和所构建的,因此可以通过 `Pair` 来表示这个输入输出关系,并输入状态变量,以此构建 `StateFlux`,例如雪层状态方程.

```julia
snowpack_flux = StateFlux([snowfall] => [melt], snowpack)
```

该公式以 `snowfall` 作为输入变量,并计算了 `snowpack` 作为输出变量,同时输入了状态变量 `snowpack`,以此构建了雪层状态方程.

一般而言在构建状态方程前,都倾向于将涉及的 Flux 由 `HydroFlux` 全部构建,然后通过 `StateFlux` 结合输入输出通量的 `Pair` 和状态变量得到一个直观的 `StateFlux`.但有时为了避免引入多余的变量,也可以直接构建较为复杂的 `StateFlux`.事实上这个构建方式于 `HydroFlux` 高度一致,同样是将输入输出变量均作为 `StateFlux` 的输入,然后命名状态变量和所需的参数,最后就构建状态方程, 我们重新构建雪层状态方程.

```julia
snowpack_flux = StateFlux([prcp, temp, melt], snowpack, [Tmin], expr=step_func(Tmin - temp) * prcp-melt)
```

在这个状态方程中没有直接计算 `snowfall`,而是通过 `prcp` 和 `temp` 完成计算,结合 `dsnowpack/dt=snowfall-melt` 的状态方程构建一个忽略 `snowfall` 中间变量的复杂状态方程.

#### Build Snowfall Bucket

在完成必要的 `HydroFlux` 和 `StateFlux` 构建后,就可以构建一个 `HydroBucket`,以雪层 `Bucket` 为例,需要将构建的 `HydroFlux` 和 `StateFlux` 传递给 `HydroBucket`.

首先将涉及的 `HydroFlux` 与 `StateFlux` 用 `Vector` 进行存储.

```julia
fluxes_1 = [
    HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
]
dfluxes_1 = [StateFlux([snowfall] => [melt], snowpack)]
```

接着就将这些 `Vector` 变量作为输入构建一个水文计算模块,即 `HydroBucket`.

```julia
snowpack_bucket = HydroBucket(name=:surface, fluxes=fluxes_1, dfluxes=dfluxes_1)
```

通过 `HydroBucket` 就可以完成对一个水文计算模块的构建,在 `HydroBucket` 的构建中,程序会执行一些自动构建方法,以此表达计算模块中的常微分方程和其余变量的计算函数,详情请参考实施细节.

#### Build Exphydro Model

按照相同的逻辑,就可以构建模型中的土壤计算模块,接着将两个模块拼接到一块构建 `Vector` 类型,并输入至 `HydroModel` 中,就完成了 `ExpHydro` 模型的构建,完整构建代码如下

```julia
#* import packages
using HydroModels

#* define variables and parameters
@variables temp lday pet prcp 
@variables snowfall rainfall melt evap baseflow surfaceflow flow
@variables snowpack soilwater
@parameters Tmin Tmax Df Smax Qmax f

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

#* define snowpack bucket
fluxes_1 = [
    HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
]
dfluxes_1 = [StateFlux([snowfall] => [melt], snowpack),]
snowpack_bucket = HydroBucket(name=:surface, fluxes=fluxes_1, dfluxes=dfluxes_1)

#* define soilwater bucket
fluxes_2 = [
    HydroFlux([soilwater, pet] => [evap], [Smax], exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
    HydroFlux([soilwater] => [baseflow], [Smax, Qmax, f], exprs=[step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))]),
    HydroFlux([soilwater] => [surfaceflow], [Smax], exprs=[max(0.0, soilwater - Smax)]),
    HydroFlux([baseflow, surfaceflow] => [flow], exprs=[baseflow + surfaceflow]),
]
dfluxes_2 = [StateFlux([rainfall, melt] => [evap, flow], soilwater)]
soilwater_bucket = HydroBucket(name=:soil, fluxes=fluxes_2, dfluxes=dfluxes_2)

#* define the Exp-Hydro model
exphydro_model = HydroModel(name=:exphydro, components=[snowpack_bucket, soilwater_bucket])
```

这样就完成了一个 `ExpHydro` 模型的构建,打印这个类型,可以看到模型包含的基本信息.

```txt
HydroModel: exphydro
   Components: surface, soil
   Inputs: [temp, lday, prcp]
   States: [snowpack, soilwater]
   Outputs: [pet, snowfall, rainfall, melt, evap, baseflow, surfaceflow, flow]
   Parameters: [Tmin, Tmax, Df, Smax, Qmax, f]
   Components:
       Fluxes: 0 fluxes
       Buckets: 2 buckets
       Routes: 0 routes
```

这些信息中展示了模型的名称,包含的组建名称,输入变量,状态变量,输出变量和参数变量,以及包含的 `Flux`,`Bucket` 和 `Route` 的数量.

#### Embedding Neural Network into ExpHydro Model

`ExpHydro` 模型近年来作为一个基础水文模型被大量应用于神经网络耦合水文模型的研究,例如 `M50`,`M100`,`ENN` 等.这些模型都是将 `ExpHydro` 模型作为基础模型,然后通过神经网络替换模型的部分计算公式（如蒸发和产流）,以此实现提升模型预测精度的目的.

下面,我们将介绍如何使用 `HydroModels.jl` 将神经网络嵌入至 `ExpHydro` 模型中.

#### Define Neural Network by Lux.jl

`HydroModels.jl` 支持通过 `Lux.jl` 来定义神经网络,这里以一个简单的全连接神经网络为例.

```julia
using Lux

# define the ET NN and Q NN
ep_nn = Lux.Chain(
    Lux.Dense(3 => 16, tanh),
    Lux.Dense(16 => 16, leakyrelu),
    Lux.Dense(16 => 1, leakyrelu),
    name=:epnn
)

q_nn = Lux.Chain(
    Lux.Dense(2 => 16, tanh),
    Lux.Dense(16 => 16, leakyrelu),
    Lux.Dense(16 => 1, leakyrelu),
    name=:qnn
)
```

首先根据 `Lux.jl` 构建了两个全连接神经网络模型 `ep_nn` 和 `q_nn`,接着就可以通过 `HydroModel.jl` 提供的 `NeuralFlux` 来构建神经网络公式.

```julia
@variables norm_snw norm_slw norm_temp norm_prcp
@variables log_evap_div_lday log_flow

ep_nn_flux = NeuralFlux([norm_snw, norm_slw, norm_temp] => [log_evap_div_lday], ep_nn)
q_nn_flux = NeuralFlux([norm_slw, norm_prcp] => [log_flow], q_nn)
```

`NeuralFlux` 是 `HydroFlux` 的一个衍生类型,他同样是接受输入输出变量的 `Pair`,但不同的是他不需要提供模型参数和计算公式,他需要将神经网络作为参数传递给 `NeuralFlux`, 以此构建神经网络公式.
值得注意的是,神经网络的输入输出变量维度需要与 `NeuralFlux` 的输入输出变量维度一致.

然后完成 `NeuralFlux` 的构建后, 就可以将其与其余 `HydroFlux` 一同传递给 `HydroBucket`,以此构建神经网络耦合的水文模型.

```julia
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

```

这个土壤模块的公式包括了标准化公式,蒸发公式和产流公式,接着就可以将这些公式传递给 HydroBucket,结合状态方程完成土壤模块的构建.

```julia
state_expr = rainfall + melt - step_func(soilwater) * lday * exp(log_evap_div_lday) - step_func(soilwater) * exp(log_flow)
soil_dfluxes = [StateFlux([soilwater, rainfall, melt, lday, log_evap_div_lday, log_flow], soilwater, expr=state_expr)]
soilwater_bucket = HydroBucket(name=:soil, fluxes=soil_fluxes, dfluxes=soil_dfluxes)
m50_model = HydroModel(name=:m50, components=[snowpack_bucket, soilwater_bucket])
```

## Run a ExpHydro Model

在完成模型的构建后,由于模型是 `callable` 的, 因此可以将输入和参数输入至模型中,并得到输出结果.

### Prepare Input Data

模型输入的观测数据限于 `AbstractMatrix` 类型的数据,维度为输入变量数目乘以时间步数,例如输入变量有 3 个,时间步数为 1000,那么输入的观测数据维度为 (3, 1000).值得注意的是,输入变量的所在维度的索引会对模型计算结果产生影响,因为模型无法自动识别 Matrix 各行对应的输入变量,因此使用者可以先准备好一个 `Dict`,`NamedTuple` 类型,然后根据 `HydroModels.get_input_names(model)` 来获取模型输入变量的名称,将这个名称结合 `Dict` 或 `NamedTuple` 构建 `Matrix` 类型的输入数据.

```julia
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
input_arr = Matrix(reduce(hcat, collect(input[HydroModels.get_input_names(exphydro_model)]))')
```

### Prepare Parameters and Initial States

对于参数的准备, `HydroModels.jl` 是接受 `ComponentVector` 类型的参数 `(using ComponentArrays)`,通过参数名称和参数值存储在 `ComponentVector` 中,例如

```julia
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
```

对于初始状态的准备, HydroModels.jl 同样是接受 `ComponentVector` 类型的初始状态,通过状态名称和状态值存储在 `ComponentVector` 中,例如

```julia
init_states = ComponentVector(snowpack=0.0, soilwater=1303.004248)
```

接着将参数和初始状态存储至同一个 `ComponentVector` 中（有时初始状态同样是作为待优化的参数,因此需要将其与参数一并存储,实现统一管理）.

```julia
pas = ComponentVector(params=params, initstates=init_states)
```

### Run the ExpHydro Model

最后将输入数据和参数输入至模型中,就可以得到输出结果.

```julia
result = exphydro_model(input_arr, pas)
```

模型计算输出同样是一个 `AbtractMatrix` 的类型,维度为输出变量数目（包括状态变量）乘以时间步数,例如输出变量有 8 个,时间步数为 1000,那么输出结果的维度为 (8, 1000).如果你想导出为 `DataFrame` 类型,可以使用 `HydroModels.get_output_names(model)` 来获取输出变量的名称,然后结合 `Dict` 或 `NamedTuple` 构建 `DataFrame` 类型的输出数据.

```julia
states_and_output_names = vcat(HydroModels.get_state_names(exphydro_model), HydroModels.get_output_names(exphydro_model))
output = NamedTuple{Tuple(states_and_output_names)}(eachslice(result, dims=1))
df = DataFrame(output)
```

计算结果如下.

```txt
10000×10 DataFrame
   Row │ snowpack  soilwater  pet      snowfall  rainfall  melt      evap      baseflow   surfaceflow  flow      
       │ Float64   Float64    Float64  Float64   Float64   Float64   Float64   Float64    Float64      Float64   
───────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │    0.0         0.0      3.1        0.0   1305.22  1.13779   0.868735  0.0212194          0.0  0.0212194
     2 │    0.0         0.0      4.24       0.0   1308.28  1.51019   1.15578   0.0223371          0.0  0.0223371
     3 │    0.0         0.0      8.02       0.0   1315.03  1.63204   1.25547   0.0250095          0.0  0.0250095
     4 │    0.0         0.0     15.27       0.0   1329.34  1.21771   0.946937  0.0317802          0.0  0.0317802
     5 │    0.0         0.0      8.48       0.0   1336.99  1.02779   0.803845  0.0361228          0.0  0.0361228
   ⋮   │    ⋮          ⋮         ⋮        ⋮         ⋮         ⋮         ⋮          ⋮           ⋮           ⋮
  9996 │  254.264       0.0      0.0       -0.0   1521.37  0.221762  0.197362  0.791897           0.0  0.791897
  9997 │  266.324      12.06     0.0       -0.0   1520.37  0.238907  0.21248   0.778688           0.0  0.778688
  9998 │  277.874      11.55     0.0       -0.0   1519.33  0.288362  0.256291  0.765307           0.0  0.765307
  9999 │  279.714       1.84     0.0       -0.0   1518.29  0.311022  0.27624   0.752073           0.0  0.752073
 10000 │  279.854       0.14     0.0       -0.0   1517.38  0.176836  0.156967  0.740711           0.0  0.740711
```

### Running with Config

当然模型计算不仅限于这些内容,为支持计算过程中的一些必要设置,`HydroModel` 在计算时能够接受 `config` 参数.

例如在默认情况下,模型求解常微分方程时仅使用了 Euler 方法,但有时需要使用更精确的求解方法,这时就需要使用 `config` 参数来设置求解方法.

为演示这个功能需要额外导入 `OrdinaryDiffEq.jl` 和 `HydroModelTools.jl`,然后设置求解方法.

```julia
using OrdinaryDiffEq
using HydroModelTools
using DataInterpolations

config = (solver=ODESolver(alg=Tsit5(), abstol=1e-3, reltol=1e-3), interp=LinearInterpolation)
output = exphydro_model(input_arr, pas, config=config)
```

`config` 参数的设置方式为 `NamedTuple` 类型,其中 `solver` 参数用于设置求解方法, `interp` 参数用于设置插值方法.`ODESolver` 是 `HydroModelTools.jl` 提供的求解器包装类,对求解结果进行了一些数据转换,它接受 `OrdinaryDiffEq.jl` 的求解方法作为参数,并设置求解方法的参数.

### Run the Neural Network Enhanced Model

在前面我们构建了一个用神经网络耦合的 `ExpHydro` 模型,现在我们就可以使用这个模型来计算结果.

首先需要准备模型的参数,在 `Exphydro` 模型参数基础上我们还需要额外准备神经网络的参数.

```julia
using StableRNGs

ep_nn_params = Vector(ComponentVector(LuxCore.initialparameters(ep_nn)))
q_nn_params = Vector(ComponentVector(LuxCore.initialparameters(q_nn)))
```

基于 `Lux.jl` 的模型能够将参数与模型结构进行解耦,这也是我们选择使用 `Lux.jl` 构建神经网络的首要原因,在使用 `LuxCore.initialparameters` 后,我们将参数（`NamedTuple`）转换为 `Vector` 类型,然后与模型参数和初始状态一同构建模型参数.

```julia
nn_pas = ComponentVector(epnn=ep_nn_params, qnn=q_nn_params)
pas = ComponentVector(params=params, initstates=init_states, nns=nn_pas)
```

值得注意的是,需要使用 `nns` 这个键名存储神经网络的参数,保证与模型参数和初始状态独立,同时神经网络参数 `nn_pas` 需要使用模型名称和参数作为键值对构建参数.

最后将模型参数输入至模型中,就可以得到输出结果.

```julia
output = m50_model(input_arr, pas)
```

## Conclusion

本教程详细介绍了 `HydroModels.jl` 的基本使用方法,包括：
- 模型的构建过程
- 参数和初始状态的设置
- 模型的运行和计算

通过本教程,用户可以快速掌握 `HydroModels.jl` 的核心功能,并开始构建自己的水文模型.

`HydroModels.jl` 还提供了更多高级功能,包括：
- 半分布式和分布式水文模型的构建与计算
- 神经网络耦合的参数自适应模型
- 支持单位线汇流的 `GR4J` 模型
- 模型参数优化
- 基于 `Wrapper` 的模型扩展功能

这些高级特性将在后续教程中详细介绍.
