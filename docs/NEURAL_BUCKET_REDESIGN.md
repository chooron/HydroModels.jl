# NeuralBucket 重新设计说明

## 设计理念

### ODE 与 RNN 的结构相似性

水文模型的 ODE 结构本质上与 RNN 高度相似：

**传统 HydroBucket (ODE):**
```
dS/dt = f(S, input, params)
S_t = S_{t-1} + dt * f(S_{t-1}, input_t, params)
output_t = g(S_t, input_t, params)
```

**NeuralBucket (RNN-like):**
```
flux_t = flux_network(S_{t-1}, input_t)
S_t = state_network(S_{t-1}, flux_t)
output_t = output_network(flux_t)
```

两者都涉及：
1. 跨时间步维护状态
2. 从状态和输入计算输出
3. 基于计算结果更新状态

---

## 核心结构

### NeuralBucket 组件

```julia
struct NeuralBucket <: Lux.AbstractLuxLayer
    name::Symbol
    flux_network::AbstractLuxLayer      # 计算通量
    state_network::AbstractLuxLayer     # 更新状态  
    output_network::AbstractLuxLayer    # 计算输出
    n_inputs::Int
    n_states::Int
    n_outputs::Int
    hru_types::Vector{Int}
    infos::HydroInfos
end
```

### 三个神经网络的作用

1. **Flux Network** - 通量计算网络
   - 输入: `[states; inputs]` (n_states + n_inputs)
   - 输出: `fluxes` (n_fluxes)
   - 作用: 类似于 HydroBucket 中的 `@hydroflux`

2. **State Network** - 状态更新网络
   - 输入: `[states; fluxes]` (n_states + n_fluxes)
   - 输出: `new_states` (n_states)
   - 作用: 类似于 HydroBucket 中的 `@stateflux`（状态微分）

3. **Output Network** - 输出计算网络
   - 输入: `fluxes` (n_fluxes)
   - 输出: `outputs` (n_outputs)
   - 作用: 从通量中提取最终输出

---

## 与 HydroBucket 的对比

| 特性 | HydroBucket | NeuralBucket |
|------|------------|--------------|
| **建模方式** | 符号表达式 (ODE) | 神经网络 (RNN-like) |
| **通量计算** | `@hydroflux` 宏 | `flux_network` |
| **状态更新** | `@stateflux` 宏 | `state_network` |
| **求解器** | 需要 (Mutable/Immutable/ODE) | 不需要 (离散步进) |
| **状态意义** | 物理状态 | 可学习的表示 |
| **参数** | 物理参数 | 神经网络权重 |
| **可解释性** | ✅ 高 | ❌ 低 |
| **灵活性** | ⚠️ 受物理约束 | ✅ 高度灵活 |
| **Lux 集成** | ❌ 否 | ✅ 是 |

---

## 使用方法

### 方法 1: 使用便捷函数

```julia
using HydroModels, Lux, ComponentArrays, Random

# 创建 NeuralBucket
bucket = create_neural_bucket(
    name = :rainfall_runoff,
    n_inputs = 2,           # 输入数量
    n_states = 2,           # 状态数量
    n_outputs = 1,          # 输出数量
    n_fluxes = 3,           # 中间通量数量
    hidden_size = 16,       # 隐藏层大小
    inputs = [:prcp, :temp],
    states = [:snowpack, :soilwater],
    outputs = [:runoff]
)

# 初始化参数
rng = Random.default_rng()
ps_lux = LuxCore.initialparameters(rng, bucket)

# 转换为 HydroModels 格式
params = ComponentVector(
    nns = (
        rainfall_runoff_flux = ps_lux.flux,
        rainfall_runoff_state = ps_lux.state,
        rainfall_runoff_output = ps_lux.output
    )
)

# 准备输入数据
input = rand(Float32, 2, 100)  # (n_inputs × n_timesteps)

# 运行
config = HydroConfig()
output = bucket(input, params, config)  # (n_outputs × n_timesteps)
```

### 方法 2: 自定义网络架构

```julia
# 定义自定义网络
flux_net = Chain(
    Dense(4 => 32, relu),    # n_states + n_inputs
    Dense(32 => 16, tanh),
    Dense(16 => 3),          # n_fluxes
    name = :flux_net
)

state_net = Chain(
    Dense(5 => 16, relu),    # n_states + n_fluxes
    Dense(16 => 2),          # n_states
    name = :state_net
)

output_net = Chain(
    Dense(3 => 1),           # n_fluxes → n_outputs
    name = :output_net
)

# 创建 NeuralBucket
bucket = NeuralBucket(;
    name = :custom,
    flux_network = flux_net,
    state_network = state_net,
    output_network = output_net,
    n_inputs = 2,
    n_states = 2,
    n_outputs = 1,
    inputs = [:prcp, :temp],
    states = [:snowpack, :soilwater],
    outputs = [:runoff]
)
```

---

## 单时间步计算 (Lux 接口)

```julia
# 初始化 Lux 状态
st = LuxCore.initialstates(rng, bucket)

# 单时间步输入
x_t = rand(Float32, n_inputs)

# 调用 Lux 接口: (x_t, ps, st) → (y_t, st_new)
y_t, st_new = bucket(x_t, ps_lux, st)

# st_new.hydro_state 包含更新后的水文状态
```

---

## 序列计算 (HydroModels 接口)

### 单节点

```julia
# 输入: (n_inputs × n_timesteps)
input = rand(Float32, 2, 100)

# 运行
output = bucket(input, params, config)  # (n_outputs × n_timesteps)

# 带初始状态
initstates = ComponentVector(
    states = (
        snowpack = [0.5f0],
        soilwater = [10.0f0]
    )
)
output = bucket(input, params, config; initstates=initstates)
```

### 多节点

```julia
# 输入: (n_inputs × n_nodes × n_timesteps)
input_multi = rand(Float32, 2, 5, 100)

# 运行 (每个节点独立维护状态)
output_multi = bucket(input_multi, params, config)  # (n_outputs × n_nodes × n_timesteps)
```

---

## 计算流程

### 单个时间步的内部流程

```julia
# 时间步 t
function single_step(bucket, x_t, hydro_state, ps)
    # 1. 计算通量
    flux_input = [hydro_state; x_t]
    fluxes = flux_network(flux_input, ps.flux)
    
    # 2. 更新状态
    state_input = [hydro_state; fluxes]
    new_hydro_state = state_network(state_input, ps.state)
    
    # 3. 计算输出
    output = output_network(fluxes, ps.output)
    
    return output, new_hydro_state
end
```

### 完整序列的计算

```julia
# 循环迭代 (类似 RNN)
for t in 1:n_timesteps
    x_t = input[:, t]
    y_t, hydro_state = single_step(bucket, x_t, hydro_state, ps)
    outputs[t] = y_t
end
```

---

## 应用场景

### 1. 纯数据驱动建模

```julia
# 直接从观测数据学习水文过程
bucket = create_neural_bucket(
    name = :data_driven,
    n_inputs = 5,
    n_states = 3,
    n_outputs = 1,
    hidden_size = 32,
    inputs = [:prcp, :temp, :humid, :wind, :rad],
    states = [:state1, :state2, :state3],
    outputs = [:discharge]
)
```

### 2. 混合建模

```julia
# 物理过程 (HydroBucket)
physical_bucket = @hydrobucket :physical begin
    fluxes = begin
        # 物理通量...
    end
    dfluxes = begin
        # 状态微分...
    end
end

# 神经网络修正 (NeuralBucket)
correction_bucket = create_neural_bucket(
    name = :correction,
    n_inputs = 2,
    n_states = 1,
    n_outputs = 1,
    inputs = [:model_output, :observation],
    states = [:correction_state],
    outputs = [:correction]
)

# 在 HydroModel 中组合使用
```

### 3. 学习复杂状态更新规则

```julia
# 使用深度网络学习非线性状态转换
bucket = NeuralBucket(;
    name = :complex,
    flux_network = Chain(
        Dense(6 => 64, relu),
        Dense(64 => 32, tanh),
        Dense(32 => 10)
    ),
    state_network = Chain(
        Dense(14 => 64, relu),
        Dense(64 => 32, relu),
        Dense(32 => 4)
    ),
    output_network = Dense(10 => 1),
    n_inputs = 2,
    n_states = 4,
    n_outputs = 1,
    inputs = [:input1, :input2],
    states = [:s1, :s2, :s3, :s4],
    outputs = [:output]
)
```

---

## 训练 NeuralBucket

### 使用 Optimization.jl

```julia
using Optimization, OptimizationOptimisers

# 定义损失函数
function loss_fn(ps_vec, data)
    inputs, targets = data
    
    # 重构参数
    ps_axes = getaxes(ComponentVector(ps_lux))
    ps_lux_reconstructed = ComponentVector(ps_vec, ps_axes)
    params = ComponentVector(
        nns = (
            bucket_flux = ps_lux_reconstructed.flux,
            bucket_state = ps_lux_reconstructed.state,
            bucket_output = ps_lux_reconstructed.output
        )
    )
    
    # 前向传播
    predictions = bucket(inputs, params, config)
    
    # 计算损失
    return sum((predictions .- targets).^2)
end

# 训练
optf = OptimizationFunction(loss_fn, Optimization.AutoZygote())
optprob = OptimizationProblem(optf, ComponentVector(ps_lux), (train_inputs, train_targets))
result = solve(optprob, Adam(0.001), maxiters=1000)
```

---

## 与 Zygote 的兼容性

NeuralBucket 完全支持 Zygote 自动微分：

```julia
using Zygote

# 包装为可微函数
function forward(ps_vec)
    ps_lux = ComponentVector(ps_vec, ps_axes)
    params = ComponentVector(nns = (...))
    output = bucket(input, params, config)
    return sum(output)
end

# 计算梯度
grads = gradient(forward, ComponentVector(ps_lux))
```

---

## 优势

### 1. 结构优雅

- ✅ 直接模拟 ODE 的三阶段过程
- ✅ 每个网络有明确的职责
- ✅ 类似于传统 HydroBucket 的结构

### 2. 灵活性

- ✅ 可以使用任意 Lux 层
- ✅ 可以添加 Dropout, BatchNorm, etc.
- ✅ 可以使用不同的激活函数

### 3. Lux 集成

- ✅ 原生参数管理
- ✅ 无缝自动微分
- ✅ GPU 加速支持

### 4. 可组合性

- ✅ 可以与 HydroBucket 混合使用
- ✅ 可以作为 HydroModel 的一部分
- ✅ 接口完全兼容

---

## 便捷函数

### `create_neural_bucket`

创建标准 NeuralBucket with 隐藏层：

```julia
bucket = create_neural_bucket(
    name = :my_bucket,
    n_inputs = 2,
    n_states = 2,
    n_outputs = 1,
    n_fluxes = 3,           # 默认 = n_states
    hidden_size = 16,       # 隐藏层大小
    inputs = [...],
    states = [...],
    outputs = [...],
    flux_activation = tanh,
    state_activation = identity,
    output_activation = identity
)
```

### `create_simple_neural_bucket`

创建最简单的线性 NeuralBucket (用于调试):

```julia
simple_bucket = create_simple_neural_bucket(
    name = :simple,
    n_inputs = 2,
    n_states = 1,
    n_outputs = 1,
    inputs = [:prcp, :temp],
    states = [:storage],
    outputs = [:flow]
)
```

---

## 最佳实践

### 1. 网络大小选择

- **隐藏层大小**:
  - 小问题: 8-16
  - 中等问题: 16-32
  - 大问题: 32-64

- **状态数量**:
  - 简单问题: 1-2
  - 中等问题: 2-4
  - 复杂问题: 4-8

### 2. 激活函数选择

- **Flux Network**: `tanh` 或 `relu`
- **State Network**: `identity` 或 `relu`
- **Output Network**: `identity` (根据需要)

### 3. 初始化

```julia
# 使用固定随机种子确保可重复性
rng = Random.default_rng()
Random.seed!(rng, 42)
ps = LuxCore.initialparameters(rng, bucket)
```

### 4. 数据标准化

```julia
# 标准化输入
mean_input = mean(train_inputs, dims=2)
std_input = std(train_inputs, dims=2)
normalized_input = (inputs .- mean_input) ./ std_input
```

---

## 常见问题

### Q1: NeuralBucket 和之前的 RNN-based bucket 有什么区别？

**A**: 完全不同！

- **之前的设计** (错误理解): 直接嵌入 LSTM/GRU cell
- **现在的设计** (正确): 使用三个独立的神经网络模拟 HydroBucket 的 ODE 结构

### Q2: 为什么需要三个网络？

**A**: 模拟 HydroBucket 的三个阶段：
1. Flux network → `@hydroflux` (计算通量)
2. State network → `@stateflux` (更新状态)
3. Output network → 提取输出

### Q3: 状态数量如何选择？

**A**: 取决于问题复杂度：
- 简单问题 (如单层土壤): 1-2 states
- 中等问题 (如雪+土壤): 2-4 states
- 复杂问题 (如多层系统): 4-8 states

### Q4: 可以不用 hidden layer 吗？

**A**: 可以！使用 `create_simple_neural_bucket` 创建线性变换版本，适合简单问题或调试。

### Q5: 如何与传统 HydroBucket 混合使用？

**A**: 在 `HydroModel` 中组合：

```julia
model = @hydromodel :hybrid begin
    physical_bucket  # 传统 HydroBucket
    neural_bucket    # NeuralBucket
end
```

---

## 总结

### 核心设计理念

✅ **ODE ≈ RNN**: 利用结构相似性  
✅ **三阶段**: 通量 → 状态 → 输出  
✅ **Lux 集成**: 原生参数管理和自动微分  
✅ **接口一致**: 与 HydroModels 无缝集成  

### 适用场景

✅ 纯数据驱动的水文建模  
✅ 混合物理-机器学习模型  
✅ 学习复杂的状态转换规则  
✅ 迁移学习应用  

### 下一步

1. 参考 `docs/notebook/test_neural_bucket_ode_like.jl` 查看完整示例
2. 根据具体问题调整网络架构
3. 使用 Optimization.jl 训练模型
4. 与传统方法对比评估

---

**版本**: HydroModels.jl v0.5  
**最后更新**: 2025-01-09  
**设计**: 基于 ODE-RNN 结构相似性

