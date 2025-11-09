# NeuralBucket Guide - RNN-style Neural Hydrological Modeling

## 概述

`NeuralBucket` 是 HydroModels.jl v0.5 新增的神经网络建模组件，它将 Lux.jl 的循环神经网络（RNN）与水文建模框架无缝集成。

## 主要特点

### 1. 基于 RNN 的序列建模

- ✅ 支持 LSTM、GRU 和自定义 RNN 单元
- ✅ 自动维护隐藏状态，处理序列数据
- ✅ 不需要定义微分方程
- ✅ 默认使用离散时间步进

### 2. 与 Lux.jl 完全集成

- ✅ 继承 `Lux.AbstractLuxLayer`
- ✅ 实现标准 Lux 接口
- ✅ 支持 Lux 参数和状态管理
- ✅ 兼容 Lux 生态系统工具

### 3. 与 HydroModels 接口一致

- ✅ 使用 `ComponentVector` 参数格式
- ✅ 支持 `HydroConfig` 配置系统
- ✅ 支持单节点和多节点计算
- ✅ 可集成到 `HydroModel` 中

---

## 快速开始

### 创建 LSTM-based NeuralBucket

```julia
using HydroModels, Lux, ComponentArrays, Random

# 使用便捷函数创建 LSTM bucket
bucket = create_lstm_bucket(
    name = :rainfall_runoff,
    n_inputs = 2,      # 输入变量数（如降水、温度）
    n_hidden = 32,     # LSTM 隐藏层大小
    n_outputs = 1,     # 输出变量数（如径流）
    inputs = [:prcp, :temp],
    outputs = [:runoff]
)

# 初始化参数
rng = Random.default_rng()
ps = LuxCore.initialparameters(rng, bucket)

# 准备参数（HydroModels 格式）
params = ComponentVector(
    nns = (
        rainfall_runoff_lstm = ps.cell,
        rainfall_runoff_output = ps.output
    )
)

# 准备输入数据 (n_inputs × n_timesteps)
input = rand(Float32, 2, 100)

# 运行模型
config = HydroConfig()
output = bucket(input, params, config)
```

---

## 详细使用方法

### 1. 使用便捷函数

#### LSTM Bucket

```julia
lstm_bucket = create_lstm_bucket(
    name = :my_lstm,
    n_inputs = 3,
    n_hidden = 64,
    n_outputs = 2,
    inputs = [:var1, :var2, :var3],
    outputs = [:out1, :out2]
)
```

#### GRU Bucket

```julia
gru_bucket = create_gru_bucket(
    name = :my_gru,
    n_inputs = 2,
    n_hidden = 32,
    n_outputs = 1,
    inputs = [:prcp, :temp],
    outputs = [:flow]
)
```

### 2. 自定义 RNN 单元

```julia
using Lux

# 创建自定义 RNN 单元
custom_cell = Chain(
    LSTMCell(3 => 64),
    LSTMCell(64 => 32),
    name = :custom_lstm
)

# 创建输出层
output_layer = Chain(
    Dense(32 => 16, relu),
    Dense(16 => 2),
    name = :custom_output
)

# 构建 NeuralBucket
bucket = NeuralBucket(
    name = :custom_bucket,
    cell = custom_cell,
    output_layer = output_layer,
    inputs = [:input1, :input2, :input3],
    outputs = [:output1, :output2]
)
```

### 3. 多层 LSTM

```julia
# 深度 LSTM 网络
deep_lstm = Chain(
    LSTMCell(2 => 128),
    LSTMCell(128 => 64),
    LSTMCell(64 => 32),
    name = :deep_lstm
)

output_layer = Dense(32 => 1, name = :output)

bucket = NeuralBucket(
    name = :deep_model,
    cell = deep_lstm,
    output_layer = output_layer,
    inputs = [:prcp, :pet],
    outputs = [:flow]
)
```

---

## 多节点计算

NeuralBucket 支持多节点（空间分布式）计算：

```julia
# 创建 bucket
bucket = create_lstm_bucket(
    name = :multi_site,
    n_inputs = 2,
    n_hidden = 16,
    n_outputs = 1,
    inputs = [:prcp, :temp],
    outputs = [:flow],
    hru_types = collect(1:5)  # 5 个节点
)

# 准备多节点输入 (n_inputs × n_nodes × n_timesteps)
n_nodes = 5
n_timesteps = 100
input = rand(Float32, 2, n_nodes, n_timesteps)

# 初始化参数
ps = LuxCore.initialparameters(rng, bucket)
params = ComponentVector(
    nns = (
        multi_site_lstm = ps.cell,
        multi_site_output = ps.output
    )
)

# 运行（每个节点独立维护隐藏状态）
output = bucket(input, params, config)  # (n_outputs × n_nodes × n_timesteps)
```

---

## 与传统 HydroBucket 的对比

| 特性 | HydroBucket | NeuralBucket |
|------|------------|--------------|
| 建模方式 | 微分方程（ODE） | 循环神经网络（RNN） |
| 求解器 | 需要（Mutable/Immutable/ODE） | 不需要（离散步进） |
| 状态管理 | 显式状态变量 | RNN 隐藏状态 |
| 参数 | 物理参数 | 神经网络权重 |
| 可解释性 | 高（物理过程） | 低（黑盒） |
| 灵活性 | 受物理约束 | 高度灵活 |
| 数据需求 | 中等 | 大量数据 |
| 训练方式 | 参数优化 | 梯度下降 |

---

## 应用场景

### 1. 纯数据驱动建模

当缺乏物理过程知识或物理模型难以构建时：

```julia
# 直接从数据学习
lstm_model = create_lstm_bucket(
    name = :data_driven,
    n_inputs = 5,      # 多个气象变量
    n_hidden = 64,
    n_outputs = 1,     # 径流预测
    inputs = [:prcp, :temp, :humid, :wind, :rad],
    outputs = [:discharge]
)
```

### 2. 混合建模（物理+神经网络）

结合物理过程模型和神经网络：

```julia
# 物理过程 bucket
@variables prcp temp snowpack soilwater
@parameters Smax Qmax

physical_bucket = @hydrobucket :physical begin
    fluxes = begin
        # 物理过程...
    end
    dfluxes = begin
        @stateflux soilwater ~ ...
    end
end

# 神经网络 bucket（用于残差修正）
correction_bucket = create_lstm_bucket(
    name = :correction,
    n_inputs = 3,  # 物理模型输出 + 观测
    n_hidden = 16,
    n_outputs = 1,  # 修正量
    inputs = [:model_output, :prcp, :temp],
    outputs = [:correction]
)

# 可以在后处理中结合两者
```

### 3. 参数估计网络

使用神经网络动态估计模型参数：

```julia
# 参数估计 bucket
param_estimator = create_lstm_bucket(
    name = :param_est,
    n_inputs = 4,      # 气象条件
    n_hidden = 32,
    n_outputs = 3,     # 估计的参数
    inputs = [:prcp, :temp, :pet, :soil_moisture],
    outputs = [:param1, :param2, :param3]
)

# 然后在物理模型中使用估计的参数
```

### 4. 实时序列预测

RNN 自然适合序列预测任务：

```julia
forecaster = create_gru_bucket(
    name = :forecast,
    n_inputs = 6,      # 当前状态 + 预报输入
    n_hidden = 128,
    n_outputs = 24,    # 未来24小时预报
    inputs = [:current_state, :forecast_prcp, :forecast_temp, ...],
    outputs = [:hour_1, :hour_2, ..., :hour_24]
)
```

---

## 训练 NeuralBucket

NeuralBucket 参数可以通过标准的深度学习训练流程进行优化：

```julia
using Optimization, OptimizationOptimisers

# 定义损失函数
function loss_fn(params_vec, data)
    inputs, targets = data
    
    # 重构参数
    params = ComponentVector(
        nns = (
            bucket_lstm = params_vec.cell,
            bucket_output = params_vec.output
        )
    )
    
    # 前向传播
    predictions = bucket(inputs, params, config)
    
    # 计算损失（如 MSE）
    return sum((predictions .- targets).^2)
end

# 准备数据
train_data = (train_inputs, train_targets)

# 初始化参数
ps_init = LuxCore.initialparameters(rng, bucket)
ps_vec = ComponentVector(ps_init)

# 设置优化器
opt = Adam(0.001)

# 定义优化问题
optf = OptimizationFunction(loss_fn, Optimization.AutoZygote())
optprob = OptimizationProblem(optf, ps_vec, train_data)

# 训练
result = solve(optprob, opt, maxiters=1000)
```

---

## 与 Zygote 自动微分的兼容性

NeuralBucket 完全支持 Zygote 自动微分：

```julia
using Zygote

# 定义可微分的前向传播
function forward(params_vec)
    params = ComponentVector(...)
    output = bucket(input, params, config)
    return sum(output)  # 标量输出用于梯度计算
end

# 计算梯度
grads = gradient(forward, ps_vec)
```

---

## 性能优化建议

### 1. 批处理

```julia
# 将多个时间序列打包为批次
batch_input = cat(seq1, seq2, seq3, dims=3)  # (n_inputs × n_nodes × n_timesteps)
```

### 2. 预分配内存

```julia
# 对于重复调用，可以预分配输出数组
output_buffer = similar(input, (n_outputs, n_timesteps))
```

### 3. GPU 加速

```julia
using CUDA

# 将数据和参数移至 GPU
input_gpu = cu(input)
params_gpu = ComponentVector(cu.(ps))

# 注意：需要 Lux 的 GPU 支持
```

---

## 最佳实践

### 1. 隐藏层大小选择

- 小数据集: 16-32
- 中等数据集: 32-64
- 大数据集: 64-128+

### 2. 避免过拟合

- 使用 Dropout（在 Lux 层中添加）
- 正则化（L2 权重衰减）
- 早停（监控验证损失）

### 3. 数据标准化

```julia
# 标准化输入数据
mean_input = mean(train_inputs, dims=2)
std_input = std(train_inputs, dims=2)
normalized_input = (inputs .- mean_input) ./ std_input
```

### 4. 模型保存和加载

```julia
using JLD2

# 保存参数
@save "model_params.jld2" ps

# 加载参数
@load "model_params.jld2" ps
```

---

## 常见问题

### Q1: NeuralBucket 和 NeuralFlux 有什么区别？

**A**: 
- `NeuralFlux`: 无状态的神经网络，用于替换单个通量计算
- `NeuralBucket`: 有状态的 RNN，维护隐藏状态处理序列

### Q2: 可以在 HydroModel 中混合使用 HydroBucket 和 NeuralBucket 吗？

**A**: 理论上可以，但需要注意接口匹配。NeuralBucket 通常用作独立组件或后处理器。

### Q3: NeuralBucket 支持双向 RNN 吗？

**A**: 可以通过自定义 cell 实现：
```julia
bidirectional_cell = Chain(
    BidirectionalRNN(LSTMCell(n_in => n_hidden)),
    name = :birnn
)
```

### Q4: 如何处理不同长度的序列？

**A**: 可以使用填充（padding）或掩码（masking）技术，需要在自定义 cell 中实现。

---

## 示例代码库

完整的示例代码请参考：
- `docs/notebook/test_neural_bucket.jl` - 基础测试和示例
- `docs/notebook/train_neural_bucket.jl` - 训练示例（待创建）

---

## 总结

`NeuralBucket` 为 HydroModels.jl 带来了强大的深度学习建模能力，特别适合：

✅ 数据驱动的水文建模  
✅ 混合物理-机器学习模型  
✅ 实时序列预测  
✅ 参数估计和校准  
✅ 迁移学习应用  

它与 Lux.jl 生态系统完全集成，同时保持了与 HydroModels.jl 的接口一致性。

---

**版本**: v0.5  
**最后更新**: 2025-01-09

