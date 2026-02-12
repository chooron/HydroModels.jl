# HydroModels 参数初始化和优化功能总结

## 完成的功能

### 1. `get_initial_params` 函数 (src/utils.jl)

自动从组件中提取并初始化参数：

**特性：**
- 从组件的 `infos` 元数据中获取参数名称
- 传统参数：初始化为 1.0（简单且稳定）
- 神经网络参数：通过 Lux 随机初始化
- 使用多重派发支持不同组件类型
- 支持 RNG 种子以确保可重现性

**设计理念：**
- 保持简单：不从 component 中提取参数边界（避免增加类型复杂度）
- `infos` 只存储 Symbol，无法获取原始 Num 对象
- 用户可以在 OptimizationProblem 中提供自定义 `initial_params`

**示例：**
```julia
# 简单 flux
@variables prcp q
@parameters k
flux = @hydroflux q ~ k * prcp

params = get_initial_params(flux)
# 返回: ComponentVector(params=(k=1.0,))

# 神经网络 flux
@variables prcp temp et
nn = Chain(Dense(2 => 10, tanh), Dense(10 => 1), name=:et_net)
neural_flux = @neuralflux et ~ nn([prcp, temp])

params = get_initial_params(neural_flux)
# 返回: ComponentVector(nns=(et_net=ComponentVector(...),))

# 混合组件
params = get_initial_params(hybrid_model)
# 返回: ComponentVector(params=(k=1.0, b=1.0), nns=(et_net=...,))
```

### 2. `get_nn_initial_params` 函数 (src/utils.jl)

仅提取神经网络参数：

```julia
nn_params = get_nn_initial_params(model)
# 返回: (et_net = ComponentVector(...), ...)
```

### 3. 简化的 OptimizationProblem 接口 (ext/HydroModelsOptimizationExt)

**自动初始化参数**，支持自定义初始值：

```julia
# 默认初始化（推荐）
prob = OptimizationProblem(
    model, input, target;
    fixed_params=ComponentVector(params=(k=0.5,)),
    lb_pas=[0.1, 0.0],
    ub_pas=[5.0, 10.0]
)

# 自定义初始参数
initial_params = ComponentVector(
    params=(k=2.0, b=1.5),
    nns=(et_net=custom_nn_params,)
)
prob = OptimizationProblem(
    model, input, target;
    initial_params=initial_params,
    lb_pas=[0.1, 0.0],
    ub_pas=[5.0, 10.0]
)
```

### 4. 层级参数固定

支持 ComponentVector 的层级结构：

```julia
# 固定传统参数
fixed_params=ComponentVector(params=(k=0.5, b=2.0))

# 固定神经网络参数
nn_params = get_nn_initial_params(model)
fixed_params=ComponentVector(nns=nn_params)

# 固定特定神经网络层
fixed_params=ComponentVector(
    params=(k=0.5,),
    nns=(
        et_net=ComponentVector(
            layer_1=nn_params.et_net.layer_1  # 只固定第一层
        ),
    )
)
```

## 测试结果

所有测试通过：

```
[Test 1] HydroFlux with parameters
  ✓ Parameters initialized to 1.0

[Test 2] NeuralFlux with neural network
  ✓ Has neural network parameters (41 params)

[Test 3] get_nn_initial_params
  ✓ Extracted NN params: true

[Test 4] Reproducibility with RNG
  ✓ Reproducible: true

[Test 5] Neural network reproducibility
  ✓ Neural network params reproducible
```

## 主要改进

1. **更简洁的 API**：不需要手动传入参数
2. **简单稳定**：传统参数默认为 1.0，避免复杂的边界提取
3. **多重派发**：使用 Julia 的多重派发而非类型判断
4. **从 infos 获取**：直接从组件元数据获取参数名称
5. **可重现性**：支持 RNG 种子
6. **灵活性**：支持自定义 `initial_params`

## 设计决策

### 为什么不从 component 提取参数边界？

1. **类型复杂度**：`infos` 只存储 Symbol，无法获取原始 Num 对象
2. **梯度计算**：在类型中存储额外信息会影响梯度计算性能
3. **简单性**：默认值 1.0 简单且稳定
4. **灵活性**：用户可以通过 `initial_params` 提供自定义初始值

### 推荐的工作流程

```julia
# 1. 创建模型
model = HydroModel([...])

# 2. 获取初始参数（可选，用于检查）
params = get_initial_params(model)

# 3. 自定义初始值（可选）
# 方式A：使用默认值
prob = OptimizationProblem(model, input, target; lb_pas=..., ub_pas=...)

# 方式B：提供自定义初始值
initial_params = ComponentVector(params=(k=2.0, b=1.5))
prob = OptimizationProblem(model, input, target;
    initial_params=initial_params, lb_pas=..., ub_pas=...)

# 4. 求解
sol = solve(prob, optimizer)
```

## 文件修改

- `src/utils.jl`: 添加 `get_initial_params` 和 `get_nn_initial_params`（简化版）
- `src/bucket.jl`: 添加 `neural_fluxes` 字段到 HydroBucket
- `ext/HydroModelsOptimizationExt/HydroModelsOptimizationExt.jl`: 添加 `initial_params` 支持
- `ext/HydroModelsOptimizationExt/README.md`: 更新文档
- `ext/HydroModelsOptimizationExt/examples.jl`: 更新示例
- `test/miscellaneous/test_get_initial_params.jl`: 添加测试
