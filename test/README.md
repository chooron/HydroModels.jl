# HydroModels.jl 测试指南

## 测试结构

测试文件组织如下：

```
test/
├── test_helpers.jl          # 共享的测试辅助函数和常量
├── run_tests.jl             # 测试运行脚本
├── runtests.jl              # 标准 Julia 测试入口
├── base/                    # 基础组件测试
│   ├── run_hydro_flux.jl    # Flux 组件测试
│   ├── run_neural_flux.jl   # 神经网络 Flux 测试
│   ├── run_single_bucket.jl # 单节点 Bucket 测试
│   ├── run_single_lumped.jl # 单节点完整模型测试
│   ├── run_multi_bucket.jl  # 多节点 Bucket 测试
│   ├── run_multi_lumped.jl  # 多节点完整模型测试
│   ├── run_unithydro.jl     # 单位线测试
│   ├── run_hydro_route.jl   # 路由组件测试
│   └── run_spatial_model.jl # 空间模型测试
├── gradient/                # 梯度计算测试
│   ├── run_hybrid_gradient.jl
│   └── run_hydro_gradient.jl
├── miscellaneous/           # 其他测试
└── models/                  # 测试模型定义
    ├── exphydro.jl
    ├── exphydrov2.jl
    └── m50.jl
```

## 运行测试

### 方法 1: 使用标准 Julia 测试

```julia
using Pkg
Pkg.test("HydroModels")
```

### 方法 2: 使用测试运行脚本

运行所有测试：
```bash
julia test/run_tests.jl
```

运行特定测试组：
```bash
julia test/run_tests.jl basic      # 基础组件测试
julia test/run_tests.jl single     # 单节点模型测试
julia test/run_tests.jl multi      # 多节点模型测试
julia test/run_tests.jl routing    # 路由组件测试
julia test/run_tests.jl gradient   # 梯度计算测试
julia test/run_tests.jl misc       # 其他测试
```

### 方法 3: 直接运行单个测试文件

```julia
using Pkg
Pkg.activate(".")

using Test
using HydroModels
using CSV, DataFrames, ComponentArrays, Graphs

# 加载测试辅助函数
include("test/test_helpers.jl")

# 定义 step_func
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

# 运行特定测试
include("test/base/run_single_bucket.jl")
```

## 测试辅助函数

`test_helpers.jl` 提供了以下常用函数：

### 数据加载

```julia
# 加载测试数据
input_ntp, input_mat, df = load_test_data(:exphydro, 1:100)
input_ntp, input_mat, df = load_test_data(:gr4j, 1:100)
input_ntp, input_mat, df = load_test_data(:m50, 1:100)
```

### 多节点数据准备

```julia
# 创建多节点输入 (variables × nodes × time)
input_arr = create_multinode_input(input_mat, num_nodes)

# 创建多节点参数
node_params = create_multinode_params(single_params, num_nodes)

# 创建多节点初始状态
node_states = create_multinode_states(single_states, num_nodes)
```

### 配置创建

```julia
# 创建测试配置
config = create_test_config(
    solver = MutableSolver,
    interpolator = Val(ConstantInterpolation),
    timeidx = 1:100,
    min_value = 1e-6,
    parallel = false
)
```

### 标准参数和状态

```julia
# ExpHydro 模型参数
EXPHYDRO_PARAMS = (
    Df = 2.674548848,
    Tmax = 0.175739196,
    Tmin = -2.092959084,
    Smax = 1709.46,
    f = 0.0167,
    Qmax = 18.47
)

# ExpHydro 初始状态
EXPHYDRO_STATES = (
    snowpack = 0.0,
    soilwater = 1303.0
)

# GR4J 模型参数
GR4J_PARAMS = (
    x1 = 320.11,
    x2 = 2.42,
    x3 = 69.63,
    x4 = 1.39
)

# GR4J 初始状态
GR4J_STATES = (
    soilwater = 235.97,
    routingstore = 45.47
)
```

## 接口变更说明

### v0.6.0 主要变更

1. **统一的 2D/3D 组件**：使用 `htypes` 参数区分单节点和多节点
   ```julia
   # 单节点 (2D: variables × time)
   bucket_2d = @hydrobucket :my_bucket begin
       fluxes = begin
           flux1
       end
   end

   # 多节点 (3D: variables × nodes × time)
   bucket_3d = @hydrobucket :my_bucket begin
       fluxes = begin
           flux1
       end
       htypes = [1, 1, 2, 2, 3]  # 参数共享索引
   end
   ```

2. **新的插值器接口**：
   - `ConstantInterpolation`: 常数插值（向上取整）
   - `LinearInterpolation`: 线性插值
   - 向后兼容别名：`DirectInterpolation` → `ConstantInterpolation`

3. **配置系统**：
   ```julia
   config = HydroConfig(
       solver = MutableSolver,
       interpolator = Val(ConstantInterpolation),
       timeidx = 1:100,
       min_value = 1e-6
   )
   ```

## 数据要求

测试需要以下数据文件：
- `data/exphydro/01013500.csv`
- `data/gr4j/sample.csv`
- `data/m50/01013500.csv`

确保这些文件存在于项目根目录的 `data/` 文件夹中。

## 故障排除

### 常见问题

1. **数据文件未找到**
   - 确保数据文件在正确的位置
   - 检查文件路径是否正确

2. **维度不匹配错误**
   - 单节点组件需要 2D 输入 (variables × time)
   - 多节点组件需要 3D 输入 (variables × nodes × time)
   - 检查 `htypes` 参数是否正确设置

3. **参数错误**
   - 单节点：`params = ComponentVector(params = (p1 = value1, ...))`
   - 多节点独立参数：每个参数是长度为 num_nodes 的向量
   - 多节点共享参数：每个参数是长度为 num_types 的向量，通过 htypes 索引

## 贡献测试

添加新测试时：
1. 使用 `test_helpers.jl` 中的辅助函数
2. 遵循现有测试的结构和命名约定
3. 添加清晰的测试描述和断言
4. 确保测试可以独立运行
