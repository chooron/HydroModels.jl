# YAML Executor for HydroModels

YAML Executor 提供了一种通过 YAML 配置文件直接执行水文模型的便捷方式，无需编写 Julia 代码即可完成模型配置、数据加载和执行。

## 功能特性

- ✅ 从 YAML 文件加载完整的模型配置
- ✅ 支持 CSV 格式的驱动数据
- ✅ 支持 NetCDF 格式的驱动数据（通过 Rasters 扩展）
- ✅ 灵活的变量名称映射
- ✅ 支持从文件或直接在 YAML 中指定参数值
- ✅ 自动处理数据排序和格式转换

## 快速开始

### 1. 基本用法

```julia
using HydroModels
using YAML

# 从 YAML 配置文件执行模型
output = execute_from_yaml("model_config.yaml")

# 或者获取所有组件
output, model, config, data = execute_from_yaml("model_config.yaml", return_components=true)
```

### 2. YAML 配置文件结构

```yaml
version: "1.0"
schema: "hydromodels"

# 参数定义
parameters:
  Smax:
    description: "Maximum soil moisture"
    units: "mm"
    default: 100.0
    bounds: [50.0, 500.0]

# 模型组件
components:
  - type: HydroBucket
    name: soil
    fluxes:
      - formula: "baseflow ~ Qmax * exp(-f * (Smax - soilwater))"
    state_fluxes:
      - formula: "soilwater ~ Prcp - evap - baseflow"

# 模型组合
model:
  type: HydroModel
  name: simple_model
  components: [soil]

# 求解器配置
config:
  solver: MutableSolver
  interpolator: DirectInterpolation
  min_value: 1.0e-6

# 驱动数据配置
data:
  type: csv
  path: "forcing.csv"
  variables:
    Prcp: Precipitation    # 模型变量: 文件列名
    PET: PotentialET
  time_column: date

# 参数值配置
parameters_data:
  type: csv
  path: "params.csv"

# 初始状态
initial_states:
  soilwater: 50.0
```

## 数据格式

### CSV 格式

#### 驱动数据 (forcing.csv)
```csv
date,Precipitation,PotentialET,Temperature
1,5.2,3.1,20.5
2,0.0,3.5,21.2
3,12.5,2.8,19.8
...
```

#### 参数数据 (params.csv)
```csv
Smax,Qmax,f,Df
100.0,10.0,2.0,0.5
```

### NetCDF 格式

对于空间分布式模型，可以使用 NetCDF 格式：

```yaml
data:
  type: netcdf
  path: "forcing_data.nc"
  variables:
    Prcp: precipitation
    PET: pet
  positions:
    - row: 10
      col: 20
    - row: 15
      col: 25
```

对于集总式模型（空间平均）：

```yaml
data:
  type: netcdf
  path: "forcing_data.nc"
  variables:
    Prcp: precipitation
    PET: pet
```

## 参数配置方式

### 方式 1: 从 CSV 文件加载

```yaml
parameters_data:
  type: csv
  path: "parameters.csv"
```

### 方式 2: 直接在 YAML 中指定

```yaml
parameters_data:
  type: values
  values:
    Smax: 120.0
    Qmax: 8.5
    f: 2.5
    Df: 0.6
```

### 方式 3: 使用默认值

如果不提供 `parameters_data`，将使用 `parameters` 部分定义的 `default` 值。

## 变量名称映射

YAML Executor 支持灵活的变量名称映射，允许数据文件中的列名与模型变量名不同：

```yaml
data:
  type: csv
  path: "forcing.csv"
  variables:
    temp: Temperature      # 模型需要 temp，文件中是 Temperature
    prcp: Precipitation    # 模型需要 prcp，文件中是 Precipitation
    pet: PET              # 模型需要 pet，文件中是 PET
```

## 完整示例

查看 `examples/yaml_executor_usage.jl` 获取完整的使用示例，包括：

1. 创建示例数据
2. 配置 YAML 文件
3. 执行模型
4. 分析和可视化结果

运行示例：

```julia
include("examples/yaml_executor_usage.jl")
```

## 配置选项

### Solver 选项
- `MutableSolver`: 可变求解器（默认）
- `ImmutableSolver`: 不可变求解器
- `ODESolver`: ODE 求解器
- `DiscreteSolver`: 离散求解器

### Interpolator 选项
- `DirectInterpolation` / `ConstantInterpolation`: 常数插值
- `LinearInterpolation` / `EnzymeInterpolation`: 线性插值

## 注意事项

1. **变量顺序**: 驱动数据矩阵的列顺序由 `variables` 字典的键顺序决定
2. **数据类型**: 所有数值数据会自动转换为 `Float64`
3. **缺失值**: CSV 中的缺失值会导致错误，请确保数据完整
4. **NetCDF 支持**: 需要安装 `Rasters` 和 `NCDatasets` 包

## 依赖包

基础功能：
- `HydroModels`
- `YAML`
- `CSV`
- `DataFrames`
- `ComponentArrays`

NetCDF 支持（可选）：
- `Rasters`
- `NCDatasets`

## 测试

运行测试：

```julia
include("test/base/run_yaml_executor.jl")
```

## 故障排除

### 错误: "Column 'XXX' not found in CSV file"
- 检查 CSV 文件的列名是否正确
- 确保 `variables` 映射中的文件列名与实际列名匹配

### 错误: "Variable XXX not found in data"
- 检查 `variables` 映射是否包含所有模型需要的变量
- 确保变量名与公式中使用的变量名一致

### 错误: "No parameter values provided"
- 提供 `parameters_data` 部分，或
- 在 `parameters` 定义中提供 `default` 值

## 贡献

欢迎提交问题和改进建议！
