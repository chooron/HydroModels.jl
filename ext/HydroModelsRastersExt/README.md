# HydroModelsRastersExt

这是HydroModels的栅格数据处理扩展，专门用于：
- 读取和处理DEM（数字高程模型）数据
- 执行D8流向分析
- 生成流向矩阵作为HydroRoute的输入
- 支持栅格数据的重采样和匹配

## 依赖

此扩展需要安装以下包：

```julia
using Pkg
Pkg.add("Rasters")
Pkg.add("ArchGDAL")  # 可选，用于更多格式支持
```

## 使用示例

### 1. 从DEM文件计算流向矩阵

```julia
using HydroModels
using Rasters

# 读取DEM数据
dem = Raster("path/to/your/dem.tif")

# 计算D8流向
flow_dir = compute_d8_flow_direction(dem)

# 提取有效位置
positions = extract_positions_from_raster(dem)

# 创建聚合函数
aggr_func = HydroModels.build_aggr_func(flow_dir, positions)
```

### 2. 匹配栅格数据到DEM

```julia
# 读取另一个栅格数据（如降雨、土地利用等）
rainfall = Raster("path/to/rainfall.tif")

# 将rainfall重采样到DEM的分辨率和范围
rainfall_matched = match_raster_to_dem(rainfall, dem)
```

### 3. 创建完整的HydroRoute

```julia
using HydroModels

# 从DEM获取基础信息
dem_info = create_hydroroute_from_dem(dem, name=:my_watershed)

# 定义水文通量
@hydroroute :watershed_route begin
    fluxes = begin
        @hydroflux outflow ~ k * storage
    end
    dfluxes = begin
        @stateflux storage ~ inflow - outflow
    end
    hru_types = dem_info.hru_types
    aggr_func = dem_info.aggr_func
end
```

## D8流向编码

本扩展使用的D8流向编码与HydroRoute兼容：

```
   32  64  128
    ↖  ↑  ↗
16 ← • → 1
    ↙  ↓  ↘
   8   4   2
```

- 1: 东 (E)
- 2: 东南 (SE)
- 4: 南 (S)
- 8: 西南 (SW)
- 16: 西 (W)
- 32: 西北 (NW)
- 64: 北 (N)
- 128: 东北 (NE)

## 功能特性

### compute_d8_flow_direction
- 自动处理NoData值
- 考虑对角线距离的坡度计算
- 处理边界条件
- 返回与HydroRoute兼容的流向矩阵

### extract_positions_from_raster
- 提取有效像元位置
- 支持自定义验证函数
- 返回位置索引列表

### match_raster_to_dem
- 支持多种重采样方法（最近邻、双线性、三次）
- 自动匹配坐标系统
- 处理分辨率和范围差异

## 注意事项

1. **坐标系统**：确保DEM和其他栅格数据使用相同的坐标系统，或使用`match_raster_to_dem`进行转换
2. **内存使用**：大型DEM文件可能占用大量内存，考虑分块处理
3. **NoData处理**：算法会自动跳过NoData值，流向矩阵中用0表示
4. **边界效应**：边缘像元可能没有完整的8个邻居，算法会适当处理

## 扩展计划

未来可能添加的功能：
- [ ] D-infinity流向算法
- [ ] 流量累积计算
- [ ] 河网提取
- [ ] 子流域划分
- [ ] 并行计算支持
- [ ] 洼地填充算法

