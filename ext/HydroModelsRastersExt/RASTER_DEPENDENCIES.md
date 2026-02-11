# 栅格处理扩展依赖库详解

本文档详细说明了 `HydroModelsRastersExt` 扩展所需的依赖库及其用途。

## 核心依赖

### 1. Rasters.jl ⭐ (必需)

**安装:**
```julia
using Pkg
Pkg.add("Rasters")
```

**版本要求:** 0.11.x 或 0.12.x

**功能特性:**
- 🗺️ **多格式支持**: GeoTIFF, NetCDF, GRIB, HDF5等
- 📐 **坐标系统**: 完整的CRS和投影支持
- 🔄 **重采样**: 最近邻、双线性、三次等多种方法
- ⚡ **高性能**: 懒加载和内存映射
- 🔧 **易用性**: 与DimensionalData.jl集成，直观的API

**主要用途:**
```julia
# 读取栅格文件
dem = Raster("dem.tif")

# 访问数据和元数据
data = dem.data          # 获取数组数据
x_dim = dims(dem, X)     # 获取X维度
crs = crs(dem)           # 获取坐标参考系统

# 重采样
resampled = resample(dem, to=target_dem, method=:bilinear)

# 裁剪
cropped = dem[X(Between(0, 100)), Y(Between(0, 100))]

# 保存
write("output.tif", dem)
```

**官方文档:** https://rafaqz.github.io/Rasters.jl/stable/

---

### 2. ArchGDAL.jl (可选，推荐)

**安装:**
```julia
using Pkg
Pkg.add("ArchGDAL")
```

**功能特性:**
- 📚 **格式最全**: 支持100+种栅格和矢量格式
- 🌍 **完整GIS工具**: GDAL的全功能Julia封装
- 🔧 **底层控制**: 适合需要精细控制的场景
- 🔗 **与Rasters.jl配合**: 可作为Rasters的后端

**主要用途:**
```julia
using ArchGDAL as AG

# 读取栅格
AG.read("dem.tif") do dataset
    # 获取栅格信息
    width = AG.width(dataset)
    height = AG.height(dataset)
    
    # 读取波段数据
    band = AG.getband(dataset, 1)
    data = AG.read(band)
    
    # 获取地理变换参数
    geotransform = AG.getgeotransform(dataset)
end

# 坐标转换
source_crs = AG.importEPSG(4326)  # WGS84
target_crs = AG.importEPSG(3857)  # Web Mercator
```

**官方文档:** https://yeesian.com/ArchGDAL.jl/stable/

---

## 辅助依赖（根据需求添加）

### 3. GeoInterface.jl

**安装:**
```julia
Pkg.add("GeoInterface")
```

**用途:**
- 统一的地理空间数据接口
- 不同GIS包之间的互操作性

**示例:**
```julia
using GeoInterface

# 检查几何类型
GeoInterface.geomtype(dem)

# 获取坐标
GeoInterface.coordinates(geometry)
```

---

### 4. Proj.jl

**安装:**
```julia
Pkg.add("Proj")
```

**用途:**
- 坐标系统转换
- 投影变换
- 大地测量计算

**示例:**
```julia
using Proj

# 创建转换
trans = Proj.Transformation("EPSG:4326", "EPSG:3857")

# 转换坐标
x_new, y_new = trans(lon, lat)
```

---

### 5. CoordinateTransformations.jl

**安装:**
```julia
Pkg.add("CoordinateTransformations")
```

**用途:**
- 通用坐标变换
- 仿射变换
- 组合变换

---

### 6. Images.jl / ImageFiltering.jl

**安装:**
```julia
Pkg.add("Images")
Pkg.add("ImageFiltering")
```

**用途:**
- 图像处理算法（可用于栅格数据）
- 滤波和卷积操作
- 形态学操作

**示例:**
```julia
using ImageFiltering

# 对DEM应用高斯滤波（平滑）
smoothed = imfilter(dem_data, Kernel.gaussian(3))

# 计算坡度（使用Sobel算子）
grad_x = imfilter(dem_data, Kernel.sobel()[1])
grad_y = imfilter(dem_data, Kernel.sobel()[2])
```

---

## D8流向分析相关

### 目前状态
Julia生态系统中**没有专门的D8流向分析包**，但我们可以使用以下组合：

#### 自己实现D8算法（推荐）✅
```julia
# 在HydroModelsRastersExt中已实现
flow_dir = compute_d8_flow_direction(dem)
```

**优势:**
- ✅ 完全控制算法细节
- ✅ 与HydroRoute完美集成
- ✅ 可自定义扩展（如D-infinity）
- ✅ 性能优化空间大

#### 使用Python互操作（备选）
如果需要更复杂的地形分析，可以通过PyCall调用Python包：

```julia
using PyCall

# 调用Python的richdem
richdem = pyimport("richdem")
dem_py = richdem.LoadGDAL("dem.tif")
flow_dir = richdem.FlowDirection(dem_py, method="D8")
```

**可用的Python包:**
- **richdem**: 高性能DEM分析
- **pysheds**: 流域分析
- **whitebox**: 地形分析工具集

---

## 扩展功能的可选依赖

### 7. NCDatasets.jl
用于NetCDF格式（气象数据常用）:
```julia
Pkg.add("NCDatasets")
```

### 8. GeoStats.jl
地统计分析和插值:
```julia
Pkg.add("GeoStats")
```

### 9. GeoDataFrames.jl
地理数据框（类似GeoPandas）:
```julia
Pkg.add("GeoDataFrames")
```

---

## 快速开始安装脚本

### 最小安装（仅核心功能）
```julia
using Pkg
Pkg.add("Rasters")
```

### 推荐安装（完整GIS功能）
```julia
using Pkg
packages = ["Rasters", "ArchGDAL", "GeoInterface", "Proj"]
Pkg.add(packages)
```

### 完整安装（包括辅助工具）
```julia
using Pkg
packages = [
    "Rasters",
    "ArchGDAL", 
    "GeoInterface",
    "Proj",
    "Images",
    "ImageFiltering",
    "NCDatasets"
]
Pkg.add(packages)
```

---

## 性能考虑

### 大型DEM处理

1. **使用懒加载**
```julia
# Rasters.jl支持懒加载
dem = Raster("large_dem.tif", lazy=true)
```

2. **分块处理**
```julia
# 按块处理大型DEM
for block in blocks(dem, size=(1000, 1000))
    flow_dir_block = compute_d8_flow_direction(block)
    # 处理块...
end
```

3. **使用多线程**
```julia
# 启动Julia时使用多线程
# julia --threads 8

using Base.Threads
@threads for i in 1:nblocks
    # 并行处理各块
end
```

---

## 故障排除

### 常见问题

**1. GDAL未安装**
```julia
# 解决方法：使用Conda安装GDAL
using Pkg
Pkg.add("Conda")
using Conda
Conda.add("gdal")
```

**2. 坐标系统不匹配**
```julia
# 检查CRS
crs(dem)

# 转换CRS
dem_reprojected = reproject(dem, crs=EPSG(4326))
```

**3. NoData值处理**
```julia
# 替换NoData
dem_clean = replace_missing(dem, NaN)

# 或使用掩码
mask = .!ismissing.(dem)
```

---

## 总结

| 功能 | 推荐包 | 优先级 |
|------|--------|--------|
| 读取栅格 | Rasters.jl | ⭐⭐⭐ 必需 |
| 格式转换 | ArchGDAL.jl | ⭐⭐ 推荐 |
| 坐标转换 | Proj.jl | ⭐⭐ 推荐 |
| D8流向分析 | 自己实现 | ⭐⭐⭐ 已实现 |
| 图像处理 | Images.jl | ⭐ 可选 |
| 统计分析 | GeoStats.jl | ⭐ 可选 |

---

## 参考资源

- [Rasters.jl文档](https://rafaqz.github.io/Rasters.jl/stable/)
- [ArchGDAL.jl文档](https://yeesian.com/ArchGDAL.jl/stable/)
- [Julia地理空间生态系统](https://github.com/JuliaGeo)
- [GDAL官方文档](https://gdal.org/)
- [D8流向算法论文](https://doi.org/10.1016/0098-3004(84)90020-7)

---

**最后更新:** 2024-11-10

