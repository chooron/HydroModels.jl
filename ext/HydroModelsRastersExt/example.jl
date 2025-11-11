# HydroModelsRastersExt 使用示例
# 
# 本文件展示如何使用栅格扩展进行DEM流向分析和HydroRoute构建

using HydroModels
using Rasters
using ComponentArrays

# ===================================================================
# 示例 1: 基础的D8流向分析
# ===================================================================

println("=" ^ 60)
println("示例 1: 从DEM计算D8流向")
println("=" ^ 60)

# 创建一个简单的测试DEM（实际使用时从文件读取）
function create_test_dem()
    # 创建一个小型的测试DEM (10x10)
    dem_data = [
        100.0 99.0 98.0 97.0 96.0 95.0 94.0 93.0 92.0 91.0;
        99.0 98.0 97.0 96.0 95.0 94.0 93.0 92.0 91.0 90.0;
        98.0 97.0 96.0 95.0 94.0 93.0 92.0 91.0 90.0 89.0;
        97.0 96.0 95.0 94.0 93.0 92.0 91.0 90.0 89.0 88.0;
        96.0 95.0 94.0 93.0 92.0 91.0 90.0 89.0 88.0 87.0;
        95.0 94.0 93.0 92.0 91.0 90.0 89.0 88.0 87.0 86.0;
        94.0 93.0 92.0 91.0 90.0 89.0 88.0 87.0 86.0 85.0;
        93.0 92.0 91.0 90.0 89.0 88.0 87.0 86.0 85.0 84.0;
        92.0 91.0 90.0 89.0 88.0 87.0 86.0 85.0 84.0 83.0;
        91.0 90.0 89.0 88.0 87.0 86.0 85.0 84.0 83.0 82.0
    ]
    
    # 创建Raster对象（带坐标信息）
    dem = Raster(
        dem_data,
        dims=(
            X(1.0:1.0:10.0),  # X坐标从1到10，步长1
            Y(1.0:1.0:10.0)   # Y坐标从1到10，步长1
        )
    )
    
    return dem
end

# 创建测试DEM
test_dem = create_test_dem()
println("DEM大小: ", size(test_dem))
println("DEM范围: X=$(dims(test_dem, X)), Y=$(dims(test_dem, Y))")

# 计算D8流向
flow_dir = compute_d8_flow_direction(test_dem)
println("\n流向矩阵 (前5x5):")
println(flow_dir[1:5, 1:5])

# 解释流向编码
println("\nD8流向编码:")
println("  32  64  128")
println("   ↖  ↑  ↗")
println("16 ← • → 1")
println("   ↙  ↓  ↘")
println("  8   4   2")

# ===================================================================
# 示例 2: 提取位置并创建聚合函数
# ===================================================================

println("\n" * "=" ^ 60)
println("示例 2: 提取位置并创建聚合函数")
println("=" ^ 60)

# 提取所有有效位置
positions = extract_positions_from_raster(test_dem)
println("有效像元数量: ", length(positions))
println("前5个位置: ", positions[1:5])

# 创建聚合函数（用于HydroRoute）
aggr_func = HydroModels.build_aggr_func(flow_dir, positions)
println("\n聚合函数已创建，类型: ", typeof(aggr_func))

# 测试聚合函数
test_outflow = ones(Float64, length(positions))  # 假设每个位置流出1单位水量
routed_inflow = aggr_func(test_outflow)
println("测试流出向量长度: ", length(test_outflow))
println("路由后流入向量长度: ", length(routed_inflow))

# ===================================================================
# 示例 3: 从真实DEM文件读取并处理
# ===================================================================

println("\n" * "=" ^ 60)
println("示例 3: 从真实文件读取DEM")
println("=" ^ 60)

# 注释掉的代码 - 实际使用时取消注释
# dem_file = "path/to/your/dem.tif"
# if isfile(dem_file)
#     # 读取DEM
#     dem = Raster(dem_file)
#     println("DEM信息:")
#     println("  大小: ", size(dem))
#     println("  分辨率: ", step.(dims(dem)))
#     println("  范围: ", bounds(dem))
#     
#     # 计算流向
#     flow_direction = compute_d8_flow_direction(dem)
#     
#     # 保存流向矩阵
#     write("flow_direction.tif", Raster(flow_direction, dims=dims(dem)))
#     println("\n流向矩阵已保存到 flow_direction.tif")
# else
#     println("DEM文件不存在: $dem_file")
# end

println("（需要提供真实DEM文件路径）")

# ===================================================================
# 示例 4: 匹配不同分辨率的栅格数据
# ===================================================================

println("\n" * "=" ^ 60)
println("示例 4: 栅格数据重采样")
println("=" ^ 60)

# 创建一个不同分辨率的测试栅格（如降雨数据）
rainfall_data = rand(5, 5) .* 100  # 5x5的降雨数据
rainfall_raster = Raster(
    rainfall_data,
    dims=(
        X(1.0:2.0:9.0),  # 分辨率为2（比DEM粗）
        Y(1.0:2.0:9.0)
    )
)

println("降雨栅格大小: ", size(rainfall_raster))
println("DEM大小: ", size(test_dem))

# 将降雨数据重采样到DEM的分辨率
rainfall_matched = match_raster_to_dem(rainfall_raster, test_dem, method=:bilinear)
println("重采样后降雨栅格大小: ", size(rainfall_matched))

# ===================================================================
# 示例 5: 构建完整的HydroRoute（概念示例）
# ===================================================================

println("\n" * "=" ^ 60)
println("示例 5: 构建HydroRoute（概念示例）")
println("=" ^ 60)

# 使用create_hydroroute_from_dem获取基础组件
dem_info = create_hydroroute_from_dem(test_dem, name=:test_watershed)

println("\n从DEM提取的信息:")
println("  流向矩阵大小: ", size(dem_info.flow_dir))
println("  有效位置数量: ", length(dem_info.positions))
println("  HRU类型: ", length(dem_info.hru_types))

# 构建完整的HydroRoute需要定义fluxes和dfluxes
# 这里是一个概念性示例

println("\n要构建完整的HydroRoute，需要定义:")
println("  1. fluxes: 水文通量方程（如出流 = k * 蓄水量）")
println("  2. dfluxes: 状态变化方程（如蓄水量变化 = 流入 - 流出）")
println("  3. 参数: 模型参数（如流速系数k）")

# 示例代码结构（需要实际的flux定义）:
# @hydroroute :my_route begin
#     fluxes = begin
#         @hydroflux outflow ~ k * storage
#     end
#     dfluxes = begin
#         @stateflux storage ~ inflow - outflow
#     end
#     hru_types = dem_info.hru_types
#     aggr_func = dem_info.aggr_func
# end

# ===================================================================
# 示例 6: 流向统计和可视化
# ===================================================================

println("\n" * "=" ^ 60)
println("示例 6: 流向统计")
println("=" ^ 60)

# 统计各个方向的像元数量
d8_codes = [1, 2, 4, 8, 16, 32, 64, 128]
d8_names = ["东(E)", "东南(SE)", "南(S)", "西南(SW)", "西(W)", "西北(NW)", "北(N)", "东北(NE)"]

println("\n流向分布:")
for (code, name) in zip(d8_codes, d8_names)
    count = sum(flow_dir .== code)
    percentage = count / length(flow_dir) * 100
    println("  $name ($code): $count 像元 ($(round(percentage, digits=2))%)")
end

# 统计无流向的像元（边界或洼地）
no_flow_count = sum(flow_dir .== 0)
println("  无流向(0): $no_flow_count 像元")

println("\n" * "=" ^ 60)
println("所有示例完成！")
println("=" ^ 60)

# ===================================================================
# 附加功能: 自定义验证函数
# ===================================================================

println("\n附加示例: 使用自定义验证函数提取位置")

# 只提取高程大于95的位置
high_positions = extract_positions_from_raster(test_dem, x -> !ismissing(x) && x > 95.0)
println("高程 > 95的位置数量: ", length(high_positions))
println("占总位置的比例: $(round(length(high_positions)/length(positions)*100, digits=2))%")

