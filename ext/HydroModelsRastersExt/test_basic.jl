# 基础测试脚本
# 测试HydroModelsRastersExt的基本功能

println("测试 HydroModelsRastersExt 基本功能")
println("=" ^ 60)

# 测试1: 模块加载
println("\n测试1: 检查依赖...")
try
    using HydroModels
    println("✓ HydroModels加载成功")
catch e
    println("✗ HydroModels加载失败: $e")
end

try
    using Rasters
    println("✓ Rasters.jl加载成功")
catch e
    println("✗ Rasters.jl加载失败: $e")
    println("  请运行: using Pkg; Pkg.add(\"Rasters\")")
end

# 测试2: 创建测试DEM
println("\n测试2: 创建测试DEM...")
try
    # 创建一个简单的倾斜平面DEM
    dem_data = Float64[
        10.0 9.5 9.0 8.5 8.0;
        9.5 9.0 8.5 8.0 7.5;
        9.0 8.5 8.0 7.5 7.0;
        8.5 8.0 7.5 7.0 6.5;
        8.0 7.5 7.0 6.5 6.0
    ]
    
    dem = Raster(
        dem_data,
        dims=(X(1.0:1.0:5.0), Y(1.0:1.0:5.0))
    )
    
    println("✓ DEM创建成功")
    println("  大小: $(size(dem))")
    println("  数值范围: $(extrema(dem_data))")
catch e
    println("✗ DEM创建失败: $e")
end

# 测试3: D8流向分析
println("\n测试3: D8流向分析...")
try
    flow_dir = compute_d8_flow_direction(dem)
    println("✓ D8流向计算成功")
    println("  流向矩阵:")
    display(flow_dir)
    
    # 检查流向编码是否合理
    valid_codes = [0, 1, 2, 4, 8, 16, 32, 64, 128]
    all_valid = all(in(valid_codes), flow_dir)
    if all_valid
        println("✓ 所有流向编码都有效")
    else
        println("✗ 发现无效的流向编码")
    end
catch e
    println("✗ D8流向计算失败: $e")
    println("  可能原因: 扩展未正确加载")
    println("  请确保在加载HydroModels后加载Rasters")
end

# 测试4: 位置提取
println("\n测试4: 位置提取...")
try
    positions = extract_positions_from_raster(dem)
    println("✓ 位置提取成功")
    println("  有效像元数: $(length(positions))")
    println("  预期数量: $(length(dem_data))")
    
    if length(positions) == length(dem_data)
        println("✓ 位置数量正确")
    else
        println("✗ 位置数量不匹配")
    end
catch e
    println("✗ 位置提取失败: $e")
end

# 测试5: 聚合函数创建
println("\n测试5: 聚合函数创建...")
try
    aggr_func = HydroModels.build_aggr_func(flow_dir, positions)
    println("✓ 聚合函数创建成功")
    
    # 测试聚合函数
    test_input = ones(Float64, length(positions))
    result = aggr_func(test_input)
    
    println("✓ 聚合函数运行成功")
    println("  输入长度: $(length(test_input))")
    println("  输出长度: $(length(result))")
    
    if length(result) == length(positions)
        println("✓ 输出维度正确")
    else
        println("✗ 输出维度不正确")
    end
catch e
    println("✗ 聚合函数测试失败: $e")
end

# 测试6: 重采样功能
println("\n测试6: 栅格重采样...")
try
    # 创建一个低分辨率的栅格
    low_res_data = Float64[
        10.0 8.0;
        8.0 6.0
    ]
    
    low_res = Raster(
        low_res_data,
        dims=(X(1.0:2.0:3.0), Y(1.0:2.0:3.0))
    )
    
    # 重采样到DEM的分辨率
    resampled = match_raster_to_dem(low_res, dem, method=:bilinear)
    
    println("✓ 重采样成功")
    println("  原始大小: $(size(low_res))")
    println("  重采样后大小: $(size(resampled))")
    
    if size(resampled) == size(dem)
        println("✓ 重采样尺寸匹配")
    else
        println("✗ 重采样尺寸不匹配")
    end
catch e
    println("✗ 重采样失败: $e")
end

# 测试总结
println("\n" * "=" ^ 60)
println("测试完成!")
println("=" ^ 60)
println("\n如果所有测试都通过，扩展已正确安装和配置。")
println("否则，请检查错误信息并安装所需的依赖包。")

