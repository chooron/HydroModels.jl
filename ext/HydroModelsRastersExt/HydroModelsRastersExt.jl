module HydroModelsRastersExt

using Rasters
using HydroModels

"""
    compute_d8_flow_direction(dem::Raster)

根据DEM数据计算D8流向矩阵。

D8编码规则（与HydroRoute兼容）：
- 1: 东 (E)
- 2: 东南 (SE)
- 4: 南 (S)
- 8: 西南 (SW)
- 16: 西 (W)
- 32: 西北 (NW)
- 64: 北 (N)
- 128: 东北 (NE)

# 参数
- `dem::Raster`: 输入的DEM栅格数据

# 返回
- `flow_direction::Matrix{Int}`: D8流向矩阵，使用上述编码
"""
function compute_d8_flow_direction(dem::Raster)
    # 获取DEM数据矩阵
    dem_data = dem.data
    rows, cols = size(dem_data)
    
    # 初始化流向矩阵
    flow_dir = zeros(Int, rows, cols)
    
    # D8方向编码和偏移量
    # 顺序：E, SE, S, SW, W, NW, N, NE
    d8_codes = [1, 2, 4, 8, 16, 32, 64, 128]
    d8_offsets = [
        (0, 1),   # E
        (1, 1),   # SE
        (1, 0),   # S
        (1, -1),  # SW
        (0, -1),  # W
        (-1, -1), # NW
        (-1, 0),  # N
        (-1, 1)   # NE
    ]
    
    # 获取像元大小（用于计算斜向距离）
    cellsize_x = abs(step(dims(dem, X)))
    cellsize_y = abs(step(dims(dem, Y)))
    
    # 对每个像元计算流向
    for i in 1:rows
        for j in 1:cols
            # 跳过NoData值
            if ismissing(dem_data[i, j]) || isnan(dem_data[i, j])
                flow_dir[i, j] = 0
                continue
            end
            
            center_elev = dem_data[i, j]
            max_slope = -Inf
            flow_direction = 0
            
            # 检查8个邻居
            for (idx, (di, dj)) in enumerate(d8_offsets)
                ni, nj = i + di, j + dj
                
                # 检查边界
                if ni < 1 || ni > rows || nj < 1 || nj > cols
                    continue
                end
                
                # 跳过NoData邻居
                if ismissing(dem_data[ni, nj]) || isnan(dem_data[ni, nj])
                    continue
                end
                
                # 计算坡度
                neighbor_elev = dem_data[ni, nj]
                
                # 计算距离（对角线距离更长）
                distance = if abs(di) + abs(dj) == 2
                    sqrt(cellsize_x^2 + cellsize_y^2)  # 对角线
                else
                    max(cellsize_x, cellsize_y)  # 正交
                end
                
                slope = (center_elev - neighbor_elev) / distance
                
                # 更新最大坡度方向
                if slope > max_slope
                    max_slope = slope
                    flow_direction = d8_codes[idx]
                end
            end
            
            flow_dir[i, j] = flow_direction
        end
    end
    
    return flow_dir
end

"""
    extract_positions_from_raster(raster::Raster, valid_mask::Function=x->!ismissing(x))

从栅格数据中提取有效像元的位置索引。

# 参数
- `raster::Raster`: 输入栅格
- `valid_mask::Function`: 判断像元是否有效的函数，默认为非缺失值

# 返回
- `positions::Vector{Tuple{Int,Int}}`: 有效像元的(row, col)位置列表
"""
function extract_positions_from_raster(raster::Raster, valid_mask::Function=x->!ismissing(x))
    data = raster.data
    rows, cols = size(data)
    positions = Tuple{Int,Int}[]
    
    for i in 1:rows
        for j in 1:cols
            if valid_mask(data[i, j])
                push!(positions, (i, j))
            end
        end
    end
    
    return positions
end

"""
    match_raster_to_dem(raster::Raster, dem::Raster; method=:bilinear)

将栅格数据重采样以匹配DEM的分辨率和范围。

# 参数
- `raster::Raster`: 需要匹配的栅格数据
- `dem::Raster`: 参考DEM
- `method::Symbol`: 重采样方法，可选 :bilinear, :nearest, :cubic

# 返回
- `matched_raster::Raster`: 匹配后的栅格数据
"""
function match_raster_to_dem(raster::Raster, dem::Raster; method=:bilinear)
    # 使用Rasters.jl的resample功能
    # 将raster重采样到dem的维度
    matched = resample(raster, to=dem, method=method)
    return matched
end

"""
    create_hydroroute_from_dem(dem::Raster; name::Symbol=:dem_route)

从DEM直接创建HydroRoute对象。

# 参数
- `dem::Raster`: 输入DEM栅格
- `name::Symbol`: 路由名称

# 返回
- `route::HydroRoute`: 配置好的HydroRoute对象
"""
function create_hydroroute_from_dem(dem::Raster; name::Symbol=:dem_route)
    # 计算流向矩阵
    flow_dir = compute_d8_flow_direction(dem)
    
    # 提取有效位置
    positions = extract_positions_from_raster(dem)
    
    # 创建聚合函数
    aggr_func = HydroModels.build_aggr_func(flow_dir, positions)
    
    # 提取HRU类型（假设每个有效像元都是同一类型）
    hru_types = ones(Int, length(positions))
    
    # 注意：这里需要用户提供fluxes和dfluxes定义
    # 这只是一个框架示例
    println("警告：需要提供flux和dflux定义来完整构建HydroRoute")
    println("流向矩阵大小: $(size(flow_dir))")
    println("有效位置数量: $(length(positions))")
    
    return (flow_dir=flow_dir, positions=positions, aggr_func=aggr_func, hru_types=hru_types)
end

# 导出主要函数
export compute_d8_flow_direction
export extract_positions_from_raster
export match_raster_to_dem
export create_hydroroute_from_dem

end

