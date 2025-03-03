push!(LOAD_PATH, "../src/")
using Documenter
using HydroModels

# English Documentation
makedocs(
    sitename = "HydroModels.jl",
    authors = "xin jing",
    format = Documenter.HTML(
        # 启用 pretty URLs，移除 .html 后缀
        # 设置文档的规范 URL
        canonical = "https://chooron.github.io/HydroModels.jl/dev",
        # 配置侧边栏
        collapselevel = 2,
        sidebar_sitename = true,
        # Add icons
        assets = ["assets/icon.ico"],
        edit_link = "main"
    ),
    # 配置模块
    modules = [HydroModels],
    clean = true,
    doctest = false,
    linkcheck = true,
    source = "src",
    build = "build_en", 
    warnonly = true,
    # 配置页面结构
    pages = [
        "Home" => "index.md",
        "Get Started with HydroModels.jl" => "get_start_en.md",
        "Basic Concepts" => "concepts_en.md",
        "Model Implementations" => [
            "construct the ExpHydro Model" => "implements/build_exphydro_model_en.md",
            "construct the M50 Model" => "implements/build_m50_model_en.md",
            "construct the discharge route model" => "implements/build_discharge_route_en.md",
        ],
        "Extend Contents" => [
            "Why not using ModelingToolkit.jl directly" => "extent/why_not_MTK_en.md",
            "Framework Comparision" => "extent/framework_comparision_en.md",
        ]
    ]
)

# # Chinese Documentation
# makedocs(
#     sitename = "HydroModels.jl",
#     authors = "xin jing",
#     format = Documenter.HTML(
#         # 启用 pretty URLs，移除 .html 后缀
#         # 设置文档的规范 URL
#         canonical = "https://chooron.github.io/HydroModels.jl/dev",
#         # 配置侧边栏
#         collapselevel = 2,
#         sidebar_sitename = true
#     ),
#     # 配置模块
#     modules = [HydroModels],
#     clean = true,
#     doctest = false,
#     linkcheck = true,
#     source = "src",
#     build = "build_zh",
#     lang = "zh",
#     # 配置页面结构
#     pages = [
#         "主页" => "index.md",
#         "入门指南" => "getting_started.md",
#         "基本概念" => [
#             "蓄水模型" => "concepts/bucket_model.md",
#             "通量模型" => "concepts/flux_model.md",
#             "汇流模型" => "concepts/route_model.md",
#             "求解器模型" => "concepts/solver_model.md",
#             "分布式模型" => "concepts/spatial_model.md",
#             "集总式模型" => "concepts/lumped_model.md",
#         ],
#         "模型构建" => [
#             "构建 ExpHydro 模型" => "implements/build_exphydro_model_zh.md",
#             "构建 M50 模型" => "implements/build_m50_model_zh.md",
#             "构建河道汇流模型" => "implements/build_discharge_route_zh.md",
#         ],
#         "模型运行" => [
#             "运行蓄水模型" => "tutorials/run_a_bucket.md",
#             "运行 ExpHydro 模型" => "tutorials/run_a_exphydro_model.md"
#         ]
#     ]
# )

# 部署配置
deploydocs(
    repo = "github.com/chooron/HydroModels.jl",
    devbranch = "main",
    push_preview = true,
    target = "build_en"  # 确保这里指定了正确的构建目录
)

# deploydocs(
#     repo = "github.com/chooron/HydroModels.jl",
#     devbranch = "main",
#     push_preview = true,
#     target = "build_zh"  # 确保这里指定了正确的构建目录
# )