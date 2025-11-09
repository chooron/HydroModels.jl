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
        "tutorials" => [
            "Framework Design" => "tutorials/framework_design.md",
            "Neural Network Embedding" => "tutorials/neuralnetwork_embeding.md",
            "Distribute Modeling" => "tutorials/distribute_modeling.md",
            "Optimal Parameters" => "tutorials/optimimal_parameters.md"
        ],
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


# 部署配置
deploydocs(
    repo = "github.com/chooron/HydroModels.jl",
    devbranch = "main",
    push_preview = true,
    target = "build_en"  # 确保这里指定了正确的构建目录
)