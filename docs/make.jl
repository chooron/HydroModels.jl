push!(LOAD_PATH, "../src/")
using Documenter
using HydroModels

makedocs(
    sitename = "HydroModels.jl",
    authors = "xin jing",
    format = Documenter.HTML(
        # 启用 pretty URLs，移除 .html 后缀
        # 设置文档的规范 URL
        canonical = "https://chooron.github.io/HydroModels.jl",
        # 设置资源文件
        assets = ["assets/icons.ico"],
        # 配置侧边栏
        collapselevel = 2,
        sidebar_sitename = true
    ),
    # 配置模块
    modules = [HydroModels],
    clean = true,
    doctest = false,
    linkcheck = true,
    source = "src",
    build = "build",
    # 配置页面结构
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Run a Bucket Model" => "tutorials/run_a_bucket.md",
            "Run ExpHydro Model" => "tutorials/run_a_exphydro_model.md"
        ],
        "Extent Content" => [
            "Modeling Framework Comparisons" => "docs/src/extent/framework comparision - en.md",
        ],
        
    ],
    # 其他选项
    checkdocs = :none,
)

# 部署配置
deploydocs(
    repo = "github.com/chooron/HydroModels.jl",
    devbranch = "main",
    push_preview = true,
    target = "build"  # 确保这里指定了正确的构建目录
)