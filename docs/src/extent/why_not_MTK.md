# HydroModels.jl与ModelingToolkit.jl

## 为什么不直接使用[ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl)

HydroModels.jl的设计理念事实上是与ModelingToolkit.jl完全一致的,两者都是使用Symbolics.jl进行符号计算的框架,ModelingToolkit.jl是一个更通用符号系统,支持多种问题的构建,同时针对特定领域[ModelingToolkitStandardLibrary.jl](https://github.com/SciML/ModelingToolkitStandardLibrary.jl)提供了基础组件的支持.
但是我之前在ModelingToolkit.jl的使用中存在一些问题:

- 尽管ModelingToolkit.jl同样能够根据符号编程表达水文模型的常微分方程,但对水文模型中单位线计算,汇流过程计算等计算过程的支持却并不是特别理想
- 在同一种模块中,比如蒸发计算公式和土壤计算模块,不同区域所需要构建的模块会因为公式的不同存在差异,需要分别构建计算公式支撑模型模型构建
- 水文模型通常包括多个常微分方程,使用一个ODESystem构建会相对混乱,使用多个ODESystem会相对复杂
- 同时对于多节点输入,空间汇流过程计算,ODESystem或无法直接支持.
- 对于神经网络模型的嵌入,尽管存在ModelingToolkitNeuralNets.jl,然而嵌入性却没有想象的那么简单
- 基于ModelingToolkit.jl对于自动微分的支持,尤其是Zygote.jl我在使用时也存在问题.

因此我决定自己构建一个模型库,参考了符号编程的模型构建方式,使其更能够支撑水文模型的一些建模需求.

## 未来的工作

不可否认ModelingToolkit.jl是一个非常好的框架,在满足了水文模型的需求之后,我会尽量将HydroModels.jl中的模块能够成为ModelingToolkitStandardLibrary.jl中相似的模块,并继承ModelingToolkit.jl提供的AbstractSystem类的支持,提供一致的功能支持,从而兼容更多SciML生态的计算功能.