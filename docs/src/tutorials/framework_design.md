# Framework Design

## Basic Structure

介绍这个框架中主要是由Flux,Element和Model三个核心模块所构成的,然后需要分别介绍这三个模块:

Flux适用于做什么的:需要分别介绍HydroFlux和StateFlux, NeuralFlux放在神经网络嵌入研究中介绍

Element就是包含了多个Flux,他是将这些Flux整合到一起进行执行常微分方程求解与各种通量一同计算的计算模块
包括Bucket和Route两种模块:
- Bucket模块是水文模型中的某个计算模块,例如融雪计算模块或土壤含水的计算模块,通常根据水量平衡原理结合模块的输入输出通量构建常微分方程进行求解,如...,一般支持单节点和多节点水文模型的计算
- Route模块是水文模型在产流计算后的河网汇流计算,根据Flux建立汇流过程的计算公式,并根据aggresive函数描述多节点之间的节点汇流过程
  
Model是将上述类型全部整合在一起的类型,他使得这些模块能够一同进行计算和参数优化

model将这些整合在一起的方式是基于Flux和element这些模块具有相同的调用一致性(所有模块都具有callable能力和相同的输入接口)

## 执行过程

1.数据需求

输入数据类型是特征*时间的矩阵输入;
模型参数通常以ComponentVector形式存储
初始状态

2.运行设置

flux和element的运行设置

就是设置一些计算所需的设置
solver: 可以使用DifferentialEquations.jl中的求解器
interp: 可以使用DataInterpolations.jl中的插值方法
timeidx: 是设置求解输出的时间索引

model的运行设置

model是将flux和element整合起来的东西,需要针对不同flux或element提供运行设置,config=config
