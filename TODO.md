## TODO

- [ ] 添加测试文件(模型优化与求解)
- [X] 添加说明文档
- [ ] 完善checkers和wrappers
- [X] unithydrograph的macro化表达
- [ ] Multi-Nodes在支持多节点梯度计算时存在问题,仅能使用ZygoteVJP,无法利用EnzymeVJP,但是Multi-Nodes的应用是绝对必要的,因此需要其他的解决方法
- [ ] 添加gradient的test

## Changes

1. Added macro-based construction methods to reduce the complexity of hydrological model development
2. Added CUDA support for multi-node computations
3. Implemented neural network structures based on the Recurrence computation process from Lux.jl, making it easier to support multi-node hydrological model calculations (still under development)

## Upcoming Features and Improvements

- [ ] Add test files for model optimization and solving
- [ ] Create comprehensive documentation
- [ ] Enhance checkers and wrappers
- [X] Develop macro-based expressions for unit hydrographs
- [ ] Address issues with multi-node gradient calculations (currently limited to ZygoteVJP, unable to utilize EnzymeVJP, but multi-node applications are essential, requiring alternative solutions)
- [ ] 添加BMI接口支持
- [ ] 新增config类型，使其支持读取，存储功能
- [ ] 仅有参数作为输入时，HydroFlux计算的结果的时间维度只有1
- [ ] unithydro需要支持更多的计算
- [ ] 需要添加文档，使其能够明白如何判断水文模型是否正确
- [ ] ifelse, abs在toexpr2中还不能支持，此外认为需要专门针对判断语句构建一个flux
- [ ] 通过HydroModels.jl实现敏感性分析（结构与参数上）
- [X] 如何实现在hydroflux中同时定义多个计算公式并权重相加
- [ ] 如何实现同级的hydrobucket权重叠加，考虑将bucket的func进行拼接
- [ ] 集成计算，bucket与model
- [ ] bucket需要分为多维和单维的
- [ ] bucket中优化参数超过100个时会使用ReverseDiffAdjoint这个敏感性算法，这时候优化会出现错误，难以解决

## BUG
- [X] bucket 中如果有neuralflux,打印的时候会报错
- [X] neuralbucket仍有问题
- [X] 还是得拆分开, 把DifferentialEquations独立出来
- [X] unit hydro 也存在问题
- [X] @neuralfux 构建时可以接收一个归一化的参数
- [X] styidx这个参数可能重要性不大
- [X] HydroBucket用于多节点计算过少,考虑抛弃这个功能,将多节点计算主要放在Lux-based的基础上
- [ ] 统一接口为 input, params, config; input, params两种
