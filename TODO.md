## TODO

- [ ] Add more examples
- [X] Using the NamedTupleAdaptor, and remove the convert_to_ntp kwargs in component running
- [ ] Write more error messages
- [ ] macro building
- [X] The package is too heavy, which needs to be splitted into small packages
- [ ] Flux output need constraint
- [ ] Add NoneMeanFlux
- [ ] **GPU support (limited by scalar index)**
    - [X] DataInterpolations.jl is not support for now, supprot in v7.1.0
- [ ] meta should store variables and parameters
- [X] compact problem
- [ ] ensemble model
- [ ] base merge for component meta
- [ ] Neural Flux make the gradient computation slow a lot 
    - [ ] @btime gradient计算: NeuralFlux (84 us) 略大于 常规计算 (69 us)
    - [ ] @btime gradient计算: TotalModel (90s) 显著大于 常规计算 (23 s)
    - [ ] 原因:
    1. 多个bucket遍历插值
- [X] 针对multiple nodes模型参数的最佳设置方式, 需要考虑计算效率, 考虑采用 ComponentVector(params=(p1=[2,3,4], p2=[3,3,4]))
- [ ] 添加gradient的test
- [ ] 添加参数检查, 参数默认设置的功能
- [ ] 添加说明文档
- [ ] 如何适应多节点计算
- [ ] 参考Lux.jl对于Dense的实现或能够替换函数的构建,这样或许能够使模型支持GPU的计算