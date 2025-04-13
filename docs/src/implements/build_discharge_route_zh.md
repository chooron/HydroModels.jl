# 构建河道汇流模型

## [discharge route model](https://gmd.copernicus.org/articles/14/7795/2021/)汇流模型介绍

河道汇流模型是描述水流在河道中运动过程的数学模型，主要用于计算河道中的流量变化。下面给出了一个简单的河道汇流模型的数学表达式:

```math
\begin{aligned}
\Delta S_{rf,n} &= Q_{rf,n-1} - Q_{rf,n} && (1) \\
Q_{rf,n} &= (S_{rf,n} + Q_{rf,n-1}) \cdot \frac{1}{\text{LAG}_{rf} + 1} && (2) \\
Q_{rf,0} &= \sum Q_{rf,up} && (3) \\
Q_{rf} &= Q_{rf,c_{max}} + (R_{sw} + R_{gw}) \cdot 0.001 \cdot A_{gc} && (4)
\end{aligned}
```

其中,$n$表示河段编号，取值范围为$\{1,...,c_{max}\}$;$\Delta S_{rf,n}$表示第 n 段河道的蓄水量变化;$Q_{rf,n}$表示第 n 段河道的出流量;$S_{rf,n}$表示第 n 段河道的蓄水量;$\text{LAG}_{rf}$ 表示河道汇流滞时参数;$Q_{rf,up}$表示上游河道的入流量;$R_{sw}$和$R_{gw}$分别表示地表径流和地下径流;$A_{gc}$表示栅格面积;$c_{max}$表示河道最大分段数.

公式(1)描述了河道蓄水量的变化，公式(2)表示河道出流量的计算，公式(3)给出了河道入流量的计算方法，公式(4)则是最终河道出流量的计算公式。
这种构建方式是根据每个 hydrological response unit(HRU)的上游输入确定当前的 HRU 入流量,这个确定方式可以由两种路由函数来表示. 当前 HRU 的出流量则是以当前的河道蓄水量$S_{rf,n}$计算得到,并于 HRU 的 R*{sw}和 R*{gw}进行相加(经面积转换)得到最终的出流量$Q_{rf}$.

## HydroModels.jl 框架下的河道汇流模型构建

### 计算公式的改造

HydroModels.jl 框架中为这种基于河道蓄水状态计算时段出流量的方法以 HydroRoute 类型进行表示,通过构建 ODE 方程的方式对这个问题进行连续性求解.与原公式实现有所不同的是,,并认为在计算$Q_{rf,n}$时不存在瞬时变化,所以将公式(2)转换为如下式子:

```math
\begin{aligned}
Q_{rf,n} &= S_{rf,n} \cdot \frac{1}{\text{LAG}_{rf} + 1} && (5)
\end{aligned}
```

式中,去掉了$Q_{rf,n-1}$项.因为认为,在持续性变化中,$S_{rf,n}$是一直处于变化状态的,所以在计算$Q_{rf,n}$时不需要考虑$Q_{rf,n-1}$的变化.此外于$S_{rf,n}$与$Q_{rf,n}$的量纲存在差异,同时如果将式中如果考虑$Q_{rf,n}$会额外增加一个状态变量,这个变量是一个瞬时值,也不适合作为状态变量来使用. 在 HydroRoute 构造中,其首要的就是针对出流计算公式进行表达,与 HydroBucket 一样,HydroRoute 可以接受 HydroFlux 以此表示$Q_{rf,n}$的计算公式,如下所示:

```julia
# 定义涉及的参数和变量
@variables q q_routed s_river
@parameters lag
# 定义
rflux = HydroFlux([q, s_river] => [q_routed], [lag], exprs=[s_river / (1 + lag) + q])
```

### 流量累积函数表示

获取 HRU 上游流量输入的方式是构成汇流模型的另一个部分,通常可以设置 Aggregation 函数来获取每个 HRU 的上游入流流量. 这个 Aggregation 函数可以由 Grid-Based 和 Vector-Based 两种方式来表示. 在 HydroModels.jl 框架中分别对应了两种 HydroRoute 的构造方式: GridRoute 和 VectorRoute.

#### 表示 Vector-Based 的河道汇流模型

Vector-Based 的河道汇流模型是以流域子流域划分结果及其连接的拓扑关系为基础,通过构建 adjacency 矩阵得到流量累积函数:

```julia
# 构建河网拓扑结构
network = DiGraph(9)
add_edge!(network, 1, 2)
add_edge!(network, 2, 5)
add_edge!(network, 3, 5)
add_edge!(network, 4, 5)
add_edge!(network, 5, 8)
add_edge!(network, 6, 9)
add_edge!(network, 7, 8)
add_edge!(network, 8, 9)
vroute = HydroModels.VectorRoute(rfunc=rflux, rstate=s_river, network=network)
```

代码中是采用 Graphs.jl 的 DiGraph 数据结构来表示流域拓扑,作为模型的输入参数,连同构建的 rflux 和 s_river,一同完成 Vector-Based 的河道汇流模型的构造.

#### 表示 Grid-Based 的河道汇流模型

Grid-Based 的河道汇流模型是 HRU 之间的连同关系所表示的, 同常是使用 d8 流向矩阵来表示其连同关系,以此 Grid-Based 的河道汇流模型的构造便是根据 d8 流向矩阵和对应的 HRU 坐标来进行构造的:

```julia
# 构建河网拓扑结构
flwdir = [1 4 8; 1 4 4; 1 1 2]
positions = [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]
groute = HydroModels.GridRoute(rfunc=rflux, rstate=s_river, flwdir=flwdir, positions=positions)
```

代码中是采用矩阵数据表示流向矩阵和HRU 坐标,作为模型的输入参数,连同构建的 rflux 和 s_river,一同完成 Grid-Based 的河道汇流模型的构造.