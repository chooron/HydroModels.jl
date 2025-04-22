using DataInterpolations
using Enzyme

matrices = rand(3, 3, 10)
t = collect(1:10)

# 创建线性插值函数
multi_interps = LinearInterpolation(reshape(matrices, :, 10), t)

function tmp_func3(inputs, p)
    prcp = inputs[1]
    temp = inputs[2]
    lday = inputs[3]
    Tmin = p[1]
    Df = p[2]
    Tmax = p[3]
    snowfall = prcp .* (temp .- Tmin) .* Df
    melt = lday .* (temp .- Tmax)
    return stack([snowfall .- melt], dims=1) |> sum
end

# 使用Enzyme计算梯度
function f(x, p)
    inputs = eachslice(reshape(multi_interps(x), 3, 3), dims=1)
    prcp = inputs[1]
    temp = inputs[2]
    lday = inputs[3]
    Tmin = p[1]
    Df = p[2]
    Tmax = p[3]
    snowfall = prcp .* (temp .- Tmin) .* Df
    melt = lday .* (temp .- Tmax)
    return stack([snowfall .- melt], dims=1) |> sum
end

ps = [2.3, 3.4, 5.0]
ps_d = zeros(eltype(ps), length(ps))
# 计算梯度
grad = Enzyme.autodiff(Reverse, f, Active, Active(3.5), Duplicated(ps, ps_d))[1]

    # prcp = inputs[1]
    # temp = inputs[2]
    # lday = inputs[3]
    # Tmin = p[1]
    # Df = p[2]
    # Tmax = p[3]
    # snowfall = prcp .* (temp .- Tmin) .* Df
    # melt = lday .* (temp .- Tmax)
    # return stack([snowfall .- melt], dims=1) |> sum