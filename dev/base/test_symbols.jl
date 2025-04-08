using SymbolicUtils.Code
using SymbolicUtils.Code: NameState,NaNMathFuns,arguments
using SymbolicUtils
using SymbolicUtils: operation
using Symbolics
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

@variables a b c d
@variables p1 p2 p3


assign_list = [
    Assignment(c, a * p1 + b * p2),
    Assignment(d, c * p1 + b * p3),
]
toexpr(assign_list[1])
var_names = [:a, :b]
ps_names = [:p1, :p2, :p3]
formula_exprs = [toexpr(a * p1 + b * p2), toexpr(c * p1 + b * p3)]
formula_targets = [:c, :d]

@generated function f_1(i, p)
    def_calls1 = [:($var = view(i, $var_idx)) for (var, var_idx) in zip(var_names, 1:length(var_names))]
    def_calls2 = [:($ps = view(p, $var_idx)) for (ps, var_idx) in zip(ps_names, 1:length(ps_names))]
    compute_calls = [:($target = @. $expr) for (target, expr) in zip(formula_targets, formula_exprs)]
    return_call = :(return [$(formula_targets...)])
    calls = reduce(vcat, [def_calls1, def_calls2, compute_calls, return_call])

    return Expr(:block,
        calls...
    )
end

function build_f1(input_names, param_names, outputs, exprs)
    def_calls1 = [:($nm = i[$idx]) for (idx, nm) in enumerate(input_names)]
    def_calls2 = [:($nm = p[$idx]) for (idx, nm) in enumerate(param_names)]
    compute_calls = [:($nm = @. $(expr)) for (nm, expr) in zip(outputs, exprs)]
    return_calls = :(return [$(outputs...)])
    return :(function (i, p)
        begin
            $(def_calls1...)
            $(def_calls2...)
            $(compute_calls...)
            $(return_calls)
        end
    end)
end

function build_f2()
    return :(function (i, p)
        a_ = i[1]
        b_ = i[2]
        p1_ = p[1]
        p2_ = p[2]
        p3_ = p[3]
        c_ = @__dot__((+)((*)(a_, p1_), (*)(b_, p2_)))
        d_ = @__dot__((+)((tanh)((*)(b_, p3_)), (sin)((*)(c_, p1_))))
        return [c_, d_]
    end)
end

input_names = [:a_, :b_]
param_names = [:p1_, :p2_, :p3_]
outputs = [:c_, :d_]

function SymbolicUtils.Code.function_to_expr(op::Union{typeof(tanh),typeof(sin)}, O, st)
    (get(st.rewrites, :nanmath, false) && op in NaNMathFuns) || return nothing
    name = nameof(op)
    fun = GlobalRef(NaNMath, name)
    args = map(Base.Fix2(toexpr, st), arguments(O))
    expr = :(@. $(Expr(:call, fun)))
    append!(expr.args, args)
    return expr
end

exprs = toexpr.([a * p1 + b * p2, sin(c * p1) + tanh(b * p3)],
    Ref(NameState(Dict(a => :a_, b => :b_, c => :c_,
        p1 => :p1_, p2 => :p2_, p3 => :p3_))))

f1_expr = build_f1(input_names, param_names, outputs, exprs)
bf1 = @RuntimeGeneratedFunction(f1_expr)

cc = toexpr(a + b)
input = ones(2, 10)
params = [2, 3, 4]
bf1(eachslice(input, dims=1), params)

# f2_expr = :(function (i, p)
#     i1 = i[1]
#     i2 = i[2]
#     p1, p2, p3 = p[1], p[2], p[3]
#     c = @. (sin(i1 * p1) + tanh(i2 * p2))
#     d = @. (c * p1 + i2 * p3)
#     [c, d]
# end)

# f3_expr = quote
#     function (i, p)
#         a_ = i[1]
#         b_ = i[2]
#         p1_ = p[1]
#         p2_ = p[2]
#         p3_ = p[3]
#         c_ = @__dot__((+)((*)(a_, p1_), (*)(b_, p2_)))
#         d_ = @__dot__((+)((tanh)((*)(b_, p3_)), (sin)((*)(c_, p1_))))
#         return [c_, d_]
#     end
# end

# bf2 = @RuntimeGeneratedFunction(f2_expr)
# bf3 = @RuntimeGeneratedFunction(f3_expr)
# input = ones(2, 10)
# params = [2, 3, 4]
# bf1(eachslice(input, dims=1), params)
# bf2(eachslice(input, dims=1), params)
# bf3(eachslice(input, dims=1), params)


# eval(:(d = @. ((+)((tanh)((*)(input1, 2)), (sin)((*)(input2, 3))))))
# input1 = eachslice(input, dims=1)[1]
# input2 = eachslice(input, dims=1)[2]
# # eval(:(@__dot__((+)((tanh)((*)(input1, 2)), (sin)((*)(input2, 3))))))

# # @btime Zygote.gradient((p) -> sum(reduce(hcat, f_1(input, p))), [2, 3, 4])
# # @btime Zygote.gradient((p) -> sum(reduce(hcat, bf1(eachslice(input, dims=1), p))), [2, 3, 4])
# # @btime Zygote.gradient((p) -> sum(reduce(hcat, bf1.(eachslice(input, dims=2), Ref(p)))), [2, 3, 4])
# # @btime Zygote.gradient((p) -> sum(reduce(hcat, cf1.(eachslice(input, dims=2), Ref(p)))), [2, 3, 4])
# # a_ = [1, 2, 3]
# # toexpr(tanh(a + b * c))

# # eval(:(@__dot__((+)((tanh)((*)([1,2,3], 1)), (sin)((*)([1,2,3], 2))))))

