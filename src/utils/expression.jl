# git from https://github.com/JuliaSymbolics/SymbolicUtils.jl/blob/master/src/code.jl

"""
Convert a symbolic expression to an expression, special handling of broadcast.

# Arguments
- `x`: The symbolic expression to convert.

# Returns
- `Expr`: The converted expression.
"""
toexprv2(x) = toexprv2(x, SymbolicUtils.Code.LazyState())

function function_to_exprv2(op, O, st)
    (get(st.rewrites, :nanmath, false) && op in SymbolicUtils.Code.NaNMathFuns) || return nothing
    name = nameof(op)
    fun = GlobalRef(SymbolicUtils.Code.NaNMath, name)
    args = map(Base.Fix2(toexprv2, st), arguments(O))
    expr = Expr(:call, fun)
    append!(expr.args, args)
    return expr
end

function function_to_exprv2(op::Union{typeof(*),typeof(+)}, O, st)
    out = get(st.rewrites, O, nothing)
    out === nothing || return out
    args = map(Base.Fix2(toexprv2, st), sorted_arguments(O))
    if length(args) >= 3 && symtype(O) <: Number
        x, xs = Iterators.peel(args)
        foldl(xs, init=x) do a, b
            Expr(:call, :broadcast, op, a, b) # add broadcast
        end
    else
        expr = Expr(:call, :broadcast, op) # add broadcast
        append!(expr.args, args)
        expr
    end
end

function function_to_exprv2(op::typeof(^), O, st)
    args = arguments(O)
    if args[2] isa Real && args[2] < 0
        args[1] = Term(inv, Any[args[1]])
        args[2] = -args[2]
    end
    if isequal(args[2], 1)
        return toexprv2(args[1], st)
    end
    if get(st.rewrites, :nanmath, false) === true && !(args[2] isa Integer)
        op = SymbolicUtils.Code.NaNMath.pow
        return toexprv2(Term(op, args), st)
    end
    return nothing
end

function function_to_exprv2(::typeof(ifelse), O, st)
    args = arguments(O)
    :($(toexprv2(args[1], st)) ? $(toexprv2(args[2], st)) : $(toexprv2(args[3], st)))
end

function function_to_exprv2(x::BasicSymbolic, O, st)
    issym(x) ? get(st.rewrites, O, nothing) : nothing
end

toexprv2(O::Expr, st) = O

function substitute_name(O, st)
    if (issym(O) || iscall(O)) && haskey(st.rewrites, O)
        st.rewrites[O]
    else
        O
    end
end

function toexprv2(O, st)
    if issym(O)
        O = substitute_name(O, st)
        return issym(O) ? nameof(O) : toexprv2(O, st)
    end
    O = substitute_name(O, st)

    !iscall(O) && return O
    op = operation(O)
    expr′ = function_to_exprv2(op, O, st)
    if expr′ !== nothing
        return expr′
    else
        !iscall(O) && return O
        args = arguments(O)
        return Expr(:call, :broadcast, toexprv2(op, st), map(x -> toexprv2(x, st), args)...)
    end
end

export toexprv2
