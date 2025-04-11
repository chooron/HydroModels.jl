using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)


func1_expr = quote
    function (i::string)
        return i
    end,
    function (i::int)
        return i + 1
    end
end

func1 = @RuntimeGeneratedFunction(func1_expr)