"""
    UHFunction{uhtype} <: Function

Represents a unit hydrograph function for routing water through a hydrological system.

# Fields
- `uhtype::Symbol`: A symbol indicating the type of unit hydrograph function. Currently, only `:UH_1_HALF` and `:UH_2_FULL` are supported.

# Methods
- `(uh::UHFunction{uhtype})(t, lag)`: Computes the unit hydrograph value at time `t` given the lag time `lag`.
- `get_uh_tmax(uh::UHFunction{uhtype}, lag)`: Returns the maximum time required for computing the unit hydrograph with the given lag time `lag`.

"""
struct UHFunction{uhtype}
    function UHFunction(uhtype::Symbol)
        return new{uhtype}()
    end
end

function (uh::UHFunction{:UH_1_HALF})(t, lag)
    if t - lag > 0
        typeof(lag)(1)
    else
        (t / lag)^2.5
    end
end

get_uh_tmax(::UHFunction{:UH_1_HALF}, lag) = ceil(lag)

function (uh::UHFunction{:UH_2_FULL})(t, lag)
    if t - lag * 2 > 0
        typeof(lag)(1)
    elseif t - lag > 0
        (1 - 0.5 * abs(2 - t / lag)^2.5)
    else
        (0.5 * abs(t / lag)^2.5)
    end
end

get_uh_tmax(::UHFunction{:UH_2_FULL}, lag) = 2 * ceil(lag)


function uh_3_half(lag; kw...)
    sf = get(kw, :smooth_func, ifelse_func)

    timeidx = 1:ceil(lag)
    ff = 1 / (0.5 * lag^2)
    value = @.(sf(lag - timeidx) * ff * (0.5 * timeidx^2 - 0.5 * (timeidx - 1)^2) +
               sf(timeidx - lag) * ff * (0.5 * lag^2 - 0.5 * (timeidx - 1)^2))
    return value
end

