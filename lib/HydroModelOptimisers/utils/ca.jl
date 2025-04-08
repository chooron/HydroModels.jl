# merge ComponentArray:https://github.com/jonniedie/ComponentArrays.jl/issues/186
#* componentarray update
function update_ca(ca::ComponentArray{T1}, ca2::ComponentArray{T2}) where {T1,T2}
    ax = getaxes(ca)
    ax2 = getaxes(ca2)
    vks = valkeys(ax[1])
    vks2 = valkeys(ax2[1])
    _p = Vector{T2}()
    for vk in vks
        if length(getaxes(ca[vk])) > 0
            _p = vcat(_p, collect(update_ca(ca[vk], vk in vks2 ? getproperty(ca2, vk) : ComponentVector())))
        else
            if vk in vks2
                _p = vcat(_p, ca2[vk])
            else
                _p = vcat(_p, ca[vk])
            end
        end
    end
    ComponentArray(_p, ax)
end

#* merge componentarray
function merge_ca(ca1::ComponentVector{T1}, ca2::ComponentVector{T2}) where {T1,T2}
    ax = getaxes(ca1)
    ax2 = getaxes(ca2)
    vks = valkeys(ax[1])
    vks2 = valkeys(ax2[1])
    idxmap = indexmap(ax[1])
    _p = Vector{promote_type(T1,T2)}()
    sizehint!(_p, length(ca1)+length(ca2))
    for vk in vks
        if vk in vks2
            _p = vcat(_p, ca2[vk])
        else
            _p = vcat(_p, ca1[vk])
        end
    end
    new_idxmap = Vector{Pair{Symbol, Int64}}([])
    sizehint!(new_idxmap, length(ca2))
    max_val = maximum(idxmap)
    for vk in vks2
        if !(vk in vks)
            _p = vcat(_p, ca2[vk])
            new_idxmap = vcat(new_idxmap, [getval(vk)=>max_val+1])
            max_val += 1
        end
    end
    merged_ax = Axis(merge(idxmap, new_idxmap))
    ComponentArray(_p, merged_ax)
end