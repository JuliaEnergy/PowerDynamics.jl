"""
    get_pd_branch(arc<:branch, branches::Vector{Branch})

Add methods for each supported subtype of `Branch`. Returns an instance
of `AbstractLine`. There are potentially more than one branch objects
"""
function get_pd_branch end

inorder(arc::Arc) = get_name(arc.from) < get_name(arc.to)

function get_pd_branch(arc::Arc, lines::AbstractVector{Line})
    src = get_name(arc.from)
    dst = get_name(arc.to)

    y = zero(ComplexF64)
    shunt_src_X = zero(Float64)
    shunt_dst_X = zero(Float64)

    for line in lines
        y += 1/(line.r + im*line.x)
        # minus becaus of 1/im
        shunt_src_X -= 1/line.b.from
        shunt_dst_X -= 1/line.b.to
    end


    if inorder(arc)
        return PiModelLine(; from=src, to=dst, y, y_shunt_km=-1/shunt_src_X, y_shunt_mk=-1/shunt_dst_X)
    else
        return PiModelLine(; from=dst, to=srf, y, y_shunt_km=-1/shunt_dst_X, y_shunt_mk=-1/shunt_src_X)
    end
end

function get_pd_branch(arc::Arc,
                       trafos::AbstractVector{<:Union{TapTransformer, Transformer2W}})
    src = get_name(arc.from)
    dst = get_name(arc.to)

    length(trafos)==1 || error("Multiple transformers at one Arc are not supportet atm 😔")
    trafo = trafos[1]

    y = 1/(trafo.r + im*trafo.x)

    if inorder(arc)
        return Transformer(; from=src, to=dst, y, t_ratio=trafo.rate)
    else
        return Transformer(; from=dst, to=src, y, t_ratio=1/trafo.rate)
    end
end
