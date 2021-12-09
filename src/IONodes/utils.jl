using BlockSystems

export BusNode, BlockPara
export get_block, get_parameters

"""
    struct BlockPara{P<:Dict}
    BlockPara(block, para; strict=true)

Composite type holdes an `IOBlock` and a parameter `Dict`.
If `strict=true` the construtor asserts that all internal parameters
of the `IOBlock` are given in the parameter dict.
"""
struct BlockPara{P<:Dict}
    block::IOBlock
    para::P
    function BlockPara(block::IOBlock, para::P; strict=true) where {P<:Dict}
        if strict && !(Set(BlockSystems.namespace_iparams(block)) ⊆ keys(para))
            throw(ArgumentError("There is a missmatch between given parameter dict and iparams: $(BlockSystems.namespace_iparams(block)), $(keys(para))"))
        end
        new{P}(block, para)
    end
end

get_block(bp::BlockPara) = bp.block
get_parameters(bp::BlockPara) = bp.para


"""
    mergep(bps::BlockPara...; strict=true)

Returns tuple of `IOBlock`s ands a merged dict of all parameters.
If `strict=true` check whether the keys in the parameterdicts are unique.
"""
mergep(bps::BlockPara...; kwargs...) = mergep(bps; kwargs...)
function mergep(bps::NTuple{N, BlockPara}; strict=true) where {N}
    isempty(bps) && return (NTuple{0, IOBlock}(), Dict())

    paras = get_parameters.(bps)

    if strict && !allunique(vcat(collect.(keys.(paras))...))
        throw(ArgumentError("Keys of parameter dicts not unique: $(vcat(collect.(keys.(paras))...))"))
    end

    return get_block.(bps), merge(paras...)
end


"""
    subscript(s, i)

Append symbol or string `s` with a integer subscript.
"""
function subscript(s::T, i::Int) where T <: Union{Symbol, AbstractString}
    @assert 0≤i≤9 "subscript only supported from 0..9"
    Symbol(s, Char(0x02080 + i))
end
