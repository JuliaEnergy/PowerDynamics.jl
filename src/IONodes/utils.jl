using BlockSystems
using BlockSystems: remove_namespace
using BlockSystems.ModelingToolkit: getname, Symbolic

export BusNode, BlockPara
export get_block, get_parameters

"""
    BlockPara(block, para; strict=true)

Composite type holdes an `IOBlock` and a parameter `Dict`.
Parameter `Dict` may be provided with Symbols (`:a`) or Symbolic types (`blk.a`).
The latter will be transformed to Symbols.

If `strict=true` the construtor asserts that the given parameters match the
internal parameters of the `IOBlock`.
"""
struct BlockPara
    block::IOBlock
    para::Dict{Symbol, Float64}
    function BlockPara(blk::IOBlock, para; strict=true)
        # only allow symbols as keys
        newp = Dict(to_symbol(blk, k) => v for (k, v) in para)

        allp = Set(getname.(blk.iparams))
        givenp = keys(newp)
        additionalp = collect(setdiff(givenp, allp))
        if !isempty(additionalp)
            throw(ArgumentError("Got unknown parameters $(additionalp)!"))
        end
        if strict && !(allp == givenp)
            missingp = collect(setdiff(allp, givenp))
            throw(ArgumentError("Parameters $(missingp) missing!"))
        end
        new(blk, newp)
    end
end

get_block(bp::BlockPara) = bp.block
get_parameters(bp::BlockPara) = bp.para

"""
    to_symbol(b::AbstractIOSystem, x)

Provide `x` as Symbol or Symbolic (with or without namespace) and get a Symbol without the namespace back.
"""
to_symbol(b::AbstractIOSystem, x::Symbol)   = x
to_symbol(b::AbstractIOSystem, x::Num)      = to_symbol(b, ModelingToolkit.value(x))
to_symbol(b::AbstractIOSystem, x::Symbolic) = remove_namespace(b.name, getname(x))

"""
    mergep(bps::BlockPara...; strict=true)

Returns tuple of `IOBlock`s ands a merged dict of all parameters.
All parameter names will get namespaced.
"""
mergep(bps::BlockPara...; kwargs...) = mergep(bps; kwargs...)
function mergep(bps; strict=true)
    isempty(bps) && return (NTuple{0, IOBlock}(), Dict{Symbol,Float64}())

    pdicts = [Dict(nspc(bp.block, k)=>v for (k,v) in bp.para)
         for bp in bps]

    return get_block.(bps), merge(pdicts...)
end

nspc(blk::IOBlock, s::Symbol) = Symbol(String(blk.name) * "₊" * String(s))

"""
    subscript(s, i)

Append symbol or string `s` with a integer subscript.
"""
function subscript(s::T, i::Int) where T <: Union{Symbol, AbstractString}
    @assert 0≤i≤9 "subscript only supported from 0..9"
    Symbol(s, Char(0x02080 + i))
end
