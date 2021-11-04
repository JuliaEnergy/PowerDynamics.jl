using ModelingToolkit

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
        if strict && !(Set(namespace_iparams(block)) ⊆ keys(para))
            throw(ArgumentError("There is a missmatch between given parameter dict and iparams: $(namespace_iparams(block)), $(keys(para))"))
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


"""
    BusNode(inj::BlockPara...; name=gensym(:Bus), verbose=false)

Create an [`IONode`](@ref) based on several current injectors/draws.
Each current injector has to be provided as a `BlockPara` object, a composite of an
`IOBlock` and a parameter dict.
Each `IOBlock` has to fulfil the interface `u_r, u_i ↦ i_r, i_i`.
"""
BusNode(injs::BlockPara...; kwargs...) = BusNode(injs; kwargs...)

# special function definition for empty bus
function BusNode(; name=gensym(:Bus), verbose=false)
    # busnode without any elements
    @parameters t i_r(t) i_i(t)
    @variables u_r(t) u_i(t)

    bar = IOBlock([0 ~ i_r, 0 ~ i_i],
                  [i_r, i_i], [u_r, u_i];
                  iv=t, name)

    return IONode(bar, Dict{Symbolic, Float64}())
end

function BusNode(injpairs::NTuple{N, BlockPara}; name=gensym(:Bus), verbose=false) where {N}
    injectors, params = mergep(injpairs)

    for inj in injectors
        @assert BlockSpec([:u_r, :u_i], [:i_r, :i_i])(inj) "Injector $inj does not satisfy injector interface!"
    end

    # TODO assert that all iv ar equal, prob done in BlockSystems?
    t = get_iv(injectors[1])

    ir, ii = Num[], Num[]
    for i in 1:N
        irs = subscript(:i_r, i)
        iis = subscript(:i_i, i)
        append!(ir, @parameters $irs(t))
        append!(ii, @parameters $iis(t))
    end

    @parameters i_r(t) i_i(t)
    @variables u_r(t) u_i(t)

    @named bar = IOBlock([0 ~ (+)(ir...) + i_r,
                          0 ~ (+)(ii...) + i_i],
                         [i_r, i_i, ir..., ii...],
                         [u_r, u_i];
                         iv=t)

    connections = Pair[]
    for (i, inj) in enumerate(injectors)
        push!(connections, inj.i_r => getproperty(bar, subscript(:i_r, i)))
        push!(connections, inj.i_i => getproperty(bar, subscript(:i_i, i)))
        push!(connections, bar.u_r => inj.u_r)
        push!(connections, bar.u_i => inj.u_i)
    end

    promotions = [bar.u_r => :u_r,
                  bar.u_i => :u_i,
                  bar.i_r => :i_r,
                  bar.i_i => :i_i]

    sys = IOSystem(connections, [bar, injectors...];
                   namespace_map=promotions, autopromote=false,
                   outputs=[bar.u_r, bar.u_i],
                   name)

    connected = connect_system(sys, verbose=verbose)

    return IONode(connected, params)
end
