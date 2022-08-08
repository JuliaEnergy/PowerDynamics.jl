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
                  iv=t, name, warn=false)

    return IONode(bar, Dict{Symbol, Float64}())
end

function BusNode(injpairs; name=gensym(:Bus), verbose=false)
    injectors, params = mergep(injpairs)

    for inj in injectors
        @assert BlockSpec([:u_r, :u_i], [:i_r, :i_i])(inj) "Injector $inj does not satisfy injector interface!"
    end

    # TODO assert that all iv ar equal, prob done in BlockSystems?
    t = get_iv(injectors[1])

    ir, ii = Num[], Num[]
    for i in 1:length(injpairs)
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
                         iv=t, warn=false)

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
