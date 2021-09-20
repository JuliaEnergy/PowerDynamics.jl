export BusNode, BusLoad

BusNode(injs::BlockPara...; kwargs...) = BusNode(injs; kwargs...)

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

function BusLoad(;P, Q, name=:load)
    pv, qv = P, Q
    @parameters t P Q u_r(t) u_i(t)
    @variables i_r(t) i_i(t)
    block = IOBlock([0 ~ u_r*i_r + u_i*i_i - P,
                     0 ~ u_i*i_r - u_r*i_i - Q],
                    [u_r, u_i],
                    [i_r, i_i];
                    name, iv=t)
    para = Dict(block.P => pv, block.Q => qv)
    BlockPara(block, para)
end
