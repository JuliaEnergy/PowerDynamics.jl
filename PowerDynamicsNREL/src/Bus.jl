export BusNode, BusLoad

BusNode(injs::BlockPara...; kwargs...) = BusNode(injs; kwargs...)

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

function BusLoad(;P, Q, name=:load)
    pv, qv = P, Q
    @parameters t P Q u_r(t) u_i(t)
    @variables i_r(t) i_i(t)
    block = IOBlock([0 ~ u_r*i_r + u_i*i_i + P,
                     0 ~ u_i*i_r - u_r*i_i + Q],
                    [u_r, u_i],
                    [i_r, i_i];
                    name, iv=t)
    para = Dict(block.P => pv, block.Q => qv)
    BlockPara(block, para)
end


"""
    get_io_load(load<:ElectricLoad)

Add methods for each supported subtype of `ElectricLoad`. Returns
an `BlockPara` for the current draw/injection.
"""
function get_io_load end

function get_io_load(load::PowerLoad)
    @assert load.dynamic_injector == nothing "Whait whaat load $(load.name) has a dynamic injection: $(load.dynamic_injector)"
    BusLoad(P=load.active_power, Q=load.reactive_power, name=Symbol(get_name(load)))
end


"""
    get_io_injection(inj<:DynamicInjector)

Add methods for each supported subtype of `DynamicInjector`. Returns
an `BlockPara` for the current draw/injection.
"""
function get_io_injection end

function get_io_injection(gen::DynamicGenerator; verbose=false)
    mover = gen_mover_block(gen.prime_mover)
    verbose && @info "mover loaded:" mover
    shaft = gen_shaft_block(gen.shaft)
    verbose && @info "shaft loaded:" shaft
    machine = gen_machine_block(gen.machine)
    verbose && @info "machine loaded:" machine
    avr = gen_avr_block(gen.avr)
    verbose && @info "avr loaded:" avr
    pss = gen_pss_block(gen.pss)
    verbose && @info "pss loaded:" pss

    MetaGenerator(mover, shaft, machine, avr, pss;
                  verbose, name=Symbol(gen.name))
end
