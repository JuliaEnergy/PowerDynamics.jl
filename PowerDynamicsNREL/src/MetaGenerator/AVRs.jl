"""
    gen_avr_block(avr <: PowerSystems.AVR)

Returns composite type `BlockPara` which consists of an `IOBlock` and a parameter `Dict`.
"""
function gen_avr_block end


"""
    function avr_fixed(;name=:AVR)

Returns an IOBlock which represents the `AVRFixed` from `PowerSystems.jl`.

See [`PowerSimulationDynamics.jl` docs](https://nrel-siip.github.io/PowerSimulationsDynamics.jl/stable/component_models/avr/#Fixed-AVR-[AVRFixed])

```jldoctest
julia> PowerDynamics.avr_fixed()
IOBlock :AVR with 1 eqs
  ├ inputs:  (empty)
  ├ outputs: v_avr(t)
  ├ istates: (empty)
  └ iparams: v_fix
```
"""
function avr_fixed(;name=:AVR)
    @parameters t v_fix
    @variables v_avr(t)

    IOBlock([v_avr ~ v_fix],
            [], [v_avr], name=name)
end

function gen_avr_block(avr::AVRFixed)
    block = avr_fixed()
    @warn "Use Vf as v_fix. Is this right? May be v_ref!"
    p = Dict(block.v_fix => avr.Vf) #FIXME: Vf or V_ref for AVRFixed
    return BlockPara(block, p)
end


"""
    function avr_simple(;name=:AVR)

Returns an IOBlock which represents the `AVRSimple` from `PowerSystems.jl`.

See [`PowerSimulationDynamics.jl` docs](https://nrel-siip.github.io/PowerSimulationsDynamics.jl/stable/component_models/avr/#Simple-AVR-[AVRSimple])

```jldoctest
julia> PowerDynamics.avr_simple()
IOBlock :AVR with 1 eqs
  ├ inputs:  v_h(t)
  ├ outputs: v_avr(t)
  ├ istates: (empty)
  └ iparams: K, v_ref
```
"""
function avr_simple(;name=:AVR)
    @parameters t K v_h(t) v_ref
    @variables v_avr(t)
    dt = Differential(t)

    IOBlock([dt(v_avr) ~ K*(v_ref - v_h)],
            [v_h], [v_avr], name=name)
end

function gen_avr_block(avr::AVRSimple)
    block = avr_simple()
    p = Dict(block.K => avr.Kv,
             block.v_ref => avr.V_ref)
    return BlockPara(block, p)
end
