"""
    gen_avr_block(avr <: PowerSystems.AVR)

Returns tuple of
 - `IOBlock` representaion of `avr`
 - parameter dict for this block based on the contents of `avr`
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
  ├ outputs: v_f(t)
  ├ istates: (empty)
  └ iparams: v_fix
```
"""
function avr_fixed(;name=:AVR)
    @parameters t v_fix
    @variables v_f(t)

    IOBlock([v_f ~ v_fix],
            [], [v_f], name=name)
end

function gen_avr_block(avr::AVRFixed)
    block = avr_fixed()
    p = Dict(block.v_fix => avr.Vf) #FIXME: Vf or V_ref for AVRFixed
    return (block, p)
end


"""
    function avr_simple(;name=:AVR)

Returns an IOBlock which represents the `AVRSimple` from `PowerSystems.jl`.

See [`PowerSimulationDynamics.jl` docs](https://nrel-siip.github.io/PowerSimulationsDynamics.jl/stable/component_models/avr/#Simple-AVR-[AVRSimple])

```jldoctest
julia> PowerDynamics.avr_simple()
IOBlock :AVR with 1 eqs
  ├ inputs:  v_h(t)
  ├ outputs: v_f(t)
  ├ istates: (empty)
  └ iparams: K, v_ref
```
"""
function avr_simple(;name=:AVR)
    @parameters t K v_h(t) v_ref
    @variables v_f(t)
    dt = Differential(t)

    IOBlock([dt(v_f) ~ K*(v_ref - v_h)],
            [v_h], [v_f], name=name)
end

function gen_avr_block(avr::AVRSimple)
    block = avr_simple()
    p = Dict(block.K => avr.Kv,
             block.v_ref => avr.V_ref)
    return (block, p)
end
