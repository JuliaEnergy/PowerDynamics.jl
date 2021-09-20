"""
    gen_mover_block(machine <: PowerSystems.TurbineGov)

Returns composite type `BlockPara` which consists of an `IOBlock` and a parameter `Dict`.
"""
function gen_mover_block end


"""
    function tg_fixed(;name=:machine)

Returns an IOBlock which represents the `TGFixed` from `PowerSystems.jl`.

See [`PowerSimulationDynamics.jl` docs](https://nrel-siip.github.io/PowerSimulationsDynamics.jl/stable/component_models/turbine_gov/#Fixed-TG-[TGFixed])

```jldoctest
julia> PowerDynamics.tg_fixed()
IOBlock :mover with 1 eqs
  ├ inputs:  (empty)
  ├ outputs: τ_m(t)
  ├ istates: (empty)
  └ iparams: P_ref, η
```
"""
function tg_fixed(;name=:mover)
    @parameters t P_ref η
    @variables τ_m(t)

    IOBlock([τ_m ~ η * P_ref],
            [], [τ_m], name=name)
end

function gen_mover_block(mover::TGFixed)
    block = tg_fixed()
    p = Dict(block.η => mover.efficiency,
             block.P_ref => mover.P_ref)
    return BlockPara(block, p)
end
