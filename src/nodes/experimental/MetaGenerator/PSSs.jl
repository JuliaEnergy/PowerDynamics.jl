"""
    gen_pss_block(pss <: PowerSystems.PSS)

Returns tuple of
 - `IOBlock` representaion of `pss`
 - parameter dict for this block based on the contents of `pss`
"""
function gen_pss_block end


"""
    function pss_fixed(;name=:PSS)

Returns an IOBlock which represents the `PSSFixed` from `PowerSystems.jl`.

See [`PowerSimulationDynamics.jl` docs](https://nrel-siip.github.io/PowerSimulationsDynamics.jl/stable/component_models/pss/#Fixed-PSS-[PSSFixed])

```jldoctest
julia> PowerDynamics.pss_fixed()
IOBlock :PSS with 1 eqs
  ├ inputs:  (empty)
  ├ outputs: v_pss(t)
  ├ istates: (empty)
  └ iparams: v_fix
```
"""
function pss_fixed(;name=:PSS)
    @parameters t v_fix
    @variables v_pss(t)

    IOBlock([v_pss ~ v_fix],
            [], [v_pss], name=name)
end

function gen_pss_block(pss::PSSFixed)
    block = pss_fixed()
    p = Dict(block.v_fix => pss.V_pss)
    return (block, p)
end


"""
    function pss_simple(;name=:PSS)

Returns an IOBlock which represents the `PSSSimple` from `PowerSystems.jl`.

See [`PowerSimulationDynamics.jl` docs](https://nrel-siip.github.io/PowerSimulationsDynamics.jl/stable/component_models/pss/#Simple-PSS-[PSSSimple])

```jldoctest
julia> PowerDynamics.pss_simple()
IOBlock :PSS with 1 eqs
  ├ inputs:  ω(t), ω_ref(t), τ_e(t), P_ref(t)
  ├ outputs: v_pss(t)
  ├ istates: (empty)
  └ iparams: K_p, K_ω
```
"""
function pss_simple(;name=:PSS)
    @parameters t K_ω K_p
    @parameters ω(t) ω_ref(t) P_ref(t) τ_e(t)
    @variables v_pss(t)

    IOBlock([v_pss ~ K_ω*(ω - ω_ref) + K_p*(ω*τ_e - P_ref)],
            [ω, ω_ref, τ_e, P_ref], [v_pss], name=name)
end

function gen_pss_block(pss::PSSSimple)
    block = pss_simple()
    p = Dict(block.K_ω => pss.K_ω,
             block.K_p => pss.K_p)
    return (block, p)
end
