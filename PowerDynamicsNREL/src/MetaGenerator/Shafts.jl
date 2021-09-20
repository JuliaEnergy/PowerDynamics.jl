"""
    gen_shaft_block(shaft <: PowerSystems.Shaft)

Returns tuple of
 - `IOBlock` representaion of `shaft`
 - parameter dict for this block based on the contents of `shaft`
"""
function gen_shaft_block end


"""
    function single_mass_shaft(;name=:shaft)

Returns an IOBlock which represents the `SingleMassShaft` from `PowerSystems.jl`.

See [`PowerSimulationDynamics.jl` docs](https://nrel-siip.github.io/PowerSimulationsDynamics.jl/stable/component_models/shafts/#Rotor-Mass-Shaft-[SingleMass])

```jldoctest
julia> PowerDynamics.single_mass_shaft()
IOBlock :shaft with 2 eqs
  ├ inputs:  τ_m(t), τ_e(t)
  ├ outputs: δ(t), ω(t)
  ├ istates: (empty)
  └ iparams: Ω_b, ω_ref, H, D
```
"""
function single_mass_shaft(;name=:shaft)
    @parameters t τ_m(t) τ_e(t) D H ω_ref Ω_b
    @variables δ(t) ω(t)
    dt = Differential(t)

    IOBlock([dt(δ) ~ Ω_b * (ω - ω_ref),
             dt(ω) ~ 1/(2H) * (τ_m - τ_e - D*(ω - ω_ref))],
            [τ_m, τ_e], [δ, ω], name=name)
end

function gen_shaft_block(s::SingleMass)
    shaft = single_mass_shaft()
    p = Dict(shaft.D => s.D,
             shaft.H => s.H,
             shaft.ω_ref => 1.0,
             shaft.Ω_b => 2π*60.0)
    # FIXME where to geht ω_ref and Ω_b? -> constants in MetaGenerator
    return (shaft, p)
end
