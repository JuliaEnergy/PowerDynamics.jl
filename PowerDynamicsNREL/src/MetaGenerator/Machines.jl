"""
    gen_machine_block(machine <: PowerSystems.Machine)

Returns tuple of
 - `IOBlock` representaion of `machine`
 - parameter dict for this block based on the contents of `machine`
"""
function gen_machine_block end


"""
    function base_machine(;name=:machine)

Returns an IOBlock which represents the `BaseMachine` from `PowerSystems.jl`.

See [`PowerSimulationDynamics.jl` docs](https://nrel-siip.github.io/PowerSimulationsDynamics.jl/stable/component_models/machines/#Classical-Model-(Zero-Order)-[BaseMachine])

```jldoctest
julia> PowerDynamics.base_machine()
IOBlock :machine with 3 eqs
  ├ inputs:  v_d(t), v_q(t)
  ├ outputs: i_d(t), i_q(t), τ_e(t)
  ├ istates: (empty)
  └ iparams: X_d, R, e_q
```
"""
function base_machine(;name=:machine)
    @parameters t R X_d e_q v_d(t) v_q(t)
    @variables i_d(t) i_q(t) τ_e(t)

    IOBlock([i_d ~ (- R*v_d + e_q*X_d - v_q*X_d) / (R^2 + X_d^2),
             i_q ~ (  R*e_q - R*v_q   + v_d*X_d) / (R^2 + X_d^2),
             τ_e ~ (v_q + R*i_q)*i_q + (v_d + R*i_d)*i_d],
            [v_d, v_q], [i_d, i_q, τ_e], name=name)
end

function gen_machine_block(machine::BaseMachine)
    block = base_machine()
    p = Dict(block.R => machine.R,
             block.X_d => machine.Xd_p,
             block.e_q => machine.eq_p)
    return (block, p)
end
