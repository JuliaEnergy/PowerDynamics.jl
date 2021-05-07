using BlockSystems
using ModelingToolkit
using ModelingToolkit: getname

export MetaGenerator

"""
         i_r,i_i    u_r,u_i
            ↓          ↑
        +------------------+
δ(t) --→| rot(δ)   rot(-δ) |--→ |v_h|
        +------------------+
            ↓          ↑
         i_d,i_q    v_d,v_q

"""
function reference_frame_conversion()
    @parameters t i_r(t) i_i(t) v_d(t) v_q(t) δ(t)
    @variables u_r(t) u_i(t) i_d(t) i_q(t) v_h(t)

    IOBlock([u_r ~ sin(-δ)*v_d - cos(-δ)*v_q,
             u_i ~ cos(-δ)*v_d + sin(-δ)*v_q,
             i_d ~ sin( δ)*i_r - cos( δ)*i_i,
             i_q ~ cos( δ)*i_r + sin( δ)*i_i,
             v_h ~ √(v_d^2 + v_q^2)],
            [δ, i_r, i_i, v_d, v_q], [u_r, u_i, i_d, i_q, v_h], name=:rfc)
end

"""
    MetaGenerator(mover, shaft, machine, controller, parameters; name)

Implementation of the NREL MetaGenerator model. This function takes several
`IOBlocks`` as inputs which should satisfy the corresponding block specification (see below).

The `parameters` should be provided as a Dict of parameters as namespaced syms, i.e.
`mover.P_ref => 123`.

Block specifications: `(input, [optional]) ↦ (output)`
 - mover                      `([ω]) ↦ (τ_m)`
 - shaft                 `(τ_m, τ_e) ↦ (δ, ω)`
 - machine               `(i_d, i_q) ↦ (v_d, v_q, τ_e)`
 - controller      `([v_df, τ_e, ω]) ↦ (v_f)`

Symbols:
 - `u_r + j u_i`; `i_r + j i_i` : current and voltage in voltage reference system (θ)
 - `v_d + j v_q`; `i_d + j i_q` : current and voltage in shaft reference system (δ)
 - `v_h` : voltage amplitude
 - `ω`, `δ` : speed and relative position (to voltage angle θ) of rotating machine
 - `τ_m`, `τ_e` : mechanical and electric torque
 - `v_f` : field voltage

```
         i_r,i_i    u_r,u_i
            ↓          ↑
     +-----------------------+ v_h
 +--→| ref. frame conversion |--------------+
 |   +-----------------------+              |  (optional inputs)
 |  i_d,i_q |          ↑ v_d,v_q            | /
 |          ↓          |                   (↓)
 |   +-----------------------+    v_f +------------+
 |   |        machine        |(←)-----| Controller |
 |   +-----------------------+        +------------+
 |               | τ_e                   (↑)  (↑)
 |               +------------------------+    |
 |               ↓                             |
 | δ +-----------------------+ ω               |
 +---|         shaft         |----+------------+
     +-----------------------+    |
                 ↑ τ_m            |
                 |                |
     +-----------------------+    |
     |         mover         |(←)-+
     +-----------------------+
```


TODO: I went with the NREL convention v instead of u. Should be discussed.
"""
function MetaGenerator(mover, shaft, machine, controller, para;
                       verbose=false, name=gensym(:MGenerator))

    verbose && println("Create MetaGenerator:")
    rfc = reference_frame_conversion()

    @assert BlockSpec([], [:τ_m])(mover) "mover has not output τₘ!"
    @assert BlockSpec([:τ_m, :τ_e], [:δ, :ω])(shaft) "shaft not of type (τₑ, τₘ) ↦ (δ, ω)!"
    @assert BlockSpec([:i_d, :i_q], [:v_d, :v_q, :τ_e])(machine) "machine not of type (i_dq) ↦ (v_h, τₑ)!"
    @assert BlockSpec([], [:v_f])(controller) "controller has not output v_f!"

    # create all necessary connections
    connections = [mover.τ_m => shaft.τ_m,
                   shaft.δ => rfc.δ,
                   rfc.i_d => machine.i_d,
                   rfc.i_q => machine.i_q,
                   machine.v_d => rfc.v_d,
                   machine.v_q => rfc.v_q,
                   machine.τ_e => shaft.τ_e]

    # optional additional connections to machine
    if :v_f ∈ getname.(machine.inputs)
        push!(connections, controller.v_f => machine.v_f)
        verbose && println("  - added optional connection: controller.v_f => machine.v_f")
    end
    # optional additional connections to controller
    controller_inputs = getname.(controller.inputs)
    if :v_h ∈ controller_inputs
        push!(connections, rfc.v_h => controller.v_h)
        verbose && println("  - added optional connection: rfc.v_h => controller.v_h")
    end
    if :τ_e ∈ controller_inputs
        push!(connections, machine.τ_e => controller.τ_e)
        verbose && println("  - added optional connection: machine.τ_e => controller.τ_e")
    end
    if :ω ∈ controller_inputs
        push!(connections, shaft.ω => controller.ω)
        verbose && println("  - added optional connection: shaft.ω => controller.ω")
    end
    # optional additional connection to mover
    if :ω ∈ getname.(mover.inputs)
        push!(connections, shaft.ω => mover.ω)
        verbose && println("  - added optional connection: shaft.ω => mover.ω")
    end

    # make sure that v_r, u_i, i_r, i_i get promoted
    rfc_promotions = [rfc.u_r => :u_r, rfc.u_i => :u_i,
                      rfc.i_r => :i_r, rfc.i_i => :i_i]

    sys = IOSystem(connections, [rfc, mover, shaft, machine, controller],
                   namespace_map=rfc_promotions, outputs=[rfc.u_r, rfc.u_i],
                   autopromote=false, name=name)
    verbose && println("IOSystem: ", sys)
    # we want no further autopromote to keep all of the parameters in the namespaces
    # this is necessary because the parameterdict will be namespaced as well!
    connected = connect_system(sys, verbose=false)
    verbose && println("IOBlock: ", connected)

    inputs = Set(getname.(connected.inputs))

    @assert Set((:i_i, :i_r)) == inputs "System inputs [:i_i, :i_r] != $inputs. It seems like there are still open inputs in the subsystems!"
    @assert Set(connected.iparams) == Set(keys(para)) "parameter dict does not equal system parameters $(connected.iparams)"

    IONode(connected, para)
end
