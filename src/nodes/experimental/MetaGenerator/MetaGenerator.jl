using BlockSystems
using ModelingToolkit
using ModelingToolkit: getname
using ModelingToolkit.SymbolicUtils: Symbolic
using PowerSystems

export MetaGenerator

# include the librarys for the components
include("Movers.jl")
include("Shafts.jl")
include("Machines.jl")
include("AVRs.jl")
include("PSSs.jl")

"""
    reference_frame_conversion()

The RFC is an IOBlock to transform between the different
reference frames.
The grid "provides" the current `i_r` and `i_i` in global reference frame. This
gets converted to `i_d` and `i_q` in machine reference frame.

On its way back the calculated voltages `v_d` and `v_q` in device frame get converted
back to the global reference fram `u_r` and `u_i`.

```
         i_r,i_i    u_r,u_i
            ↓          ↑
        +------------------+
δ(t) --→| rot(δ)   rot(-δ) |--→ |v_h|
        +------------------+
            ↓          ↑
         i_d,i_q    v_d,v_q
```
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
    MetaGenerator(mover, shaft, machine, AVR, parameters; name)

Implementation of the NREL MetaGenerator model. This function takes several
`IOBlocks`` as inputs which should satisfy the corresponding block specification (see below).

The `parameters` should be provided as a Dict of parameters as namespaced syms, i.e.
`mover.P_ref => 123`.

Block specifications: `(input) ↦ (output)`
 - mover                   `() ↦ (τ_m)`
 - shaft           `(τ_m, τ_e) ↦ (δ, ω)`
 - machine         `(i_d, i_q) ↦ (v_d, v_q, τ_e)`
 - AVR                     `() ↦ (v_avr)`
 - PSS                     `() ↦ (v_pss)`

All of the blocks may specify any of optional inputs `v_h`, `v_f`, `ω`, ...
to get access to the output of another blocks.

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
 +--→| ref. frame conversion |----→(all)   (all optional)
 |   +-----------------------+                   ↓
 |  i_d,i_q |          ↑ v_d,v_q      v_pss +---------+
 |          ↓          |               +----|   PSS   |
 |   +-----------------------+    v_f  ↓    +---------+
 |   |        machine        |(←)----(add)
 |   +-----------------------+         ↑    +---------+
 |               | τ_e                 +----|   AVR   |
 +-→(all)        +-----→(all)         v_avr +---------+
 |               ↓                               ↑
 | δ +-----------------------+ ω           (all optional)
 +---|         shaft         |--→(all)
     +-----------------------+
                 ↑
                 +-----→(all)
                 | τ_m
     +-----------------------+
     |         mover         |
     +-----------------------+
```


TODO: I went with the NREL convention v instead of u. Should be discussed.
"""
function MetaGenerator(mover, shaft, machine, AVR, PSS, para;
                       verbose=false, name=gensym(:MGenerator))
    para = copy(para) # create new reference

    verbose && println("Create MetaGenerator:")
    rfc = reference_frame_conversion()
    adder = IOComponents.Adder(name=:exitation, out=:v_f, a₁=:v_avr, a₂=:v_pss)

    @assert BlockSpec([], [:τ_m])(mover) "mover has not output τₘ!"
    @assert BlockSpec([:τ_m, :τ_e], [:δ, :ω])(shaft) "shaft not of type (τₑ, τₘ) ↦ (δ, ω)!"
    @assert BlockSpec([:i_d, :i_q], [:v_d, :v_q, :τ_e])(machine) "machine not of type (i_dq) ↦ (v_h, τₑ)!"
    @assert BlockSpec([], [:v_avr])(AVR) "AVR has not output v_avr!"
    @assert BlockSpec([], [:v_pss])(PSS) "PSS has not output v_pss!"

    # for more flexibility everything with matching name gets connected
    components = [mover, shaft, machine, AVR, PSS, rfc, adder]
    # create a dict of :sym => block
    # of each available output and the block where to find it as an output
    outputs = collect(sym => comp for comp in components for sym in getname.(comp.outputs))
    # all of the outputs need to have unique names, we can promote their namespaces!
    output_promotions = map(p -> getproperty(p.second, p.first) => p.first, outputs)
    verbose && println("Available Outputs: ", first.(output_promotions))
    @assert allunique(first.(outputs)) "There are multiple components which have outputs with the same name $(first.(output_promotions))"
    outputs = Dict(outputs...)

    # if any component specifies one of the available outputs as an input, connect to it!
    connections = Pair{Symbolic, Symbolic}[]
    for comp in components
        inputs = getname.(comp.inputs) # get the inputs as :sym without namespace
        for key in keys(outputs) ∩ inputs
            src = getproperty(outputs[key], key)
            dst = getproperty(comp, key)
            push!(connections, src => dst)
            verbose && println("  - added connection: $(src) => $(dst)")
        end
    end

    # make sure that the rfc inputs get promoted
    # the rfc outputs are allready in output_promotions
    rfc_promotions = [rfc.i_r => :i_r, rfc.i_i => :i_i]
    promotions = vcat(rfc_promotions, output_promotions)
    # we want no further autopromote to keep all of the parameters in the namespaces
    # this is necessary because the parameterdict will be namespaced as well!

    sys = IOSystem(connections, [rfc, mover, shaft, machine, AVR, PSS, adder],
                   namespace_map=promotions, outputs=[rfc.u_r, rfc.u_i],
                   autopromote=false, name=name)
    verbose && println("IOSystem before connection: ", sys)
    connected = connect_system(sys, verbose=verbose)
    verbose && println("IOBlock after connection: ", connected)

    inputs = Set(getname.(connected.inputs))
    @assert inputs == Set((:i_i, :i_r)) "System inputs [:i_i, :i_r] != $inputs. It seems like there are still open inputs in the subsystems!"

    # compare the provied parameters with the needed parameters
    needed_p = Set(connected.iparams)
    provided_p = keys(para)

    missing_p = Tuple(setdiff(needed_p, provided_p))
    isempty(missing_p) || throw(ArgumentError("Parameters $missing_p are missing from the parameterset!"))

    additional_p = Tuple(setdiff(provided_p, needed_p))
    isempty(additional_p) || @warn "Provided parameters $additional_p are not used in the equations and will be dropped!"
    for k in additional_p
        delete!(para, k)
    end

    IONode(connected, para)
end

function MetaGenerator(gen::DynamicGenerator; verbose=false)
    (mover, mover_p)     = gen_mover_block(gen.prime_mover)
    verbose && @info "mover loaded:" mover mover_p
    (shaft, shaft_p)     = gen_shaft_block(gen.shaft)
    verbose && @info "shaft loaded:" shaft shaft_p
    (machine, machine_p) = gen_machine_block(gen.machine)
    verbose && @info "machine loaded:" machine machine_p
    (avr, avr_p)         = gen_avr_block(gen.avr)
    verbose && @info "avr loaded:" avr avr_p
    (pss, pss_p)         = gen_pss_block(gen.pss)
    verbose && @info "pss loaded:" pss pss_p

    para = merge(mover_p, shaft_p, machine_p, avr_p, pss_p)

    MetaGenerator(mover, shaft, machine, avr, pss, para; verbose)
end
