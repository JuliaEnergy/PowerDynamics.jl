export MetaGenerator

"""
    reference_frame_conversion()

The RFC is an IOBlock to transform between the different
reference frames.
The grid "provides" the voltage `u_r` and `u_i` in global reference frame. This
gets converted to `v_d` and `v_q` in machine reference frame.

On its way back the calculated currents `i_d` and `i_q` in device frame get converted
back to the global reference fram `i_r` and `i_i`.

```
         u_r,u_i    i_r,i_i
            ↓          ↑
        +------------------+
δ(t) --→| rot(δ)   rot(-δ) |--→ |v_h|
        +------------------+
            ↓          ↑
         v_d,v_q    i_d,i_q
```
"""
function reference_frame_conversion()
    @parameters t u_r(t) u_i(t) i_d(t) i_q(t)  δ(t)
    @variables i_r(t) i_i(t) v_d(t) v_q(t) v_h(t)

    IOBlock([v_d ~ sin( δ)*u_r - cos( δ)*u_i,
             v_q ~ cos( δ)*u_r + sin( δ)*u_i,
             i_r ~ sin(-δ)*i_d - cos(-δ)*i_q,
             i_i ~ cos(-δ)*i_d + sin(-δ)*i_q,
             v_h ~ √(v_d^2 + v_q^2)],
            [δ, u_r, u_i, i_d, i_q], [i_r, i_i, v_d, v_q, v_h], name=:rfc)
end

"""
    MetaGenerator(parameters, mover, shaft, machine, AVR, PSS; constants=[], name)

Implementation of the NREL MetaGenerator model. This function takes several
`IOBlock`s as inputs which should satisfy the corresponding block specification (see below).

The `parameters` should be provided as a Dict of parameters as namespaced syms, i.e.
`mover.P_ref => 123`.
The optional `constants` parameter can be used to pass generator wide constants
as a collection of pairs, i.e. `:ω_ref => 50`. Those constants are available
as optional inputs to all of the blocks.

Block specifications: `(input) ↦ (output)`
 - mover                   `() ↦ (τ_m)`
 - shaft           `(τ_m, τ_e) ↦ (δ, ω)`
 - machine         `(v_d, v_q) ↦ (i_d, i_q, τ_e)`
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
         u_r,u_i    i_r,i_i
            ↓          ↑
     +-----------------------+ v_h
 +--→| ref. frame conversion |----→(all)   (all optional)
 |   +-----------------------+                   ↓
 |  v_d,v_q |          ↑ i_d,i_q      v_pss +---------+
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
function MetaGenerator(block_p_pairs::Vararg{BlockPara, 5}; kwargs...)
    blocks, params = mergep(block_p_pairs)
    MetaGenerator(params, blocks...; kwargs...)
end

function MetaGenerator(para::Dict, mover, shaft, machine, AVR, PSS;
                       constants=nothing, verbose=false, name=gensym(:MGenerator))
    para = copy(para) # create new reference

    verbose && println("Create MetaGenerator:")
    rfc = reference_frame_conversion()
    adder = IOComponents.Adder(name=:exitation, out=:v_f, a₁=:v_avr, a₂=:v_pss)

    # crate an constblk IOBlock which holds additional generator wide constants
    components = if constants !== nothing
        constblk = IOComponents.Constants(Dict(constants); name=:Const)
        [mover, shaft, machine, AVR, PSS, rfc, adder, constblk]
    else
        [mover, shaft, machine, AVR, PSS, rfc, adder]
    end

    @assert BlockSpec([], [:τ_m])(mover) "mover has not output τₘ!"
    @assert BlockSpec([:τ_m, :τ_e], [:δ, :ω])(shaft) "shaft not of type (τₑ, τₘ) ↦ (δ, ω)!"
    @assert BlockSpec([:v_d, :v_q], [:i_d, :i_q, :τ_e])(machine) "machine not of type (v_dq) ↦ (i_dq, τₑ)!"
    @assert BlockSpec([], [:v_avr])(AVR) "AVR has not output v_avr!"
    @assert BlockSpec([], [:v_pss])(PSS) "PSS has not output v_pss!"

    # for more flexibility everything with matching name gets connected
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
    rfc_promotions = [rfc.u_r => :u_r, rfc.u_i => :u_i]
    promotions = vcat(rfc_promotions, output_promotions)
    # we want no further autopromote to keep all of the parameters in the namespaces
    # this is necessary because the parameterdict will be namespaced as well!

    sys = IOSystem(connections, components,
                   namespace_map=promotions, outputs=[rfc.i_r, rfc.i_i],
                   autopromote=false, name=name)
    verbose && println("IOSystem before connection: ", sys)
    connected = connect_system(sys, verbose=verbose)
    verbose && println("IOBlock after connection: ", connected)

    inputs = Set(getname.(connected.inputs))
    @assert inputs == Set((:u_r, :u_i)) "System inputs [:u_r, :u_i] != $inputs. It seems like there are still open inputs in the subsystems!"

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

    para_renamespaced = typeof(para)()
    for (k, v) in para
        newk = rename(k, renamespace(connected.name, getname(k)))
        para_renamespaced[newk] = v
    end

    return BlockPara(connected, para_renamespaced)
end
