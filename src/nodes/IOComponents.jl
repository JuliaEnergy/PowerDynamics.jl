module IOComponents

using BlockSystems
using ModelingToolkit

export LowPassFilter, DroopControl, VoltageSource, Power

"""
    LowPassFilter(;name, renamings...)

Returns a low pass filter. The name of the system and the names of the vars
can be changed with keyword arguments `name=:myname, τ=:mytau, …`.

    out'(t) = 1/τ (in(t) - out(t))

               +-----+
    input(t) --|  τ  |-- output(t)
               +-----+

    IOBlock :##LPF# with 1 eqs
    ├ inputs:  input(t)
    ├ outputs: output(t)
    ├ istates: (empty)
    └ iparams: τ
"""
function LowPassFilter(;name=gensym(:lpf), renamings...)
    @parameters t τ
    @parameters input(t)
    @variables output(t)
    D = Differential(t)

    block = IOBlock([D(output) ~ 1/τ * (- output + input)],
                [input], [output]; name)
    return isempty(renamings) ? block : rename_vars(block; renamings...)
end

"""
    DroopControl(;name, renamings...)

Returns a DroopControl. The name of the system and the names of the vars
can be changed with keyword arguments `name=:myname, K=:myK, …`.

    u = - K*(x - x_ref) + u_ref

           +-----------------+
    x(t) --| K, x_ref, u_ref |-- u(t)
           +-----------------+

    IOBlock :##droop# with 1 eqs
    ├ inputs:  x(t)
    ├ outputs: u(t)
    ├ istates: (empty)
    └ iparams: K, x_ref, u_ref
"""
function DroopControl(;name=gensym(:droop), renamings...)
    @parameters t K x_ref u_ref
    @parameters x(t)
    @variables u(t)
    D = Differential(t)

    block = IOBlock([u ~ - K * (x - x_ref) + u_ref], # output is the droop voltage v
                    [x], [u]; name)

    return isempty(renamings) ? block : rename_vars(block; renamings...)
end

"""
    VoltageSource(;name, renamings...)

Returns a VoltageSource Block. Models the complex voltage dynamic as a low pass
inspired by [Schiffer et. al.](http://eprints.whiterose.ac.uk/92371/1/schiffer_etal_automatica_2014.pdf)
for a reference frequency ω and voltage magnitude V.

       A = 1/τ ⋅ (V/√(uᵢ² + uᵣ²) - 1)
    u_r' = -ω uᵢ + A uᵣ
    u_i' =  ω uᵣ + A uᵢ

           +-----+
    ω(t) --|  τ  |-- u_r(t)
    V(t) --|     |-- u_i(t)
           +-----+

    IOBlock :##vsource# with 3 eqs
    ├ inputs:  ω(t), V(t)
    ├ outputs: u_i(t), u_r(t)
    ├ istates: A(t)
    └ iparams: τ
"""
function VoltageSource(;name=gensym(:vsource), renamings...)
    @parameters t τ
    @parameters ω(t) V(t)
    @variables u_i(t) u_r(t) A(t)
    D = Differential(t)

    block = IOBlock([A ~ 1/τ * (V/√(u_i^2 + u_r^2) - 1),
                     D(u_r) ~ -ω * u_i + A*u_r,
                     D(u_i) ~  ω * u_r + A*u_i],
                    [ω, V], [u_i, u_r]; name)
    # TODO: reduce_eqs(voltage_source) once implemented in blocksystems to reduce the
    # internal algebraic state A
    return isempty(renamings) ? block : rename_vars(block; renamings...)
end

"""
    Power(;name, renamings...)

Returns a Block which calculates the active and reactive power for a given complex input.

    P = uᵣ iᵣ + uᵢ iᵢ
    Q = uᵢ iᵣ - uᵣ iᵢ

             +-----+
    u_r(t) --|     |-- P(t)
    u_i(t) --|     |
    i_r(t) --|     |
    i_i(t) --|     |-- Q(t)
             +-----+

    IOBlock :##power# with 2 eqs
    ├ inputs:  u_i(t), u_r(t), i_i(t), i_r(t)
    ├ outputs: P(t), Q(t)
    ├ istates: (empty)
    └ iparams: (empty)
"""
function Power(;name=gensym(:power), renamings...)
    @parameters t
    @parameters u_i(t) u_r(t) i_i(t) i_r(t)
    @variables P(t) Q(t)

    block = IOBlock([P ~ u_r*i_r + u_i*i_i,
                     Q ~ u_i*i_r - u_r*i_i],
                    [u_i, u_r, i_i, i_r], [P, Q]; name)

    return isempty(renamings) ? block : rename_vars(block; renamings...)
end

end # module
