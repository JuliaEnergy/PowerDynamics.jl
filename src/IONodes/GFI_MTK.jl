using BlockSystems
using MarinePowerDynamics.IOComponents

export GFI
"""
    GFI(;τ_v,τ_P,τ_Q,K_P,K_Q,V_r,P,Q,ω_r)

Implementation of `VSIVoltagePT1` as an IONode.

                +----------+ +---------+
        +-----+ | p_filter |-| p_droop |  +----------+
    i_r-|    Q|-|   τ_P    | | P, ω_r  |--|          |
    i_i-|     | +----------+ +---------+ ω| v_source |--+-u_r
        | pow |                           |   τ_v    |  |
     +--|     | +----------+ +---------+ V|          |-+|-u_i
     |+-|    P|-| q_filter |-| q_droop |--|          | ||
     || +-----+ |   τ_Q    | | Q, V_r  |  +----------+ ||
     ||         +----------+ +---------+               ||
     |+------------------------------------------------+|
     +--------------------------------------------------+
"""
function GFI(;τ_v,τ_P,τ_Q,K_P,K_Q,V_r,P,Q,ω_r, name=gensym(:GridForming))

    para = Dict(
        :τ_v => τ_v,    # time constant voltage control delay
        :τ_P => τ_P,    # time constant active power measurement
        :τ_Q => τ_Q,    # time constant reactive power measurement
        :K_P => K_P,    # droop constant frequency droop
        :K_Q => K_Q,    # droop constant voltage droop
        :V_r => V_r,    # reference/ desired voltage
        :P   => P,      # active (real) power infeed
        :Q   => Q,      # reactive (imag) power infeed                .
        :ω_r => ω_r)    # refrence/ desired frequency

    p_filter = LowPassFilter(name=:p_filter, τ=:τ_P)
    q_filter = LowPassFilter(name=:q_filter, τ=:τ_Q)
    p_droop = DroopControl(name=:p_droop, x_ref=:P, u_ref=:ω_r, K=:K_P)
    q_droop = DroopControl(name=:q_droop, x_ref=:Q, u_ref=:V_r, K=:K_Q)
    v_source = VoltageSource(name=:v_source, τ=:τ_v)
    pow = Power(name=:pow)

    gfi = IOSystem([pow.P => p_filter.input,
                    pow.Q => q_filter.input,
                    p_filter.output => p_droop.x,
                    q_filter.output => q_droop.x,
                    p_droop.u => v_source.ω,
                    q_droop.u => v_source.V,
                    v_source.u_r => pow.u_r,
                    v_source.u_i => pow.u_i],
                   [pow, p_filter, p_droop, q_filter, q_droop, v_source],
                   name=name,
                   outputs=[v_source.u_i, v_source.u_r],
                   namespace_map=[p_droop.P => :P,
                                  q_droop.Q => :Q])
    # The symbosl P and Q are defined booth in the droops and in the pow block
    # therefore we had to provide this part of the namespace_map manually

    connected = connect_system(gfi)
    IONode(connected, para)
end
