#=
## Generate `PowerDynamics.jl`-Node with `BlockSystems`
We want to model the
[`VSIVoltagePT1`](https://juliaenergy.github.io/PowerDynamics.jl/dev/node_types/#PowerDynamics.VSIVoltagePT1)
with the help of `BlockSystems`.
We start by defining some general stuff...
=#
using BlockSystems
using ModelingToolkit

function GFI(;τ_v,τ_P,τ_Q,K_P,K_Q,V_r,P,Q,ω_r)

    para = Dict(
    :τ_v => τ_v,    # time constant voltage control delay
    :τ_P => τ_P,      # time constant active power measurement
    :τ_Q => τ_Q,      # time constant reactive power measurement
    :K_P => K_P,     # droop constant frequency droop
    :K_Q => K_Q,     # droop constant voltage droop
    :V_r => V_r,   # reference/ desired voltage
    :P   => P, # active (real) power infeed
    :Q   => Q, # reactive (imag) power infeed                .
    :ω_r => ω_r)   # refrence/ desired frequency

    
    @parameters t
    D = Differential(t)


    # #### low pass filter
    @parameters τ input(t)
    @variables filtered(t)

    lpf = IOBlock([D(filtered) ~ 1/τ * (- filtered + input)],
                [input], [filtered])

    # #### voltage source

    @parameters ω(t) v(t) τ
    @variables u_i(t) u_r(t) A(t)

    ## explicit algebraic equation for A will be reduced at connect
    voltage_source = IOBlock([A ~ 1/τ * (v/√(u_i^2 + u_r^2) - 1),
                            D(u_r) ~ -ω * u_i + A*u_r,
                            D(u_i) ~  ω * u_r + A*u_i],
                            [ω, v], [u_i, u_r])

    # #### Droop control

    @parameters K u_ref x_ref x(t)
    @variables u(t)

    droop_control = IOBlock([
        u ~ - K * (x - x_ref) + u_ref # output is the droop voltage v
        ], [x], [u])



    p_filter = IOBlock(lpf, name = :p_filter)
    q_filter = IOBlock(lpf, name = :q_filter)
    p_droop = IOBlock(droop_control, name = :p_droop)
    q_droop = IOBlock(droop_control, name = :q_droop)
    v_source = IOBlock(voltage_source, name = :v_source)


    gfi = IOSystem([p_filter.filtered => p_droop.x,
                    q_filter.filtered => q_droop.x,
                    p_droop.u => v_source.ω,
                    q_droop.u => v_source.v],
                [p_filter, q_filter, p_droop, q_droop, v_source],
                name = :GridForming,
                namespace_map = [p_filter.input => :P_in,
                                    q_filter.input => :Q_in,
                                    p_filter.filtered => :p_filtered,
                                    q_filter.filtered => :q_filtered,
                                    ## parameter names which match VSIVoltagePT1
                                    v_source.τ => :τ_v, # time constant voltage control delay
                                    p_filter.τ => :τ_P, # time constant active power measurement
                                    q_filter.τ => :τ_Q, # time constant reactive power measurement
                                    p_droop.K  => :K_P, # droop constant frequency droop
                                    q_droop.K  => :K_Q, # droop constant voltage droop
                                    q_droop.u_ref => :V_r, # reference/ desired voltage
                                    p_droop.u_ref => :ω_r, # reference/ desired frequency
                                    p_droop.x_ref => :P, # active (real) power infeed
                                    q_droop.x_ref => :Q], # reactive (imag) power infeed                .
                outputs = [v_source.u_i, v_source.u_r])


    @parameters u_i(t) u_r(t) i_i(t) i_r(t)
    @variables P_in(t) Q_in(t)
    pow = IOBlock([P_in ~ u_r*i_r + u_i*i_i,
                Q_in ~ u_i*i_r - u_r*i_i],
                [u_i, u_r, i_i, i_r], [P_in, Q_in], name=:pow)

    gfi2 = IOSystem([gfi.u_i => pow.u_i,
                    gfi.u_r => pow.u_r,
                    pow.P_in => gfi.P_in,
                    pow.Q_in => gfi.Q_in],
                    [gfi, pow], outputs=[gfi.u_i, gfi.u_r], name=:gfi_i_to_u)


    # and the system can be reduced to a single IOBlock
    connected = connect_system(gfi2)

    IONode(connected, para)
end


export GFI
