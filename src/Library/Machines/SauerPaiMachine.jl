"""
    SauerPaiMachine

Machine Model as described in the book of Sauer and Pai.

$(PowerDynamics.ref_source_file(@__FILE__, @__LINE__))
"""
@cached_component function SauerPaiMachine(; name, vf_input=true, τ_m_input=true, stator_dynamics=false)
    @named terminal = Terminal()
    # outputs
    @named δout = RealOutput() # rotor angle
    @named ωout = RealOutput() # rotor speed [pu]
    @named v_mag_out = RealOutput() # terminal voltage [pu]
    @named Pout = RealOutput() # active power [pu]
    @named Qout = RealOutput() # reactive power [pu]
    systems = [terminal, δout, ωout, v_mag_out, Pout, Qout]

    params = @parameters begin
        R_s, [description="stator resistance"]
        X_d, [description="d-axis synchronous reactance"]
        X_q, [description="q-axis synchronous reactance"]
        X′_d, [description="d-axis transient reactance"]
        X′_q, [description="q-axis transient reactance"]
        X″_d, [description="d-axis subtransient reactance"]
        X″_q, [description="q-axis subtransient reactance"]
        X_ls, [description="stator leakage reactance"]
        T′_d0, [description="d-axis transient time constant"]
        T″_d0, [description="d-axis subtransient time constant"]
        T′_q0, [description="q-axis transient time constant"]
        T″_q0, [description="q-axis subtransient time constant"]
        H, [description="inertia constant"]
        D=0, [description="direct shaft damping"]
        Sbase, [bound_to = :systembase₊Sbase, description="System power base [MVA]"]
        Vbase, [bound_to = :busbar₊Vbase, description="Bus voltage base [kV]"]
        ωbase, [bound_to = :systembase₊ωbase, description="System frequency base [rad/s]"]
        ωframe, [bound_to = :systembase₊ωframe, description="Global dq frame speed [pu]"]
        Sn, [initf_weak = Sbase, description="Machine power rating [MVA]"]
        Vn, [initf_weak = Vbase, description="Machine voltage rating [kV]"]
    end
    vars = @variables begin
        ψ_d(t), [guess=1, description="d-axis flux linkage"]
        ψ_q(t), [guess=0, description="q-axis flux linkage"]
        ψ″_d(t), [guess=1, description="flux linkage assosciated with X″_d"]
        ψ″_q(t), [guess=0, description="flux linkage assosciated with X″_q"]
        I_d(t), [guess=0, description="d-axis current"]
        I_q(t), [guess=0, description="q-axis current"]
        V_d(t), [guess=0, description="d-axis voltage"]
        V_q(t), [guess=1, description="q-axis voltage"]
        E′_d(t), [guess=0, description="transient voltage behind transient reactance in d-axis"]
        E′_q(t), [guess=1, description="transient voltage behind transient reactance in q-axis"]
        δ(t), [guess=0, description="rotor angle"]
        ω(t), [guess=1, description="rotor speed"]
        τ_e(t), [bounds=(0,Inf) ,description="electrical torque"]
        # observables
        v_mag(t), [description="terminal voltage [machine pu]"]
        v_arg(t), [description="Generator terminal angle"]
        P(t), [description="active power [machine pu]"]
        Q(t), [description="reactive power [machine pu]"]
        # inputs/parameters
        vf(t), [bounds=(0,Inf), description="field voltage"]
        τ_m(t), [bounds=(0,Inf), description="mechanical torque"]
    end

    # input/parameter switches
    if vf_input
        @named vf_in = RealInput(guess=1) # field voltage input [pu]
        push!(systems, vf_in)
        vf_source = vf_in.u
    else
        @parameters vf_set [guess=1, bounds=(0,Inf), description="field voltage"]
        push!(params, vf_set)
        vf_source = vf_set
    end
    if τ_m_input
        @named τ_m_in = RealInput() # mechanical torque [pu]
        push!(systems, τ_m_in)
        τ_m_source = τ_m_in.u
    else
        @parameters τ_m_set [guess=1, bounds=(0,Inf), description="mechanical torque"]
        push!(params, τ_m_set)
        τ_m_source = τ_m_set
    end

    γ_d1 = (X″_d - X_ls)/(X′_d - X_ls)
    γ_q1 = (X″_q - X_ls)/(X′_q - X_ls)
    γ_d2 = (X′_d-X″_d)/(X′_d-X_ls)^2 # ~ (1 - γ_d1)/(X′_d - X_ls)
    γ_q2 = (X′_q-X″_q)/(X′_q-X_ls)^2 # ~ (1 - γ_q1)/(X′_q - X_ls)
    # This transformation seems identical to Milano and PSD models
    T_to_loc(α)  = [ sin(α) -cos(α);
                     cos(α)  sin(α)]
    T_to_glob(α) = [ sin(α)  cos(α);
                    -cos(α)  sin(α)]

    eqs = [
        # Park's transformations
        [terminal.u_r, terminal.u_i] .~ T_to_glob(δ)*[V_d, V_q] * (Vn/Vbase)
        [I_d, I_q] .~ T_to_loc(δ)*[terminal.i_r, terminal.i_i] * (Ibase(Sbase, Vbase)/Ibase(Sn, Vn))

        τ_e ~ ψ_d*I_q - ψ_q*I_d
        # for static ψ, this becomes which makes sense!
        # τ_e ~  (P + R_s*(I_d^2 + I_q^2))/ω

        # δ is measured against the global dq frame → ωframe; the `1` in the damping
        # term is nominal speed (a setpoint), which is frame-invariant.
        Dt(δ) ~ ωbase*(ω - ωframe)
        2*H * Dt(ω) ~ τ_m  - τ_e - D*(ω - 1)

        T′_d0 * Dt(E′_q) ~ -E′_q - (X_d - X′_d)*(I_d - γ_d2*ψ″_d - (1-γ_d1)*I_d + γ_d2*E′_q) + vf
        T′_q0 * Dt(E′_d) ~ -E′_d + (X_q - X′_q)*(I_q - γ_q2*ψ″_q - (1-γ_q1)*I_q - γ_q2*E′_d)
        T″_d0 * Dt(ψ″_d) ~ -ψ″_d + E′_q - (X′_d - X_ls)*I_d
        T″_q0 * Dt(ψ″_q) ~ -ψ″_q - E′_d - (X′_q - X_ls)*I_q

        # this constraint essentialy forces ψ_d and ψ_q
        ψ_d ~ -X″_d*I_d + γ_d1*E′_q + (1-γ_d1)*ψ″_d
        ψ_q ~ -X″_q*I_q - γ_q1*E′_d + (1-γ_q1)*ψ″_q
        # I_d ~ (-ψ_d + γ_d1*E′_q + (1-γ_d1)*ψ″_d)/X″_d
        # I_q ~ (-ψ_q - γ_q1*E′_d + (1-γ_q1)*ψ″_q)/X″_q
        # 0 ~ -X″_d*I_d + γ_d1*E′_q + (1-γ_d1)*ψ″_d - ψ_d
        # 0 ~ -X″_q*I_q - γ_q1*E′_d + (1-γ_q1)*ψ″_q - ψ_q

        # inputs
        vf ~ vf_source
        τ_m ~ τ_m_source

        # observables
        v_mag ~ sqrt(V_d^2 + V_q^2)
        v_arg ~ atan(V_q, V_d)
        P ~ V_d*I_d + V_q*I_q
        Q ~ V_q*I_d - V_d*I_q

        #outputs
        Pout.u ~ P
        Qout.u ~ Q
        v_mag_out.u ~ v_mag
        δout.u ~ δ
        ωout.u ~ ω
    ]
    if stator_dynamics
        # stator equations, in the *rotor* frame: the cross-coupling speed is the
        # rotor speed ω (a state), not the global frame. The bare 1/ωbase is the
        # flux base Vbase/ωbase, i.e. the unit of the differentiated quantity.
        append!(eqs, [1/ωbase * Dt(ψ_d) ~ R_s*I_d + ω * ψ_q + V_d,
                      1/ωbase * Dt(ψ_q) ~ R_s*I_q - ω * ψ_d + V_q])
    else
        # static fomulation
        append!(eqs, [0 ~ R_s*I_d + ω * ψ_q + V_d,
                      0 ~ R_s*I_q - ω * ψ_d + V_q])
    end
    System(eqs, t, vars, params; name, systems)
end

#=
MWE why whe need to have ψ force by constraint rather than free
in a nutshell: because the input i will force psi anyway, and if
and then the voltage depends on the derivative of I

@mtkmodel SauerPaiMachine2 begin
    @parameters begin
        I_d, [guess=0, description="d-axis current"]
        I_q, [guess=0, description="q-axis current"]
    end
    @variables begin
        ψ_d(t), [description="d-axis flux linkage"]
        ψ_q(t), [description="q-axis flux linkage"]
        V_d(t), [guess=0, description="d-axis voltage"]
        V_q(t), [guess=1, description="q-axis voltage"]
    end
    @equations begin
        Dt(ψ_d) ~ I_d + ψ_q + V_d
        Dt(ψ_q) ~ I_q - ψ_d + V_q
        0 ~ -I_d - ψ_d
        0 ~ -I_q - ψ_q
    end
end
VertexModel(Library.SauerPaiMachine2(name=:machine), [:I_d, :I_q], [:V_d, :V_q])
=#
