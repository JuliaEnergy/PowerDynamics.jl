@mtkmodel SauerPaiMachine begin
    @structural_parameters begin
        vf_input = true
        τ_m_input = true
    end
    @components begin
        terminal=Terminal()
        # inputs
        if vf_input
            vf_in = RealInput(guess=1) # field voltage input [pu]
        end
        if τ_m_input
            τ_m_in = RealInput() # mechanical torque [pu]
        end
        # outputs
        δout = RealOutput() # rotor angle
        ωout = RealOutput() # rotor speed [pu]
        v_mag_out = RealOutput() # terminal voltage [pu]
        Pout = RealOutput() # active power [pu]
        Qout = RealOutput() # reactive power [pu]
    end
    @parameters begin
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
        # System and machine base
        S_b, [description="System power basis in MVA"]
        V_b, [description="System voltage basis in kV"]
        ω_b, [description="System base frequency in rad/s"]
        Sn=S_b, [description="Machine power rating in MVA"]
        Vn=V_b, [description="Machine voltage rating in kV"]
        # input/parameter switches
        if !vf_input
            vf_set, [guess=1, description="field voltage"]
        end
        if !τ_m_input
            τ_m_set, [guess=1, description="mechanical torque"]
        end
    end
    @variables begin
        ψ_d(t), [description="d-axis flux linkage"]
        ψ_q(t), [description="q-axis flux linkage"]
        ψ″_d(t), [guess=1, description="flux linkage assosciated with X″_d"]
        ψ″_q(t), [guess=0, description="flux linkage assosciated with X″_q"]
        I_d(t), [guess=0, description="d-axis current"]
        I_q(t), [guess=0, description="q-axis current"]
        V_d(t), [guess=0, description="d-axis voltage"]
        V_q(t), [guess=1, description="q-axis voltage"]
        E′_d(t), [guess=1, description="transient voltage behind transient reactance in d-axis"]
        E′_q(t), [guess=0, description="transient voltage behind transient reactance in q-axis"]
        δ(t), [guess=0, description="rotor angle"]
        ω(t), [guess=1, description="rotor speed"]
        τ_e(t), [description="electrical torque"]
        # observables
        v_mag(t), [description="terminal voltage [machine pu]"]
        v_arg(t), [description="Generator terminal angle"]
        P(t), [description="active power [machine pu]"]
        Q(t), [description="reactive power [machine pu]"]
        # inputs/parameters
        vf(t), [description="field voltage"]
        τ_m(t), [description="mechanical torque"]
    end
    begin
        γ_d1 = (X″_d - X_ls)/(X′_d - X_ls)
        γ_q1 = (X″_q - X_ls)/(X′_q - X_ls)
        γ_d2 = (X′_d-X″_d)/(X′_d-X_ls)^2 # ~ (1 - γ_d1)/(X′_d - X_ls)
        γ_q2 = (X′_q-X″_q)/(X′_q-X_ls)^2 # ~ (1 - γ_q1)/(X′_q - X_ls)
        # This transformation seems identical to Milano and PSD models
        T_to_loc(α)  = [sin(α) -cos(α);  cos(α) sin(α)]
        T_to_glob(α) = [sin(α)  cos(α); -cos(α) sin(α)]
        # does not work, looks even strange. Steady satet is no steady state
        # T_to_loc(α)  = [ cos(α) sin(α); sin(α) -cos(α)];
        # T_to_glob(α) = [ cos(α) sin(α); sin(α) -cos(α)];
    end
    @equations begin
        # Park's transformations
        [terminal.u_r, terminal.u_i] .~ T_to_glob(δ)*[V_d, V_q] * Vn/V_b
        [I_d, I_q] .~ T_to_loc(δ)*[terminal.i_r, terminal.i_i] * Ibase(S_b, V_b)/Ibase(Sn, Vn)

        τ_e ~ ψ_d*I_q - ψ_q*I_d
        # for static ψ, this becomes which makes sense!
        # τ_e ~  (P + R_s*(I_d^2 + I_q^2))/ω

        Dt(δ) ~ ω_b*(ω - 1)
        2*H * Dt(ω) ~ τ_m  - τ_e - D*(ω - 1)

        # stator equations
        # 1/ω_b * Dt(ψ_d) ~ R_s*I_d + ω * ψ_q + V_d
        # 1/ω_b * Dt(ψ_q) ~ R_s*I_q - ω * ψ_d + V_q
        # static fomulation
        # 0 ~ R_s*I_d + ω * ψ_q + V_d
        # 0 ~ R_s*I_q - ω * ψ_d + V_q
        # static formualion in V_d, V_q which is the only free stuff
        V_d ~ -R_s*I_d - ω * ψ_q
        V_q ~ -R_s*I_q + ω * ψ_d

        T′_d0 * Dt(E′_q) ~ -E′_q - (X_d - X′_d)*(I_d - γ_d2*ψ″_d - (1-γ_d1)*I_d + γ_d2*E′_q) + vf
        T′_q0 * Dt(E′_d) ~ -E′_d - (X_q - X′_q)*(I_q - γ_q2*ψ″_q - (1-γ_q1)*I_q - γ_q2*E′_d)
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
        vf ~ vf_input ? vf_in.u : vf_set
        τ_m ~ τ_m_input ? τ_m_in.u : τ_m_set

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
    end
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
