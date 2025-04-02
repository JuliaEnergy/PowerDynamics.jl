#=
* This file is based on the PowerFactory Standard Model for synchronous machines (PowerFactory 2023 Technical Reference - Synchronous Machine by DIgSILENT GmbH)
=#

@mtkmodel StandardModel_pf begin
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
        #primary parameters
        R_s, [description="stator resistance"]
        X_rld, [description="coupling reactance between field and damper winding"]
        X_rlq, [description="coupling reactance between q-axis damper windings"]
        X_d, [description="d-axis synchronous reactance"]
        X_q, [description="q-axis synchronous reactance"]
        X′_d, [description="d-axis transient reactance"]
        X′_q=0.0001, [description="q-axis transient reactance, not needed for salient pole machine"]
        X″_d, [description="d-axis subtransient reactance"]
        X″_q, [description="q-axis subtransient reactance"]
        X_ls, [description="stator leakage reactance"]
        T′_d0, [description="d-axis transient time constant"]
        T″_d0, [description="d-axis subtransient time constant"]
        T′_q0=0, [description="q-axis transient time constant, not needed for salient pole machine"]
        T″_q0, [description="q-axis subtransient time constant"]
        H, [description="inertia constant"]
        S_b, [description="System power basis in MVA"]
        V_b, [description="System voltage basis in kV"]
        ω_b, [description="System base frequency in rad/s"]
        Sn, [description="Machine power rating in MVA"]
        Vn, [description="Machine voltage rating in kV"]
        cosn, [description="rated power factor - ??"]
        n_ref=1, [description="nominal speed (1 p.u.) of the machine or the speed of the local reference machine"]
        dkd, [description="Damping torque coefficient"]
        dpe, [description="Damping torque coefficient based on power"]
        D, [description="Damping constant"]
        salientpole=1, [description="salient pole or round-rotor machine"]
        xmdm=0, [description="Torque Input input signal in p.u."]
        addmt=0, [description="Additional Torque parameter in p.u."]
        pt, [guess=1, bounds=(0,Inf), description="Turbine Power input signal in p.u."]
        dpu=0, [description="dpu * n is turbine shaft friction torque in p.u.;"]
        speedvar=0, [description="speed variations considered in stator voltage equations"]
        speedvar_part=1, [description="speed variations partly neglected in stator voltage equations"]
        # input/parameter switches
        if !vf_input
            vf_set, [guess=1, bounds=(0,Inf), description="field voltage"]
        end
        if !τ_m_input
            τ_m_set, [guess=1, bounds=(0,Inf), description="mechanical torque"]
        end
        if speedvar==1
            speedvar_part=0
        end
    end
    @variables begin
        ψ_d(t), [description="d-axis flux linkage"]
        ψ_q(t), [description="q-axis flux linkage"]
        ψ″_d(t), [guess=1, description="flux linkage assosciated with X″_d"]
        ψ″_q(t), [guess=0, description="flux linkage assosciated with X″_q"]
        I_d(t), [description="d-axis current"]
        I_q(t), [description="q-axis current"]
        V_d(t), [description="d-axis voltage"]
        V″_d(t), [description="subtransient d-axis voltage"]
        V″_q(t), [description="subtransient q-axis voltage"]
        V_q(t), [description="q-axis voltage"]
        E′_d(t), [guess=0, description="transient voltage behind transient reactance in d-axis"]
        E′_q(t), [guess=1, description="transient voltage behind transient reactance in q-axis"]
        δ(t), [guess=0, description="rotor angle"]
        ϕ(t), [guess=0, description="angle between d-axis and reference voltage of network"]
        ω(t), [guess=1, description="rotor speed"]
        ψ_fd(t), [guess=1, description=" flux linkage"]
        ψ_1d(t), [guess=1, description=" flux linkage"]
        ψ_1q(t), [guess=0, description=" flux linkage"]
        ψ_2q(t), [guess=0,description=" flux linkage"]
        I_fd(t), [description="d-axis current"]
        I_1d(t), [description="d-axis current"]
        I_1q(t), [description="d-axis current"]
        I_2q(t), [description="d-axis current"]
        # observables
        v_mag(t), [description="terminal voltage [machine pu]"]
        v_arg(t), [description="Generator terminal angle"]
        P(t), [description="active power [machine pu]"]
        Q(t), [description="reactive power [machine pu]"]
        # inputs/parameters
        vf(t), [description="field voltage"]
        τ_m(t), [description="mechanical torque"]
        τ_e(t), [description="electrical torque"]
        τ_dkd(t), [description="electrical torque"]
        τ_dpe(t), [description="damping torque based on power"]
        τ_ag(t), [description="acceleration time constant in s"]
        n(t), [guess=1, description="rotor speed"]
        #secondary para,ezers
        T″_d(t)
        T″_q(t)
        T′_d(t)
        T′_q(t)
        X_ad(t), [description="Mutual (magnetising) reactance, d-axis"]
        X_aq(t), [description="Mutual (magnetising) reactance, q-axis"]
        X_1(t)
        X_2(t)
        X_3(t)
        T_1(t)
        T_2(t)
        a(t)
        b(t)
        T_σfd(t)
        T_σ1d(t)
        X_fd(t), [description="Reactance of excitation (field) winding (d-axis)"]
        X_1d(t), [description="Reactance of 1d-damper winding (d-axis)"]
        R_fd(t), [description="Resistance of excitation winding (d-axis)"]
        R_1d(t), [description="Resistance of 1d-damper winding (d-axis)"]
        X_4(t)
        X_5(t)
        X_6(t)
        T_3(t)
        T_4(t)
        c(t)
        d(t)
        T_σ2q(t)
        T_σ1q(t)
        X_2qr(t)
        X_1qr(t)
        R_2qr(t)
        R_1qr(t)
        X_1qs(t)
        R_1qs(t)
        X_1q(t), [description="Reactance of 1q-damper winding (q-axis)"]
        R_1q(t), [description="Resistance of 1q-damper winding (q-axis)"]
        X_2q(t), [description="Reactance of 2q-damper winding (q-axis)"]
        R_2q(t), [description="Resistance of 2q-damper winding (q-axis)"]
        k_fd(t)
        k_1d(t)
        X″_d_new(t)
        k_1qs(t)
        X″_qs(t)
        k_1qr(t)
        k_2qr(t)
        X″_qr(t)
        k_1q(t)
        k_2q(t)
        X″_q_new(t)
        X_det_d(t)
        X_det_q(t)
        X_fd_loop(t)
        X_1d_loop(t)
        X_1q_loop(t)
        X_2q_loop(t)

    end
    begin
        T_park(α) = [sin(α) cos(α); -cos(α) sin(α)]
        T_to_loc(α)  = [ sin(α) -cos(α);
                         cos(α)  sin(α)]
        T_to_glob(α) = [ sin(α)  cos(α);
                        -cos(α)  sin(α)]
    end
    @equations begin
        #secondary parameters calculation
        T″_d ~ T″_d0 * X″_d/X′_d
        T″_q ~ T″_q0 * X″_q/(X′_q * (1-salientpole) + salientpole * X_q)
        T′_d ~ T′_d0 * X′_d/X_d
        T′_q ~ T′_q0 * X′_q/X_q
        X_ad ~ X_d - X_ls
        X_aq ~ X_q - X_ls
        X_1 ~ X_d - X_ls + X_rld
        X_2 ~ X_1 - (X_d - X_ls)^2 / X_d
        X_3 ~ (X_2 - X_1 * X″_d/X_d) / (1 - X″_d/X_d)
        T_1 ~ X_d / X′_d * T′_d + (1 - X_d/X′_d + X_d/X″_d) * T″_d
        T_2 ~ T″_d + T′_d
        a ~ (X_2 * T_1 - X_1 * T_2) / (X_1 - X_2)
        b ~ X_3 * T″_d * T′_d / (X_3 - X_2)
        T_σfd ~ -a/2 + sqrt(a^2/4 - b)
        T_σ1d ~ -a/2 - sqrt(a^2/4 - b)
        X_fd ~ (T_σfd - T_σ1d) / ((T_1 - T_2)/(X_1 - X_2) + T_σ1d / X_3)
        X_1d ~ (T_σ1d - T_σfd) / ((T_1 - T_2)/(X_1 - X_2) + T_σfd / X_3)
        R_fd ~ X_fd / (ω_b * T_σfd)
        R_1d ~ X_1d / (ω_b * T_σ1d)
        X_4 ~ X_q - X_ls + X_rlq
        X_5 ~ X_4 - (X_q - X_ls)^2 / X_q
        X_6 ~ (X_5 - X_4 * X″_q/X_q) / (1 - X″_q/X_q)
        T_3 ~ X_q / X′_q * T′_q + (1 - X_q/X′_q + X_q/X″_q) * T″_q
        T_4 ~ T″_q + T′_q
        c ~ (X_5 * T_3 - X_4 * T_4) / (X_4 - X_5)
        d ~ X_6 * T″_q * T′_q / (X_6 - X_5)
        T_σ2q ~ -c/2 + sqrt(c^2/4 - d)
        T_σ1q ~ -c/2 - sqrt(c^2/4 - d)
        X_2qr ~ (T_σ2q - T_σ1q) / ((T_3 - T_4)/(X_4 - X_5) + T_σ1q / X_6) #round-rotor
        X_1qr ~ (T_σ1q - T_σ2q) / ((T_3 - T_4)/(X_4 - X_5) + T_σ2q / X_6) #round-rotor
        R_2qr ~ X_2qr / (ω_b * T_σ2q) #round-rotor
        R_1qr ~ X_1qr / (ω_b * T_σ1q) #round-rotor
        X_1qs ~ (X_q - X_ls) * (X″_q - X_ls) / (X_q - X″_q) #salient pole
        R_1qs ~ X″_q / X_q * (X_q - X_ls + X_1qs) / (ω_b * T″_q) #salient pole
        X_1q ~ salientpole * X_1qs + (1-salientpole) * X_1qr
        R_1q ~ salientpole * R_1qs + (1-salientpole) * R_1qr
        X_2q ~ (1-salientpole) * X_2qr
        R_2q ~ (1-salientpole) * R_2qr
        k_fd ~ (X_ad * X_1d) / ((X_ad + X_rld) * (X_1d + X_fd) + X_fd * X_1d)
        k_1d ~ (X_ad * X_fd) / ((X_ad + X_rld) * (X_1d + X_fd) + X_fd * X_1d)
        X″_d_new ~ X_ad + X_ls - (k_1d + k_fd) * X_ad
        k_1qs ~ X_aq / (X_aq + X_rlq + X_1q) #salient pole
        X″_qs ~ X_aq + X_ls - k_1qs * X_aq #salient pole
        k_1qr ~ (X_aq * X_2q)/((X_aq + X_rlq) * (X_2q + X_1q) + X_2q * X_1q) #round rotor
        k_2qr ~ (X_aq * X_1q)/((X_aq + X_rlq) * (X_2q + X_1q) + X_2q * X_1q) #round rotor
        X″_qr ~ X_aq + X_ls - (k_2qr + k_1qr) * X_aq #round rotor
        k_1q ~ salientpole * k_1qs + (1-salientpole) * k_1qr
        k_2q ~ (1-salientpole) * k_2qr
        X″_q_new ~ salientpole * X″_qs + (1-salientpole)* X″_qr
        X_det_d ~ (X_ad + X_rld) * (X_1d + X_fd) + X_fd * X_1d
        X_det_q ~ (X_aq + X_rlq) * (X_2q + X_1q) + X_2q * X_1q
        X_fd_loop ~ X_ad + X_rld + X_fd
        X_1d_loop ~ X_ad + X_rld + X_1d
        X_1q_loop ~ X_aq + X_rlq + X_1q
        X_2q_loop ~ X_aq + X_rlq + X_2q
        # Park's transformations
        [terminal.u_r, terminal.u_i] .~ T_to_glob(δ)*[V_d, V_q] * Vn/V_b
        [I_d, I_q] .~ T_to_loc(δ)*[terminal.i_r, terminal.i_i] * Ibase(S_b, V_b)/Ibase(Sn, Vn)

        #stator flux equations (55), (60), (59)
        ψ_d ~ -(X_ls + X_ad) * I_d + X_ad * I_fd + X_ad * I_1d #(55)
        ψ_q ~ -(X_ls + X_aq) * I_q + X_aq * I_2q + X_aq * I_1q
        ψ_d ~ -X″_d * I_d + ψ″_d #(60)
        ψ_q ~ -X″_q * I_q + ψ″_q
        #ψ″_d ~ k_fd * ψ_fd + k_1d * ψ_1d #(59)
        #ψ″_q ~ k_1q * ψ_1q + k_2q * ψ_2q

        #stator voltage equations (72) (73) (RMS Modell, instead of (61), (62) resp. (54))
        V_d ~ V″_d - R_s * I_d + X″_q * I_q * ((1-speedvar) + n * speedvar) #if neglected: n set to initial speed in this equation, here 1
        V_q ~ V″_q - R_s * I_q - X″_d * I_d * ((1-speedvar) + n * speedvar)
        V″_d ~ - ψ″_q * (n* (speedvar_part + speedvar) + (1-(speedvar + speedvar_part)))
        V″_q ~  ψ″_d * (n*(speedvar_part + speedvar) + (1-(speedvar_part + speedvar)))

        #electrical torque (66)
        τ_e * cosn ~ I_q * ψ_d - I_d * ψ_q

        #mechanical equation motor (103), (104), (105)
        τ_dkd ~ dkd * (n - n_ref)
        τ_dpe ~ dpe/n * (n - n_ref)
        #τ_ag ~ 2 * H * 100 / (S_b * cosn)
        τ_ag ~ 2 * H
        Dt(n) ~ (τ_m - τ_e - τ_dkd - τ_dpe) / τ_ag

        #Dt(ϕ) ~ ω_b * (n - ω_b/(2*π)) #(115)
        #δ ~ ϕ + π/2 #- phiu
        Dt(δ) ~ ω_b * (n - 1)

        #rotor flux linkage (58)
        #ψ_fd ~ -X_ad * I_d + (X_ad + X_rld + X_fd) * I_fd + (X_ad + X_rld) * I_1d
        #ψ_1d ~ -X_ad * I_d + (X_ad + X_rld) * I_fd + (X_ad + X_rld + X_1d) * I_1d
        #ψ_1q ~ -X_aq * I_q + (X_aq + X_rlq) * I_2q + (X_aq + X_rlq + X_1q) * I_1q
        #ψ_2q ~ -X_aq * I_q + (X_aq + X_rlq + X_2q) * I_2q + (X_aq + X_rlq) * I_1q

        #rotor voltage & current equations (67), (68), (70)
        I_fd ~ k_fd * I_d + (X_1d_loop * ψ_fd - (X_ad + X_rld) * ψ_1d) / X_det_d
        I_1d ~ k_1d * I_d + (X_fd_loop * ψ_1d - (X_ad + X_rld) * ψ_fd) / X_det_d
        I_1q ~ k_1q * I_q + ((X_2q_loop * ψ_1q - (X_aq + X_rlq) * ψ_2q) / X_det_q) * (1-salientpole) + salientpole * ψ_1q/X_1q
        I_2q ~ (k_2q * I_q + (X_1q_loop * ψ_2q - (X_aq + X_rlq) * ψ_1q) / X_det_q) * (1-salientpole)

        R_fd * I_fd + 1/ω_b * Dt(ψ_fd) ~ R_fd/X_ad * vf
        R_1d * I_1d + 1/ω_b * Dt(ψ_1d) ~ 0
        R_1q * I_1q + 1/ω_b * Dt(ψ_1q) ~ 0
        R_2q * I_2q + 1/ω_b * Dt(ψ_2q) ~ 0

        # inputs
        vf ~ vf_input ? vf_in.u : vf_set
        #τ_m ~ τ_m_input ? τ_m_in.u : τ_m_set
        #Alternative: (102)
        τ_m ~ pt/n - xmdm - dpu * n + addmt #xmdm Torque input signal; addmt additional torque parameter; dpu * n turbine shaft friction torque

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
