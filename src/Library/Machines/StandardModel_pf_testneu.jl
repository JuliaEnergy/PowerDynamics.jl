@mtkmodel StandardModel_pf_testneu begin
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
        X_rld, [description="coupling reactance between field and damper winding"]
        X_rlq, [description="coupling reactance between q-axis damper windings"]
        #X_d, [description="d-axis synchronous reactance"]
        #X_q, [description="q-axis synchronous reactance"]
        #X′_d, [description="d-axis transient reactance"]
        #X′_q=0.0001, [description="q-axis transient reactance, not needed for salient pole machine"]
        X″_d, [description="d-axis subtransient reactance"]
        X″_q, [description="q-axis subtransient reactance"]
        X_ls, [description="stator leakage reactance"]
        #T′_d0, [description="d-axis transient time constant"]
        #T″_d0, [description="d-axis subtransient time constant"]
        #T′_q0=0, [description="q-axis transient time constant, not needed for salient pole machine"]
        #T″_q0, [description="q-axis subtransient time constant"]
        X_ad, [description="Mutual (magnetising) reactance, d-axis"]
        X_aq, [description="Mutual (magnetising) reactance, q-axis"]
        X_1q, [description="Hilfsvariable"]
        X_det_d, [description="Hilfsvariable"]
        X_det_q, [description="Hilfsvariable"]
        X_fd_loop, [description="Hilfsvariable"]
        X_1d_loop, [description="Hilfsvariable"]
        X_1q_loop, [description="Hilfsvariable"]
        X_2q_loop, [description="Hilfsvariable"]
        k_fd, [description="Hilfsvariable"]
        k_1d, [description="Hilfsvariable"]
        k_1q, [description="Hilfsvariable"]
        k_2q, [description="Hilfsvariable"]
        #X_fd, [description="Reactance of excitation (field) winding (d-axis)"]
        #X_1d, [description="Reactance of 1d-damper winding (d-axis)"]
        #X_1q, [description="Reactance of 1q-damper winding (q-axis)"]
        #X_2q, [description="Reactance of 2q-damper winding (q-axis)"]
        R_fd, [description="Resistance of excitation winding (d-axis)"]
        R_1d, [description="Resistance of 1d-damper winding (d-axis)"]
        R_1q, [description="Resistance of 1q-damper winding (q-axis)"]
        R_2q, [description="Resistance of 2q-damper winding (q-axis)"]
        H, [description="inertia constant"]
        # System and machine base
        S_b, [description="System power basis in MVA"]
        V_b, [description="System voltage basis in kV"]
        ω_b, [description="System base frequency in rad/s"] #(106)
        Sn, [description="Machine power rating in MVA"]
        Vn, [description="Machine voltage rating in kV"]
        cosn, [description="rated power factor - ??"] #oder ist das Variable?
        n_ref=1, [description="nominal speed (1 p.u.) of the machine or the speed of the local reference machine"]
        dkd, [description="Damping torque coefficient"]
        dpe, [description="Damping torque coefficient based on power"]
        D, [description="Damping constant"]
        salientpole=1, [description="salient pole or round-rotor machine"]
        xmdm=0, [description="Torque Input input signal in p.u."]
        addmt=0, [description="Additional Torque parameter in p.u."]
        pt, [description="Turbine Power input signal in p.u."]
        dpu=0, [description="dpu * n is turbine shaft friction torque in p.u.;"]
        # input/parameter switches
        if !vf_input
            vf_set, [guess=1, bounds=(0,Inf), description="field voltage"]
        end
        if !τ_m_input
            τ_m_set, [guess=1, bounds=(0,Inf), description="mechanical torque"]
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
    end
    begin
        T_to_loc(α)  = [ sin(α) -cos(α);
                         cos(α)  sin(α)]
        T_to_glob(α) = [ sin(α)  cos(α);
                        -cos(α)  sin(α)]
    end
    @equations begin
        # Park's transformations
        [terminal.u_r, terminal.u_i] .~ T_to_glob(δ)*[V_d, V_q] * Vn/V_b
        [I_d, I_q] .~ T_to_loc(δ)*[terminal.i_r, terminal.i_i] * Ibase(S_b, V_b)/Ibase(Sn, Vn)


        #stator flux equations (55), (60) ((59) wird indirekt über ks und Definitionen von ψ_1d abgedeckt?, (56) und (57) sind nur andere Darstellungsform)
        ψ_d ~ -(X_ls + X_ad) * I_d + X_ad * I_fd + X_ad * I_1d
        ψ_q ~ -(X_ls + X_aq) * I_q + X_aq * I_2q + X_aq * I_1q
        ψ_d ~ -X″_d * I_d + ψ″_d
        ψ_q ~ -X″_q * I_q + ψ″_q
        #ψ″_d ~ k_fd * ψ_fd + k_1d * ψ_1d #(59)
        #ψ″_q ~ k_1q * ψ_1q + k_2q * ψ_2q


        #stator voltage equations (72) und (73) (RMS Modell, anstatt (61) und (62) bzw (54))
        V_d ~ V″_d - R_s * I_d + n * X″_q * I_q
        V_q ~ V″_q - R_s * I_q - n * X″_d * I_d
        V″_d ~ - n * ψ″_q
        V″_q ~ n * ψ″_d

        #electrical torque (66)
        τ_e * cosn ~ I_q * ψ_d - I_d * ψ_q #cosn?

        #mechanical equation motor (103), (104), (105)
        τ_dkd ~ dkd * (n - n_ref)
        τ_dpe ~ dpe/n * (n - n_ref)
        τ_ag ~ 2 * H * 100 / (S_b * cosn)   #τ_ag und H werden dann hier auf Pgn bezogen - ist das richtig??
        Dt(n) ~ (τ_m - τ_e - τ_dkd - τ_dpe) / τ_ag #(100), (101) wird gar nicht gebraucht (t_base?)

        #Dt(ϕ) ~ ω_b * (n - ω_b/(2*π)) #(115) wenn δ = ϕ -> stimmt nicht. ϕ ist rotor Position im Vergleich zur Referenz-Spannung des Netzes. Ich brauche aber firel, also Winkel zwischen Refernzmaschine d-Achse und Generator d-Achse
        #δ ~ ϕ + π/2 #- phiu #phiu is the voltage angle of the machine terminal m:phiu (scheint 0 zu sein #(112), passt das mit den Achsen überhaupt?; δ hier in rad; fipol ist von Generator terminal zu q-Achse, δ in Milano zur d-Achse, und da sind d- und q-Achse auch vertauscht
        Dt(δ) ~ ω_b * (n - 1)

        #rotor flux linkage (58), überflüssig
        #ψ_fd ~ -X_ad * I_d + (X_ad + X_rld + X_fd) * I_fd + (X_ad + X_rld) * I_1d
        #ψ_1d ~ -X_ad * I_d + (X_ad + X_rld) * I_fd + (X_ad + X_rld + X_1d) * I_1d
        #ψ_1q ~ -X_aq * I_q + (X_aq + X_rlq) * I_2q + (X_aq + X_rlq + X_1q) * I_1q
        #ψ_2q ~ -X_aq * I_q + (X_aq + X_rlq + X_2q) * I_2q + (X_aq + X_rlq) * I_1q

        #rotor voltage & current equations (67) + (68), und (70)
        I_fd ~ k_fd * I_d + (X_1d_loop * ψ_fd - (X_ad + X_rld) * ψ_1d) / X_det_d #X_rl = X_rld ?? -> in Doku nur X_rl
        I_1d ~ k_1d * I_d + (X_fd_loop * ψ_1d - (X_ad + X_rld) * ψ_fd) / X_det_d
        I_1q ~ k_1q * I_q + ((X_2q_loop * ψ_1q - (X_aq + X_rlq) * ψ_2q) / X_det_q) * (1-salientpole) + salientpole * ψ_1q/X_1q #k_1q * I_q + ψ_1q/X_1q für salient pole; k_1q * I_q + (X_2q_loop * ψ_1q - (X_aq + X_rlq) * ψ_2q) / X_det_q round rotor
        I_2q ~ (k_2q * I_q + (X_1q_loop * ψ_2q - (X_aq + X_rlq) * ψ_1q) / X_det_q) * (1-salientpole) #0 für salient pole; k_2q * I_q + (X_1q_loop * ψ_2q - (X_aq + X_rlq) * ψ_1q) / X_det_q  round rotor

        R_fd * I_fd + 1/ω_b * Dt(ψ_fd) ~ vf
        R_1d * I_1d + 1/ω_b * Dt(ψ_1d) ~ 0
        R_1q * I_1q + 1/ω_b * Dt(ψ_1q) ~ 0
        R_2q * I_2q + 1/ω_b * Dt(ψ_2q) ~ 0

        # inputs
        vf ~ vf_input ? vf_in.u : vf_set
        τ_m ~ τ_m_input ? τ_m_in.u : τ_m_set
        #Alternativ: (102)
        #τ_m ~ pt/n - xmdm - dpu * n + addmt #xmdm Torque input signal; addmt additional torque parameter; dpu * n turbine shaft friction torque

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
