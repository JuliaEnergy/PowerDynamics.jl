@mtkmodel SauerPaiMachine begin
    @components begin
        terminal=Terminal()
    end
    @parameters begin
        ω_s=2π*50, [description="synchronous speed"]
        R_s, [description="stator resistance"]
        X_d, [description="d-axis synchronous reactance"]
        X′_d, [description="d-axis transient reactance"]
        X_q, [description="q-axis synchronous reactance"]
        X′_q, [description="q-axis transient reactance"]
        X″_d, [description="d-axis subtransient reactance"]
        X″_q, [description="q-axis subtransient reactance"]
        X_ls, [description="stator leakage reactance"]
        T′_d0, [description="d-axis transient time constant"]
        T″_d0, [description="d-axis subtransient time constant"]
        T′_q0, [description="q-axis transient time constant"]
        T″_q0, [description="q-axis subtransient time constant"]
        H, [description="inertia constant"]
    end
    @variables begin
        ψ_d(t), [description="d-axis flux linkage"]
        ψ_q(t), [description="q-axis flux linkage"]
        # ψ_0(t), [description="zero sequence flux linkage"]
        ψ_1d(t), [description="flux linkage assosciated with X″_d"]
        ψ_2q(t), [description="flux linkage assosciated with X″_q"]
        I_d(t), [description="d-axis current"]
        I_q(t), [description="q-axis current"]
        # I_0(t), [description="zero sequence current"]
        V_d(t), [description="d-axis voltage"]
        V_q(t), [description="q-axis voltage"]
        # V_0(t), [description="zero sequence voltage"]
        E_fd(t), [description="field voltage"]
        E′_d(t), [description="transient voltage behind transient reactance in d-axis"]
        E′_q(t), [description="transient voltage behind transient reactance in q-axis"]
        δ(t), [description="rotor angle"]
        ω(t), [description="rotor speed"]
        T_m(t), [description="mechanical torque"]
        T_fw(t), [description="windage and friction torque"]
    end
    @equations begin
        1/ω_s * Dt(ψ_d) ~ R_s*I_d + ω/ω_s * ψ_q + V_d
        1/ω_s * Dt(ψ_q) ~ R_s*I_q - ω/ω_s * ψ_d + V_q
        # 1/ω_s * Dt(ψ_0) ~ R_s*I_0 + V_0
        T′_d0 * Dt(E′_q) ~ -E′_q - (X_d - X′_d)* (I_d - (X′_d-X″_d)/(X′_d-X_ls)^2 * (ψ_1d + (X′_d-X_ls)*I_d - E′_q)) + E_fd
        T″_d0 * Dt(ψ_1d) ~ -ψ_1d + E′_q - (X′_d - X_ls)*I_d
        T′_q0 * Dt(E′_d) ~ -E′_d - (X_q - X′_q)* (I_q - (X′_q-X″_q)/(X′_q-X_ls)^2 * (ψ_2q + (X′_q-X_ls)*I_q - E′_d))
        T″_q0 * Dt(ψ_2q) ~ -ψ_2q - E′_d - (X′_q - X_ls)*I_q
        Dt(δ) ~ ω - ω_s
        2*H/ω_s * Dt(ω) ~ T_m  - (ψ_d*I_q - ψ_q*I_d) - T_fw
        ψ_d ~ -X″_d*I_d + (X″_d - X_ls)/(X′_d - X_ls) * E′_q + (X′_d - X″_d)/(X′_d - X_ls) * ψ_1d
        ψ_q ~ -X″_q*I_q + (X″_q - X_ls)/(X′_q - X_ls) * E′_d + (X′_q - X″_q)/(X′_q - X_ls) * ψ_2q
        # ψ_0 ~ -X_ls*I_0+ V_0
    end
end
