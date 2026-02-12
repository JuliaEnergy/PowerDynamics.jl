"""
    DynamicRLBranch(; R, L, ω0, r_src=1, r_dst=1)

Dynamic transmission line modeled as a series R-L circuit with optional transformer ratios.

This model represents a series resistance and inductance connecting two buses.
The line current is a differential state, suitable for:
- Transmission lines with significant inductive reactance
- Transformer models with leakage impedance
- Dynamic analysis requiring explicit current dynamics

# Parameters
- `R`: Line resistance [pu]
- `L`: Line reactance [pu] at frequency ω0. Related to physical inductance by L = ω0 * L_actual.
- `ω0`: Frame angular frequency [rad/s]. Default: 2π*50 rad/s.
- `r_src`: Transformer voltage ratio at source. Default: 1.
- `r_dst`: Transformer voltage ratio at destination. Default: 1.
"""
@mtkmodel DynamicRLBranch begin
    @parameters begin
        R, [description="Line resistance [pu]"]
        L, [description="Line reactance [pu] (frequency-normalized inductance)"]
        ω0=2π*50, [description="Frame angular frequency [rad/s]"]
        r_src = 1, [description="Transformer voltage ratio at source"]
        r_dst = 1, [description="Transformer voltage ratio at destination"]
    end
    @components begin
        src = Terminal()
        dst = Terminal()
    end
    @variables begin
        i_line_r(t), [guess=0, description="Series RL line current real part (dq frame) [pu]"]
        i_line_i(t), [guess=0, description="Series RL line current imaginary part (dq frame) [pu]"]
        i_mag(t), [description="Current magnitude [pu]"]
    end
    @equations begin
        # Series RL current dynamics in rotating dq frame
        Dt(i_line_r) ~ ω0 / L * (r_src*src.u_r - r_dst*dst.u_r) - R / L * ω0 * i_line_r + ω0 * i_line_i
        Dt(i_line_i) ~ ω0 / L * (r_src*src.u_i - r_dst*dst.u_i) - R / L * ω0 * i_line_i - ω0 * i_line_r

        # Terminal currents scaled by transformer ratios (for power conservation)
        dst.i_r ~ i_line_r * r_dst
        dst.i_i ~ i_line_i * r_dst
        src.i_r ~ -i_line_r * r_src
        src.i_i ~ -i_line_i * r_src

        i_mag ~ sqrt(dst.i_r^2 + dst.i_i^2)
    end
end
