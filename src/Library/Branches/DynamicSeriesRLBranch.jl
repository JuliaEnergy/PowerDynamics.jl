"""
    DynamicSeriesRLBranch(; R, L, r_src=1, r_dst=1)

Dynamic transmission line modeled as a series R-L circuit with optional transformer ratios.

This model represents a series resistance and inductance connecting two buses.
The line current is a differential state, suitable for:
- Transmission lines with significant inductive reactance
- Transformer models with leakage impedance
- Dynamic analysis requiring explicit current dynamics

# Parameters
- `R`: Line resistance [pu]
- `L`: Line reactance [pu] at frequency `ωbase`. Related to physical inductance by L = ωbase * L_actual.
- `r_src`: Transformer voltage ratio at source. Default: 1.
- `r_dst`: Transformer voltage ratio at destination. Default: 1.

$(PowerDynamics.ref_source_file(@__FILE__, @__LINE__))
"""
@mtkmodel DynamicSeriesRLBranch begin
    @parameters begin
        R, [description="Line resistance [pu]"]
        L, [description="Line reactance [pu] (frequency-normalized inductance)"]
        ωbase, [bound_to = :systembase₊ωbase, description="System frequency base [rad/s]"]
        ωframe, [bound_to = :systembase₊ωframe, description="Global dq frame speed [pu]"]
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
        # Series RL current dynamics in the rotating dq frame, in the classical
        # reactance-parameterized form: the only ωbase is the L/ωbase coefficient of the
        # derivative (unit), the cross terms carry the frame speed ωframe [pu].
        L/ωbase * Dt(i_line_r) ~ (r_src*src.u_r - r_dst*dst.u_r) - R*i_line_r + ωframe*L*i_line_i
        L/ωbase * Dt(i_line_i) ~ (r_src*src.u_i - r_dst*dst.u_i) - R*i_line_i - ωframe*L*i_line_r

        # Terminal currents scaled by transformer ratios (for power conservation)
        dst.i_r ~ i_line_r * r_dst
        dst.i_i ~ i_line_i * r_dst
        src.i_r ~ -i_line_r * r_src
        src.i_i ~ -i_line_i * r_src

        i_mag ~ sqrt(dst.i_r^2 + dst.i_i^2)
    end
end
