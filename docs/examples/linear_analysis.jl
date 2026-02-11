
#
# This model implements the SimplusGT SynchronousMachine Type 0 model:
# "constant field flux, rotor motion of torque"
#
# This implementation stays as close as possible to the MATLAB original,
# using the same parameter names, equations, and d-q reference frame convention.

using PowerDynamics
using PowerDynamics.Library
using PowerDynamics.Library.ComposableInverter
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using OrdinaryDiffEqRosenbrock
using NetworkDynamics
using OrdinaryDiffEqNonlinearSolve
using CairoMakie

####
#### Models
####


@mtkmodel SyncMachineStatorDynamics begin
    @components begin
        terminal = Terminal()
    end

    @parameters begin
        J, [description="Inertia constant [MWs²/MVA]"]
        D, [description="Damping coefficient [pu]"]
        wL, [description="Stator inductance * base frequency [pu]"]
        R, [description="Stator resistance [pu]"]
        w0, [description="Base frequency [rad/s]"]
        psi_f, [guess=1, description="Field flux linkage [pu]"]
        T_m, [guess=1, description="Mechanical torque [pu]"]
    end

    @variables begin
        # State variables
        i_d(t), [guess=0, description="d-axis stator current [pu]"]
        i_q(t), [guess=1, description="q-axis stator current [pu]"]
        w(t), [guess=2*pi*50, description="Rotor speed deviation [rad/s]"]
        theta(t), [guess=0, description="Rotor angle [rad]"]
        # Algebraic variables
        v_d(t), [guess=1, description="d-axis terminal voltage [pu]"]
        v_q(t), [guess=0, description="q-axis terminal voltage [pu]"]
        psi_d(t), [guess=1, description="d-axis flux linkage [pu]"]
        psi_q(t), [guess=0, description="q-axis flux linkage [pu]"]
        Te(t), [guess=1, description="Electrical torque [pu]"]
    end

    begin
        # Parameter scaling (MATLAB lines 56-59)
        J_pu = J*2/w0^2
        D_pu = D/w0^2
        L = wL/w0

        T_to_loc(α)  = [ cos(α) sin(α);
                         -sin(α)  cos(α)]
        T_to_glob(α) = T_to_loc(-α)
    end

    @equations begin
        # Coordinate transformations (q-axis aligned with field flux, MATLAB convention)
        [terminal.i_r, terminal.i_i] .~ -T_to_glob(theta)*[i_d, i_q]
        [v_d, v_q] .~ T_to_loc(theta)*[terminal.u_r, terminal.u_i]

        # Flux linkages (MATLAB lines 72-73)
        psi_d ~ L*i_d
        psi_q ~ L*i_q - psi_f

        # Electrical torque (MATLAB line 74)
        Te ~ psi_f * i_d

        # State equations - Type 0 (MATLAB lines 75-78)
        Dt(i_d) ~ (v_d - R*i_d + w*psi_q)/L
        Dt(i_q) ~ (v_q - R*i_q - w*psi_d)/L
        Dt(w) ~ (Te - T_m - D_pu*w)/J_pu
        Dt(theta) ~ w - w0
    end
end

@mtkmodel DynRLLine begin
    @parameters begin
        R
        L
        ωbase
        r_src = 1, [description="Transformer ratio at source"]
        r_dst = 1, [description="Transformer ratio at destination"]
    end
    @components begin
        src = Terminal()
        dst = Terminal()
    end
    @variables begin
        i_line_r(t), [guess=0, description="Series RL current (real)"]
        i_line_i(t), [guess=0, description="Series RL current (imag)"]
        i_mag(t)
    end
    @equations begin
        # Series RL current dynamics
        Dt(i_line_r) ~ ωbase / L * (r_src*src.u_r - r_dst*dst.u_r) - R / L * ωbase * i_line_r + ωbase * i_line_i
        Dt(i_line_i) ~ ωbase / L * (r_src*src.u_i - r_dst*dst.u_i) - R / L * ωbase * i_line_i - ωbase * i_line_r

        # Terminal currents scaled by transformer ratios
        dst.i_r ~ i_line_r * r_dst
        dst.i_i ~ i_line_i * r_dst
        src.i_r ~ -i_line_r * r_src
        src.i_i ~ -i_line_i * r_src

        i_mag ~ sqrt(dst.i_r^2 + dst.i_i^2)
    end
end

@mtkmodel ConstantRShunt begin
    @parameters begin
        R
    end
    @components begin
        terminal = Terminal()
    end
    @equations begin
        terminal.i_r ~ -terminal.u_r / R
        terminal.i_i ~ -terminal.u_i / R
    end
end

@mtkmodel DynamicShunt begin
    @components begin
        terminal = Terminal()
    end
    @parameters begin
        ω0=2π*50
        B
    end
    @variables begin
        V_C_r(t), [guess=1]
        V_C_i(t), [guess=0]
        i_C_r(t), [guess=0]
        i_C_i(t), [guess=0]
    end
    @equations begin
        B/ω0 * Dt(V_C_r) ~ -i_C_r + B*V_C_i
        B/ω0 * Dt(V_C_i) ~ -i_C_i - B*V_C_r
        # Terminal voltage = capacitor voltage
        terminal.u_r ~ V_C_r
        terminal.u_i ~ V_C_i
        # Grid current: i_C points towards grid, terminal.i uses sending convention
        terminal.i_r ~ i_C_r
        terminal.i_i ~ i_C_i
    end
end

####
####  4 Bus example
####
## Network from SimplusGT UserData.json:
## Bus 1: SG Type 0 (slack), Bus 2: SG Type 0 (PV),
## Bus 3: GFM Type 20 (PV), Bus 4: GFL Type 10 (PQ)
## Lines: 1→2, 2→3, 3→1 (R=0.01, wL=0.3), 3→4 (R=0.01, wL=0.3, turns_ratio=0.99)
## Self-branch shunts: Bus 1,2: G=0.6, wC=1e-5; Bus 3: G=0.75; Bus 4: G=0.05
using PowerDynamics.Library.ComposableInverter
w0 = 2π*50

## ===== Bus 1: SG Type 0 (Slack for PF) =====
sg1_bus = let
    @named sm = SyncMachineStatorDynamics(J=3.5, D=1, wL=0.05, R=0.01, w0=w0)
    @named rshunt = ConstantRShunt(R=1/0.6)
    @named cshunt = DynamicShunt(B=1e-5, ω0=w0)
    bus = compile_bus(MTKBus([sm, rshunt, cshunt]; name=:bus1))

    @named slack = Library.VδConstraint(V=1, δ=0)
    @named rshunt_pf = ConstantRShunt(R=1/0.6)
    pfmod = compile_bus(MTKBus([slack, rshunt_pf]); name=:bus1_pfmod)
    set_pfmodel!(bus, pfmod)
    bus
end

## ===== Bus 2: SG Type 0 (PV for PF) =====
sg2_bus = let
    @named sm = SyncMachineStatorDynamics(J=3.5, D=1, wL=0.05, R=0.01, w0=w0)
    @named rshunt = ConstantRShunt(R=1/0.6)
    @named cshunt = DynamicShunt(B=1e-5, ω0=w0)
    bus = compile_bus(MTKBus([sm, rshunt, cshunt]; name=:bus2))

    @named pv = Library.PVConstraint(P=0.5, V=1)
    @named rshunt_pf = ConstantRShunt(R=1/0.6)
    pfmod = compile_bus(MTKBus([pv, rshunt_pf]); name=:bus2_pfmod)
    set_pfmodel!(bus, pfmod)
    set_voltage!(bus; mag=1, arg=0)
    bus
end

## ===== Bus 3: GFM Type 20 (PV for PF) =====
gfm_bus = let
    @named droop = ComposableInverter.DroopInverter(filter_type=:LCL)
    @named rshunt = ConstantRShunt(R=1/0.75)
    @named cshunt = DynamicShunt(B=1e-5, ω0=w0)
    bus = compile_bus(MTKBus([droop, rshunt, cshunt]; name=:bus3))

    # Parameters from JSON
    xwLf=0.05; Rf=0.01; xwCf=0.02; xwLc=0.01; Rc=0.002
    Xov=0.01; xDw=0.05; xfdroop=5; xfvdq=300; xfidq=600

    Lf = xwLf/w0; Cf = xwCf/w0
    Dw = xDw*w0; wf = xfdroop*2*pi
    w_v_odq = xfvdq*2*pi; w_i_ldq = xfidq*2*pi
    kp_i_ldq = w_i_ldq*Lf; ki_i_ldq = w_i_ldq^2*Lf/4
    kp_v_odq = w_v_odq*Cf; ki_v_odq = w_v_odq^2*Cf/4*50

    set_default!(bus, :droop₊droop₊Kp, Dw)
    set_default!(bus, :droop₊droop₊ω0, w0)
    set_default!(bus, :droop₊droop₊Kq, 0)
    set_default!(bus, :droop₊droop₊τ_q, Inf)
    set_default!(bus, :droop₊droop₊τ_p, 1/wf)
    set_default!(bus, :droop₊vsrc₊CC1_F, 0)
    set_default!(bus, :droop₊vsrc₊CC1_KI, ki_i_ldq)
    set_default!(bus, :droop₊vsrc₊CC1_KP, kp_i_ldq)
    set_default!(bus, :droop₊vsrc₊CC1_Fcoupl, 0)
    set_default!(bus, :droop₊vsrc₊VC_KP, kp_v_odq)
    set_default!(bus, :droop₊vsrc₊VC_KI, ki_v_odq)
    set_default!(bus, :droop₊vsrc₊VC_F, 0)
    set_default!(bus, :droop₊vsrc₊VC_Fcoupl, 0)
    set_default!(bus, :droop₊vsrc₊X_virt, Xov)
    set_default!(bus, :droop₊vsrc₊R_virt, 0)
    set_default!(bus, :droop₊vsrc₊X_g, xwLc)
    set_default!(bus, :droop₊vsrc₊B_c, xwCf)
    set_default!(bus, :droop₊vsrc₊Rf, Rf)
    set_default!(bus, :droop₊vsrc₊X_f, xwLf)
    set_default!(bus, :droop₊vsrc₊ω0, w0)
    set_default!(bus, :droop₊vsrc₊Rg, Rc)

    @named pv = Library.PVConstraint(P=0.5, V=1)
    @named rshunt_pf = ConstantRShunt(R=1/0.75)
    pfmod = compile_bus(MTKBus([pv, rshunt_pf]); name=:bus3_pfmod)
    set_pfmodel!(bus, pfmod)
    set_voltage!(bus; mag=1, arg=0)
    bus
end

## ===== Bus 4: GFL Type 10 with DC-link (PQ for PF) =====
gfl_bus = let
    @named gfl = ComposableInverter.SimpleGFLDC()
    @named rshunt = ConstantRShunt(R=1/0.05)
    @named cshunt = DynamicShunt(B=1e-5, ω0=w0)
    bus = compile_bus(MTKBus([gfl, rshunt, cshunt]; name=:bus4))

    # Parameters from JSON
    V_dc=2.5; C_dc=1.25; f_v_dc=5
    xwLf=0.03; Rf=0.01
    f_pll=5; f_tau_pll=300; f_i_dq=600

    Lf = xwLf/w0
    w_vdc = f_v_dc*2*pi
    kp_v_dc = V_dc*C_dc*w_vdc
    ki_v_dc = kp_v_dc*w_vdc/4

    w_pll = f_pll*2*pi
    kp_pll = w_pll
    ki_pll = w_pll^2/4
    tau_pll = 1/(f_tau_pll*2*pi)

    w_i_dq = f_i_dq*2*pi
    kp_i_dq = Lf * w_i_dq
    ki_i_dq = Lf * w_i_dq^2 / 4

    set_default!(bus, :gfl₊ω0, w0)
    set_default!(bus, :gfl₊X_f, xwLf)
    set_default!(bus, :gfl₊Rf, Rf)
    set_default!(bus, :gfl₊PLL_Kp, kp_pll)
    set_default!(bus, :gfl₊PLL_Ki, ki_pll)
    set_default!(bus, :gfl₊PLL_τ_lpf, tau_pll)
    set_default!(bus, :gfl₊CC1_KP, kp_i_dq)
    set_default!(bus, :gfl₊CC1_KI, ki_i_dq)
    set_default!(bus, :gfl₊CC1_F, 0)
    set_default!(bus, :gfl₊CC1_Fcoupl, 0)
    set_default!(bus, :gfl₊C_dc, C_dc)
    set_default!(bus, :gfl₊V_dc, V_dc)
    set_default!(bus, :gfl₊kp_v_dc, kp_v_dc)
    set_default!(bus, :gfl₊ki_v_dc, ki_v_dc)

    @named pq = Library.PQConstraint(P=0.5, Q=-0.2)
    @named rshunt_pf = ConstantRShunt(R=1/0.05)
    pfmod = compile_bus(MTKBus([pq, rshunt_pf]); name=:bus4_pfmod)
    set_pfmodel!(bus, pfmod)
    set_voltage!(bus; mag=1, arg=0)
    bus
end

## ===== Lines =====
line12 = let
    @named branch = DynRLLine(R=0.01, L=0.3, ωbase=w0)
    lm = compile_line(MTKLine(branch); name=:l12, src=:bus1, dst=:bus2)
    @named branch_pf = PiLine(R=0.01, X=0.3)
    pfmod = compile_line(MTKLine(branch_pf); name=:l12_pfmod)
    set_pfmodel!(lm, pfmod)
    lm
end

line23 = let
    @named branch = DynRLLine(R=0.01, L=0.3, ωbase=w0)
    lm = compile_line(MTKLine(branch); name=:l23, src=:bus2, dst=:bus3)
    @named branch_pf = PiLine(R=0.01, X=0.3)
    pfmod = compile_line(MTKLine(branch_pf); name=:l23_pfmod)
    set_pfmodel!(lm, pfmod)
    lm
end

line31 = let
    @named branch = DynRLLine(R=0.01, L=0.3, ωbase=w0)
    lm = compile_line(MTKLine(branch); name=:l31, src=:bus3, dst=:bus1)
    @named branch_pf = PiLine(R=0.01, X=0.3)
    pfmod = compile_line(MTKLine(branch_pf); name=:l31_pfmod)
    set_pfmodel!(lm, pfmod)
    lm
end

line34 = let
    @named branch = DynRLLine(R=0.01, L=0.3, ωbase=w0, r_dst=0.99)
    lm = compile_line(MTKLine(branch); name=:l34, src=:bus3, dst=:bus4)
    @named branch_pf = PiLine(R=0.01, X=0.3, r_dst=0.99)
    pfmod = compile_line(MTKLine(branch_pf); name=:l34_pfmod)
    set_pfmodel!(lm, pfmod)
    lm
end

## ===== Network =====
nw = Network([sg1_bus, sg2_bus, gfm_bus, gfl_bus], [line12, line23, line31, line34])
show_powerflow(solve_powerflow(nw))
s0 = initialize_from_pf!(nw; subverbose=true, tol=1e-7, nwtol=1e-7)

## ===== Eigenvalue analysis =====
eigenvalues = jacobian_eigenvals(nw, s0) ./ (2 * pi)
println("4-Bus Eigenvalues:")
display(eigenvalues)

let
    fig = Figure(size=(1200,400))
    ax1 = Axis(fig[1, 1], xlabel="Real Part [Hz]", ylabel="Imaginary Part [Hz]", title="All Eigenvalues")
    scatter!(ax1, real.(eigenvalues), imag.(eigenvalues), marker=:xcross)
    ax2 = Axis(fig[1, 2], xlabel="Real Part [Hz]", ylabel="Imaginary Part [Hz]", title="Zoomed")
    scatter!(ax2, real.(eigenvalues), imag.(eigenvalues), marker=:xcross)
    ylims!(ax2, -2550, 2550)
    ax3 = Axis(fig[1, 3], xlabel="Real Part [Hz]", ylabel="Imaginary Part [Hz]", title="Low-Frequency Modes")
    scatter!(ax3, real.(eigenvalues), imag.(eigenvalues), marker=:xcross)
    xlims!(ax3, -80, 20); ylims!(ax3, -150, 150)
    fig
end


####
#### Transfer Function Analysis
####

function bode_plot(G, title="")
    fs = 10 .^ (range(log10(1e-1), log10(1e4); length=1000))
    ωs = 2π * fs
    ss = im .* ωs

    fig = Figure(; size=(800, 600))

    output, input = 1, 1
    gains = map(s -> 20 * log10(abs(G(s)[output, input])), ss)
    phases = map(s -> angle(G(s)[output, input]) * 180 / pi, ss)
    Label(fig[1, 1], title*"Bode Plot of $input ↦ $output", fontsize=16, halign=:center, tellwidth=false)
    ax1 = Axis(fig[2, 1], xlabel="Frequency (rad/s)", ylabel="Gain (dB)", xscale=log10)
    lines!(ax1, fs, gains, color=:blue, label="Gain")
    ax2 = Axis(fig[3, 1], xlabel="Frequency (rad/s)", ylabel="Phase (deg)", xscale=log10)
    lines!(ax2, fs, phases, color=:red, label="Phase")
    fig
end

COMP = 1
vs = VIndex(COMP, [:busbar₊u_r, :busbar₊u_i])
cs = [VIndex(COMP, :busbar₊i_r), VIndex(COMP, :busbar₊i_i)]
G = NetworkDynamics.linearized_model(nw, s0; in=vs, out=cs)
bode_plot(s -> -G(s), "Component $i: u -> i Transfer function")
