
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
sg1_bus, bus1, loop1 = let
    @named sm = SyncMachineStatorDynamics(J=3.5, D=1, wL=0.05, R=0.01, w0=w0)
    @named sg1_bus = compile_bus(MTKBus(sm); current_source=true)
    set_pfmodel!(sg1_bus, pfSlack(V=1, δ=0; current_source=true, assume_io_coupling=true))

    @named shunt = DynamicRCShunt(R=1/0.6, C=1e-5, ω0=w0)
    @named bus1 = compile_bus(MTKBus(shunt))
    set_pfmodel!(bus1, pfShunt(G=0.6, B=1e-5))

    loop1 = LoopbackConnection(; potential=[:u_r, :u_i], flow=[:i_r, :i_i], src=:sg1_bus, dst=:bus1)

    sg1_bus, bus1, loop1
end

## ===== Bus 2: SG Type 0 (PV for PF) =====
sg2_bus, bus2, loop2 = let
    @named sm = SyncMachineStatorDynamics(J=3.5, D=1, wL=0.05, R=0.01, w0=w0)
    @named sg2_bus = compile_bus(MTKBus(sm); current_source=true)
    set_pfmodel!(sg2_bus, pfPV(P=0.5, V=1; current_source=true, assume_io_coupling=true))

    @named shunt = DynamicRCShunt(R=1/0.6, C=1e-5, ω0=w0)
    @named bus2 = compile_bus(MTKBus(shunt))
    set_pfmodel!(bus2, pfShunt(G=0.6, B=1e-5))

    loop2 = LoopbackConnection(; potential=[:u_r, :u_i], flow=[:i_r, :i_i], src=:sg2_bus, dst=:bus2)

    sg2_bus, bus2, loop2
end

## ===== Bus 3: GFM Type 20 (PV for PF) =====
gfm_bus, bus3, loop3 = let
    @named droop = ComposableInverter.DroopInverter(filter_type=:LCL)
    @named gfm_bus = compile_bus(MTKBus(droop); current_source=true)

    # Parameters from JSON
    xwLf=0.05; Rf=0.01; xwCf=0.02; xwLc=0.01; Rc=0.002
    Xov=0.01; xDw=0.05; xfdroop=5; xfvdq=300; xfidq=600

    Lf = xwLf/w0; Cf = xwCf/w0
    Dw = xDw*w0; wf = xfdroop*2*pi
    w_v_odq = xfvdq*2*pi; w_i_ldq = xfidq*2*pi
    kp_i_ldq = w_i_ldq*Lf; ki_i_ldq = w_i_ldq^2*Lf/4
    kp_v_odq = w_v_odq*Cf; ki_v_odq = w_v_odq^2*Cf/4*50

    set_default!(gfm_bus, :droop₊droop₊Qset, 0)
    set_default!(gfm_bus, :droop₊droop₊Kp, Dw)
    set_default!(gfm_bus, :droop₊droop₊ω0, w0)
    set_default!(gfm_bus, :droop₊droop₊Kq, 0)
    set_default!(gfm_bus, :droop₊droop₊τ_q, Inf)
    set_default!(gfm_bus, :droop₊droop₊τ_p, 1/wf)
    set_default!(gfm_bus, :droop₊vsrc₊CC1_F, 0)
    set_default!(gfm_bus, :droop₊vsrc₊CC1_KI, ki_i_ldq)
    set_default!(gfm_bus, :droop₊vsrc₊CC1_KP, kp_i_ldq)
    set_default!(gfm_bus, :droop₊vsrc₊CC1_Fcoupl, 0)
    set_default!(gfm_bus, :droop₊vsrc₊VC_KP, kp_v_odq)
    set_default!(gfm_bus, :droop₊vsrc₊VC_KI, ki_v_odq)
    set_default!(gfm_bus, :droop₊vsrc₊VC_F, 0)
    set_default!(gfm_bus, :droop₊vsrc₊VC_Fcoupl, 0)
    set_default!(gfm_bus, :droop₊vsrc₊X_virt, Xov)
    set_default!(gfm_bus, :droop₊vsrc₊R_virt, 0)
    set_default!(gfm_bus, :droop₊vsrc₊X_g, xwLc)
    set_default!(gfm_bus, :droop₊vsrc₊B_c, xwCf)
    set_default!(gfm_bus, :droop₊vsrc₊Rf, Rf)
    set_default!(gfm_bus, :droop₊vsrc₊X_f, xwLf)
    set_default!(gfm_bus, :droop₊vsrc₊ω0, w0)
    set_default!(gfm_bus, :droop₊vsrc₊Rg, Rc)

    set_pfmodel!(gfm_bus, pfPV(P=0.5, V=1; current_source=true, assume_io_coupling=true))

    @named shunt = DynamicRCShunt(R=1/0.75, C=1e-5, ω0=w0)
    @named bus3 = compile_bus(MTKBus(shunt))
    set_pfmodel!(bus3, pfShunt(G=0.75, B=1e-5))

    loop3 = LoopbackConnection(; potential=[:u_r, :u_i], flow=[:i_r, :i_i], src=:gfm_bus, dst=:bus3)

    gfm_bus, bus3, loop3
end

## ===== Bus 4: GFL Type 10 with DC-link (PQ for PF) =====
gfl_bus, bus4, loop4 = let
    @named gfl = ComposableInverter.SimpleGFLDC()
    @named gfl_bus = compile_bus(MTKBus(gfl); current_source=true)

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

    set_default!(gfl_bus, :gfl₊ω0, w0)
    set_default!(gfl_bus, :gfl₊X_f, xwLf)
    set_default!(gfl_bus, :gfl₊Rf, Rf)
    set_default!(gfl_bus, :gfl₊PLL_Kp, kp_pll)
    set_default!(gfl_bus, :gfl₊PLL_Ki, ki_pll)
    set_default!(gfl_bus, :gfl₊PLL_τ_lpf, tau_pll)
    set_default!(gfl_bus, :gfl₊CC1_KP, kp_i_dq)
    set_default!(gfl_bus, :gfl₊CC1_KI, ki_i_dq)
    set_default!(gfl_bus, :gfl₊CC1_F, 0)
    set_default!(gfl_bus, :gfl₊CC1_Fcoupl, 0)
    set_default!(gfl_bus, :gfl₊C_dc, C_dc)
    set_default!(gfl_bus, :gfl₊V_dc, V_dc)
    set_default!(gfl_bus, :gfl₊kp_v_dc, kp_v_dc)
    set_default!(gfl_bus, :gfl₊ki_v_dc, ki_v_dc)

    set_pfmodel!(gfl_bus, pfPQ(P=0.5, Q=-0.2; current_source=true))

    @named shunt = DynamicRCShunt(R=1/0.05, C=1e-5, ω0=w0)
    @named bus4 = compile_bus(MTKBus(shunt))
    set_pfmodel!(bus4, pfShunt(G=0.05, B=1e-5))

    loop4 = LoopbackConnection(; potential=[:u_r, :u_i], flow=[:i_r, :i_i], src=:gfl_bus, dst=:bus4)

    gfl_bus, bus4, loop4
end

## ===== Lines =====
line12 = let
    @named branch = DynamicRLBranch(R=0.01, L=0.3, ω0=w0)
    lm = compile_line(MTKLine(branch); name=:l12, src=:bus1, dst=:bus2)
    @named branch_pf = PiLine(R=0.01, X=0.3)
    pfmod = compile_line(MTKLine(branch_pf); name=:l12_pfmod)
    set_pfmodel!(lm, pfmod)
    lm
end

line23 = let
    @named branch = DynamicRLBranch(R=0.01, L=0.3, ω0=w0)
    lm = compile_line(MTKLine(branch); name=:l23, src=:bus2, dst=:bus3)
    @named branch_pf = PiLine(R=0.01, X=0.3)
    pfmod = compile_line(MTKLine(branch_pf); name=:l23_pfmod)
    set_pfmodel!(lm, pfmod)
    lm
end

line31 = let
    @named branch = DynamicRLBranch(R=0.01, L=0.3, ω0=w0)
    lm = compile_line(MTKLine(branch); name=:l31, src=:bus3, dst=:bus1)
    @named branch_pf = PiLine(R=0.01, X=0.3)
    pfmod = compile_line(MTKLine(branch_pf); name=:l31_pfmod)
    set_pfmodel!(lm, pfmod)
    lm
end

line34 = let
    @named branch = DynamicRLBranch(R=0.01, L=0.3, ω0=w0, r_dst=0.99)
    lm = compile_line(MTKLine(branch); name=:l34, src=:bus3, dst=:bus4)
    @named branch_pf = PiLine(R=0.01, X=0.3, r_dst=0.99)
    pfmod = compile_line(MTKLine(branch_pf); name=:l34_pfmod)
    set_pfmodel!(lm, pfmod)
    lm
end


## ===== Network =====
nw = Network([sg1_bus, bus1, sg2_bus, bus2, gfm_bus, bus3, gfl_bus, bus4],
             [loop1, loop2, loop3, loop4, line12, line23, line31, line34]; warn_order=false)
# pfs0 = NWState(powerflow_model(nw))
# uflat(pfs0) .= -0.1

pfnw = powerflow_model(nw)
pfs0 = NWState(pfnw)

pfs = solve_powerflow(nw; abstol=1e-10, reltol=1e-10)
show_powerflow(pfs)
s0 = initialize_from_pf!(nw; pfs, subverbose=true, tol=1e-7, nwtol=1e-7)
break


## ===== Eigenvalue analysis =====
eigenvalues = jacobian_eigenvals(nw, s0) ./ (2 * pi)
println("4-Bus Eigenvalues:")
display(eigenvalues)

let
    fig = Figure(size=(600,500))
    ax1 = Axis(fig[1, 1], xlabel="Real Part [Hz]", ylabel="Imaginary Part [Hz]", title="Global Pole Map")
    scatter!(ax1, real.(eigenvalues), imag.(eigenvalues), marker=:xcross)
    ax2 = Axis(fig[1, 2], xlabel="Real Part [Hz]", ylabel="Imaginary Part [Hz]", title="Zoomed In")
    scatter!(ax2, real.(eigenvalues), imag.(eigenvalues), marker=:xcross)
    xlims!(ax2, -80, 20); ylims!(ax2, -150, 150)
    fig
end


####
#### Transfer Function Analysis
####

function bode_plot(Gs, title="", labels=["Bus $i" for i in 1:length(Gs)])
    with_theme(Theme(palette = (; color = Makie.wong_colors()[[1, 6, 2, 4]]))) do
        fig = Figure(; size=(800, 600))
        Label(fig[1, 1], title*"Bode Plot", fontsize=16, halign=:center, tellwidth=false)
        ax1 = Axis(fig[2, 1], xlabel="Frequency (rad/s)", ylabel="Gain (dB)", xscale=log10)
        ax2 = Axis(fig[3, 1], xlabel="Frequency (rad/s)", ylabel="Phase (deg)", xscale=log10)

        for (G, label) in zip(Gs, labels)
            fs = 10 .^ (range(log10(1e-1), log10(1e4); length=1000))
            ωs = 2π * fs
            ss = im .* ωs
            output, input = 1, 1
            gains = map(s -> 20 * log10(abs(G(s)[output, input])), ss)
            phases = rad2deg.(unwprap_rad(map(s -> angle(G(s)[output, input]), ss)))

            lines!(ax1, fs, gains; label, linewidth=2)
            lines!(ax2, fs, simple_unwrap(phases); label, linewidth=2)
        end
        axislegend(ax1)
        fig
    end
end

Gs = map([:sg1_bus, :sg2_bus, :gfm_bus, :gfl_bus]) do COMP
    vs = VIndex(COMP, [:busbar₊u_r, :busbar₊u_i])
    cs = VIndex(COMP, [:busbar₊i_r, :busbar₊i_i])
    G = NetworkDynamics.linearize_network(nw, s0; in=vs, out=cs).G
    # s -> -G(s)
end
bode_plot(Gs, "Y_dd ")
