module ComposableInverter

using ModelingToolkit: ModelingToolkit, @named, simplify, t_nounits as t, D_nounits as Dt,
                       @component
# needed for @mtkmodel
using ModelingToolkit: @mtkmodel, @variables, @parameters, @unpack, Num, System, Equation, connect, setmetadata
using ModelingToolkitStandardLibrary.Blocks: RealInput, RealOutput

using ..PowerDynamics: Terminal

export LCFilter, LCLFilter, VC, CC1, CC2, SimplePLL, VoltageSource, CurrentSource
export DroopOuter, DroopInverter

@mtkmodel LCFilter begin
    @components begin
        terminal = Terminal()
    end
    @structural_parameters begin
        ω0
        Rf
        X_f  # Inductive reactance (pu Ω)
        B_c  # Capacitive susceptance (pu S)
    end
    @variables begin
        i_f_r(t), [guess=0]
        i_f_i(t), [guess=0]
        i_f_mag(t)
        V_C_r(t), [guess=1]
        V_C_i(t), [guess=0]
        V_C_mag(t)
        V_I_r(t)  # Input from controller
        V_I_i(t)
        i_g_r(t), [guess=0]  # Virtual grid current (alias for uniform interface)
        i_g_i(t), [guess=0]  # Virtual grid current (alias for uniform interface)
        i_g_mag(t)
    end
    @equations begin
        # Inductor: (X_f/ω0)*di_f/dt = V_I - V_C - Rf*i_f + X_f*W*i_f
        (X_f/ω0) * Dt(i_f_r) ~ V_I_r - V_C_r - Rf*i_f_r + X_f*i_f_i
        (X_f/ω0) * Dt(i_f_i) ~ V_I_i - V_C_i - Rf*i_f_i - X_f*i_f_r
        # Capacitor: (B_c/ω0)*dV_C/dt = i_f - i_g + B_c*W*V_C
        (B_c/ω0) * Dt(V_C_r) ~ i_f_r - i_g_r + B_c*V_C_i
        (B_c/ω0) * Dt(V_C_i) ~ i_f_i - i_g_i - B_c*V_C_r
        # Terminal voltage = capacitor voltage
        terminal.u_r ~ V_C_r
        terminal.u_i ~ V_C_i
        # Grid current: i_g points towards grid, terminal.i uses sending convention
        terminal.i_r ~ i_g_r
        terminal.i_i ~ i_g_i
        #observables
        i_f_mag ~ sqrt(i_f_r^2 + i_f_i^2)
        V_C_mag ~ sqrt(V_C_r^2 + V_C_i^2)
        i_g_mag ~ sqrt(i_g_r^2 + i_g_i^2)
    end
end

@mtkmodel LCLFilter begin
    @components begin
        terminal = Terminal()
    end
    @structural_parameters begin
        ω0
        Rf
        X_f  # Inverter-side inductive reactance (pu Ω)
        B_c  # Capacitive susceptance (pu S)
        Rg
        X_g  # Grid-side inductive reactance (pu Ω)
    end
    @parameters begin
        connected = 1
    end
    @variables begin
        i_f_r(t), [guess=0]
        i_f_i(t), [guess=0]
        i_f_mag(t)
        V_C_r(t), [guess=1]
        V_C_i(t), [guess=0]
        V_C_mag(t)
        i_g_r(t), [guess=0]
        i_g_i(t), [guess=0]
        i_g_mag(t)
        V_I_r(t)  # Input from controller
        V_I_i(t)
        P_g(t)
        Q_g(t)
    end
    @equations begin
        # Inverter-side inductor
        (X_f/ω0) * Dt(i_f_r) ~ V_I_r - V_C_r - Rf*i_f_r + X_f*i_f_i
        (X_f/ω0) * Dt(i_f_i) ~ V_I_i - V_C_i - Rf*i_f_i - X_f*i_f_r
        # Capacitor
        (B_c/ω0) * Dt(V_C_r) ~ i_f_r - i_g_r + B_c*V_C_i
        (B_c/ω0) * Dt(V_C_i) ~ i_f_i - i_g_i - B_c*V_C_r
        # Grid-side inductor
        (X_g/ω0) * Dt(i_g_r) ~ V_C_r - connected*terminal.u_r - Rg*i_g_r + X_g*i_g_i
        (X_g/ω0) * Dt(i_g_i) ~ V_C_i - connected*terminal.u_i - Rg*i_g_i - X_g*i_g_r
        # Current flows out through terminal
        terminal.i_r ~ connected*i_g_r
        terminal.i_i ~ connected*i_g_i
        #observables
        i_f_mag ~ sqrt(i_f_r^2 + i_f_i^2)
        V_C_mag ~ sqrt(V_C_r^2 + V_C_i^2)
        i_g_mag ~ sqrt(i_g_r^2 + i_g_i^2)
        # power between i_g and output voltage
        P_g ~ terminal.u_r*i_g_r + terminal.u_i*i_g_i
        Q_g ~ terminal.u_i*i_g_r - terminal.u_r*i_g_i
    end
end

@component function CC1(; name, X_f, F, Fcoupl=1, KP, KI)
    vars = @variables begin
        γ_d(t), [guess=0]
        γ_q(t), [guess=0]
        V_I_d(t)
        V_I_q(t)
        i_f_d(t)
        i_f_q(t)
        i_f_ref_d(t)
        i_f_ref_q(t)
        V_C_d(t)
        V_C_q(t)
    end

    γ = [γ_d, γ_q]
    i_f = [i_f_d, i_f_q]
    i_f_ref = [i_f_ref_d, i_f_ref_q]
    V_I = [V_I_d, V_I_q]
    V_C = [V_C_d, V_C_q]

    W = [0 1; -1 0]

    eqs = vcat(
        Dt.(γ) .~ i_f_ref - i_f,
        V_I .~ -Fcoupl*X_f*W*i_f + KP*(i_f_ref - i_f) + KI*γ + F*V_C
    )
    System(eqs, t; name)
end

@component function VC(; name, B_c, F, Fcoupl=1, KP, KI)
    vars = @variables begin
        γ_d(t), [guess=0]
        γ_q(t), [guess=0]
        i_f_ref_d(t)
        i_f_ref_q(t)
        V_C_d(t)
        V_C_q(t)
        V_C_ref_d(t)
        V_C_ref_q(t)
        i_g_d(t)
        i_g_q(t)
    end

    γ = [γ_d, γ_q]
    V_C = [V_C_d, V_C_q]
    i_g = [i_g_d, i_g_q]
    V_C_ref = [V_C_ref_d, V_C_ref_q]
    i_f_ref = [i_f_ref_d, i_f_ref_q]

    W = [0 1; -1 0]

    eqs = vcat(
        Dt.(γ) .~ V_C_ref - V_C,
        i_f_ref .~ -Fcoupl*B_c*W*V_C + KP*(V_C_ref - V_C) + KI*γ + F*i_g
    )
    System(eqs, t; name)
end

@component function SimplePLL(; name, Kp, Ki)
    vars = @variables begin
        δ_pll(t), [guess=0]    # PLL angle output
        ω_pll(t), [guess=0]    # PLL frequency output
        u_q(t)                  # q-axis voltage (internal)
        u_r(t)                  # measured voltage (real)
        u_i(t)                  # measured voltage (imag)
    end

    # u_q drives to zero when δ_pll tracks grid angle
    # u_q = (u_r*sin(-δ) + u_i*cos(-δ)) / |u|
    u_mag = sqrt(u_r^2 + u_i^2)

    eqs = [
        u_q ~ (u_r*sin(-δ_pll) + u_i*cos(-δ_pll)) / u_mag
        Dt(δ_pll) ~ ω_pll + Kp*u_q
        Dt(ω_pll) ~ Ki*u_q
    ]
    System(eqs, t; name)
end

@component function CC2(; name, X_g, F, Fcoupl=1, KP, KI)
    vars = @variables begin
        γ_d(t), [guess=0]
        γ_q(t), [guess=0]
        V_C_ref_d(t)
        V_C_ref_q(t)
        i_g_d(t)
        i_g_q(t)
        i_g_ref_d(t)
        i_g_ref_q(t)
        V_g_d(t)
        V_g_q(t)
    end

    γ = [γ_d, γ_q]
    i_g = [i_g_d, i_g_q]
    i_g_ref = [i_g_ref_d, i_g_ref_q]
    V_C_ref = [V_C_ref_d, V_C_ref_q]
    V_g = [V_g_d, V_g_q]

    W = [0 1; -1 0]

    eqs = vcat(
        Dt.(γ) .~ i_g_ref - i_g,
        V_C_ref .~ -Fcoupl*X_g*W*i_g + KP*(i_g_ref - i_g) + KI*γ + F*V_g
    )
    System(eqs, t; name)
end


@component function VoltageSource(; name, Vset_input=false, filter_type=:LC)
    @parameters begin
        # Filter parameters (actual electrical components in per-unit)
        Rf  = 0.01
        X_f = 0.007, [description="Inverter-side reactance (pu Ω)"]
        B_c = 0.5,  [description="Capacitive susceptance (pu S)"]
        ω0  = 2π*50
        # cc1
        CC1_KP = 0.063
        CC1_KI = 63
        CC1_F  = 0
        CC1_Fcoupl = 1
        # vc
        VC_KP = 0.952
        VC_KI = 317.4
        VC_F = 0
        VC_Fcoupl = 1
        # virtual inductance
        R_virt = 0, [description="Virtual resistance (feedforward of i_g) (pu Ω)"]
        X_virt = 0, [description="Virtual inductance (feedforward of i_g) for feedforward (pu Ω)"]
    end
    # Grid-side inductor parameters (only for LCL)
    if filter_type == :LCL
        @parameters begin
            Rg  = 0.01
            X_g = 0.007, [description="Grid-side reactance (pu Ω)"]
        end
    end
    if !Vset_input
        @parameters begin
            Vset, [guess=1, description="Voltage magnitude setpoint in pu"]
            δset, [guess=0, description="dq-frame position relative to global frame"]
        end
    end

    # Create filter subsystem based on filter_type
    if filter_type == :LC
        @named filter = LCFilter(; ω0, Rf, X_f, B_c)
    else  # :LCL
        @named filter = LCLFilter(; ω0, Rf, X_f, B_c, Rg, X_g)
    end

    @named cc1 = CC1(; X_f=X_f, F=CC1_F, Fcoupl=CC1_Fcoupl, KP=CC1_KP, KI=CC1_KI)
    systems = @named begin
        terminal = Terminal()
        vc = VC(; B_c=B_c, F=VC_F, Fcoupl=VC_Fcoupl, KP=VC_KP, KI=VC_KI)
    end
    push!(systems, filter)
    push!(systems, cc1)
    if Vset_input
        @named Vset_in = RealInput()
        @named δset_in = RealInput()
        append!(systems, [Vset_in, δset_in])
    end
    _Vset = Vset_input ? Vset_in.u : Vset
    _δ = Vset_input ? δset_in.u : δset

    eqs = [
        # output connection
        connect(terminal, filter.terminal)
        # connect cc1 measurements
        [cc1.i_f_d, cc1.i_f_q] .~ _ri_to_dq(filter.i_f_r, filter.i_f_i, _δ)
        [cc1.V_C_d, cc1.V_C_q] .~ _ri_to_dq(filter.V_C_r, filter.V_C_i, _δ)
        # connect vc measurements
        [vc.V_C_d, vc.V_C_q] .~ _ri_to_dq(filter.V_C_r, filter.V_C_i, _δ)
        [vc.i_g_d, vc.i_g_q] .~ _ri_to_dq(filter.i_g_r, filter.i_g_i, _δ)
        # connect cc1 to filter inputs
        [filter.V_I_r, filter.V_I_i] .~ _dq_to_ri(cc1.V_I_d, cc1.V_I_q, _δ)
        # connect vc to cc1 reference
        cc1.i_f_ref_d ~ vc.i_f_ref_d
        cc1.i_f_ref_q ~ vc.i_f_ref_q
        # connect vc reference to setpoint/input (q-axis aligned: Vset goes into q)
        vc.V_C_ref_d ~       - (vc.i_g_d * R_virt - vc.i_g_q * X_virt)
        vc.V_C_ref_q ~ _Vset - (vc.i_g_q * R_virt + vc.i_g_d * X_virt)
    ]
    System(eqs, t; name, systems)
end

@component function CurrentSource(; name, iset_input=false)
    @parameters begin
        # LCL Filter parameters (actual electrical components in per-unit)
        Rf  = 0.01
        X_f = 0.007, [description="Inverter-side reactance (pu Ω)"]
        B_c = 0.5,  [description="Capacitive susceptance (pu S)"]
        ω0  = 2π*50
        Rg  = 0.01
        X_g = 0.007, [description="Grid-side reactance (pu Ω)"]

        # PLL parameters
        PLL_Kp = 250
        PLL_Ki = 1000

        # CC1 (inner current loop)
        CC1_KP = 0.063
        CC1_KI = 63
        CC1_F  = 0
        CC1_Fcoupl = 1

        # VC (voltage loop)
        VC_KP = 0.952
        VC_KI = 317.4
        VC_F = 0
        VC_Fcoupl = 1

        # CC2 (outer current loop)
        CC2_KP = 0.189
        CC2_KI = 0.630
        CC2_F  = 1  # Grid voltage feedforward
        CC2_Fcoupl = 1
    end

    # Current setpoint parameters (when not using inputs)
    if !iset_input
        @parameters begin
            iset_d, [guess=0, description="dq-frame current setpoint (d-axis)"]
            iset_q, [guess=0, description="dq-frame current setpoint (q-axis)"]
        end
    end

    # Create filter subsystem (always LCL for current control)
    @named filter = LCLFilter(; ω0, Rf, X_f, B_c, Rg, X_g)

    # Create PLL subsystem
    @named pll = SimplePLL(; Kp=PLL_Kp, Ki=PLL_Ki)

    # Create controller subsystems
    @named cc1 = CC1(; X_f=X_f, F=CC1_F, Fcoupl=CC1_Fcoupl, KP=CC1_KP, KI=CC1_KI)

    systems = @named begin
        terminal = Terminal()
        vc = VC(; B_c=B_c, F=VC_F, Fcoupl=VC_Fcoupl, KP=VC_KP, KI=VC_KI)
        cc2 = CC2(; X_g=X_g, F=CC2_F, Fcoupl=CC2_Fcoupl, KP=CC2_KP, KI=CC2_KI)
    end
    push!(systems, filter)
    push!(systems, pll)
    push!(systems, cc1)

    if iset_input
        @named iset_d_in = RealInput()
        @named iset_q_in = RealInput()
        append!(systems, [iset_d_in, iset_q_in])
    end

    _iset_d = iset_input ? iset_d_in.u : iset_d
    _iset_q = iset_input ? iset_q_in.u : iset_q

    eqs = [
        # Output connection
        connect(terminal, filter.terminal)

        # PLL measures terminal voltage (grid voltage)
        pll.u_r ~ terminal.u_r
        pll.u_i ~ terminal.u_i

        # CC1 measurements (using PLL angle)
        [cc1.i_f_d, cc1.i_f_q] .~ _ri_to_dq(filter.i_f_r, filter.i_f_i, pll.δ_pll)
        [cc1.V_C_d, cc1.V_C_q] .~ _ri_to_dq(filter.V_C_r, filter.V_C_i, pll.δ_pll)

        # VC measurements
        [vc.V_C_d, vc.V_C_q] .~ _ri_to_dq(filter.V_C_r, filter.V_C_i, pll.δ_pll)
        [vc.i_g_d, vc.i_g_q] .~ _ri_to_dq(filter.i_g_r, filter.i_g_i, pll.δ_pll)

        # CC2 measurements (grid current and terminal voltage)
        [cc2.i_g_d, cc2.i_g_q] .~ _ri_to_dq(filter.i_g_r, filter.i_g_i, pll.δ_pll)
        [cc2.V_g_d, cc2.V_g_q] .~ _ri_to_dq(terminal.u_r, terminal.u_i, pll.δ_pll)

        # CC1 output → filter
        [filter.V_I_r, filter.V_I_i] .~ _dq_to_ri(cc1.V_I_d, cc1.V_I_q, pll.δ_pll)

        # VC output → CC1 reference
        cc1.i_f_ref_d ~ vc.i_f_ref_d
        cc1.i_f_ref_q ~ vc.i_f_ref_q

        # CC2 output → VC reference
        vc.V_C_ref_d ~ cc2.V_C_ref_d
        vc.V_C_ref_q ~ cc2.V_C_ref_q

        # External current setpoint → CC2 reference
        cc2.i_g_ref_d ~ _iset_d
        cc2.i_g_ref_q ~ _iset_q
    ]
    System(eqs, t; name, systems)
end

# q-axis aligned convention (q-axis aligned with phase a)
# At δ=0: q = r (real/phase-a), d = -i (90° behind)
function _ri_to_dq(_r, _i, δ)
    _d =  sin(δ)*_r - cos(δ)*_i
    _q =  cos(δ)*_r + sin(δ)*_i
    return _d, _q
end
function _dq_to_ri(_d, _q, δ)
    _r =  sin(δ)*_d + cos(δ)*_q
    _i = -cos(δ)*_d + sin(δ)*_q
    return _r, _i
end


@component function DroopOuter(; name, pq_input=false, Pset=nothing, Qset=nothing)
    @parameters begin
        Vset, [description = "Voltage magnitude setpoint [pu]", guess=1]
        ω0=2π*50, [guess=2π*50, description = "Nominal frequency [pu]"]
        Kp = 0.4, [description = "Active power droop coefficient"]
        Kq = 0.04, [description = "Reactive power droop coefficient"]
        τ_p = 0.1, [description = "Active Power filter time constant [s]"]
        τ_q = 0.1, [description = "Reactive Power filter time constant [s]"]
    end
    systems = @named begin
        V_out = RealOutput()
        δ_out = RealOutput()
        P_meas_in = RealInput()
        Q_meas_in = RealInput()
    end

    if pq_input
        @named P_in = RealInput()
        @named Q_in = RealInput()
        append!(systems, (P_in, Q_in))
    else
        @parameters begin
            Pset=Pset, [guess=1, description = "Active power setpoint [pu]"]
            Qset=Qset, [guess=0, description = "Reactive power setpoint [pu]"]
        end
    end
    _Pset = pq_input ? P_in.u : Pset
    _Qset = pq_input ? Q_in.u : Qset

    @variables begin
        Pmeas(t), [description = "Measured active power [pu]", guess = 1]
        Qmeas(t), [description = "Measured reactive power [pu]", guess = 0]
        Pfilt(t), [description = "Filtered active power [pu]", guess = 1]
        Qfilt(t), [description = "Filtered reactive power [pu]", guess = 0]
        ω(t), [description = "Frequency [pu]", guess = 0]
        δ(t), [description = "Voltage angle [rad]", guess = 0]
        V(t), [description = "Voltage magnitude [pu]", guess = 1]
    end

    eqs = [
        # Power measurement from terminal quantities
        # Pmeas need to be calculated extenrally and connected as input!
        Pmeas ~ P_meas_in.u
        Qmeas ~ Q_meas_in.u

        # First-order low-pass filtering
        τ_p * Dt(Pfilt) ~ Pmeas - Pfilt
        τ_q * Dt(Qfilt) ~ Qmeas - Qfilt

        # Droop control equations
        ω ~ ω0 - Kp * (Pfilt - _Pset)
        V ~ Vset - Kq * (Qfilt - _Qset)

        # Voltage angle dynamics
        Dt(δ) ~ ω - ω0
        V_out.u ~ V
        δ_out.u ~ δ
    ]

    return System(eqs, t; name, systems)
end

function DroopInverter(; name=:droop_inv, filter_type)
    @named droop = DroopOuter(; pq_input=false)
    @named vsrc = VoltageSource(; Vset_input=true, filter_type)
    @named terminal = Terminal()

    eqs = [
        connect(vsrc.Vset_in, droop.V_out)
        connect(vsrc.δset_in, droop.δ_out)
        connect(vsrc.terminal, terminal)
        droop.P_meas_in.u ~ vsrc.filter.V_C_r * vsrc.filter.i_g_r + vsrc.filter.V_C_i * vsrc.filter.i_g_i
        droop.Q_meas_in.u ~ vsrc.filter.V_C_i * vsrc.filter.i_g_r - vsrc.filter.V_C_r * vsrc.filter.i_g_i
    ]

    System(eqs, t; name, systems=[droop, vsrc, terminal])
end

end
