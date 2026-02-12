module ComposableInverter

using ModelingToolkit: ModelingToolkit, @named, simplify, t_nounits as t, D_nounits as Dt,
                       @component
# needed for @mtkmodel
using ModelingToolkit: @mtkmodel, @variables, @parameters, @unpack, Num, System, Equation, connect, setmetadata
using ModelingToolkitStandardLibrary.Blocks: RealInput, RealOutput
using NetworkDynamics: set_mtk_defaults!

using ..PowerDynamics: Terminal

export LCFilter, LCLFilter, LFilter, VC, CC1, CC2, SimplePLL, PLL_LPF, VoltageSource, CurrentSource, SimpleGFL, SimpleGFLDC
export DroopOuter, DroopInverter

@mtkmodel LFilter begin
    @components begin
        terminal = Terminal()
    end
    @structural_parameters begin
        ω0
        Rf
        Lf  # Filter reactance [pu] (frequency-normalized inductance)
    end
    @variables begin
        i_f_r(t), [guess=0]
        i_f_i(t), [guess=0]
        i_f_mag(t)
        V_I_r(t)
        V_I_i(t)
    end
    @equations begin
        (Lf/ω0) * Dt(i_f_r) ~ V_I_r - terminal.u_r - Rf*i_f_r + Lf*i_f_i
        (Lf/ω0) * Dt(i_f_i) ~ V_I_i - terminal.u_i - Rf*i_f_i - Lf*i_f_r
        terminal.i_r ~ i_f_r
        terminal.i_i ~ i_f_i
        i_f_mag ~ sqrt(i_f_r^2 + i_f_i^2)
    end
end

@mtkmodel LCFilter begin
    @components begin
        terminal = Terminal()
    end
    @structural_parameters begin
        ω0
        Rf
        Lf  # Filter reactance [pu] (frequency-normalized inductance)
        C   # Filter susceptance [pu] (frequency-normalized capacitance)
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
        # Inductor: (Lf/ω0)*di_f/dt = V_I - V_C - Rf*i_f + Lf*W*i_f
        (Lf/ω0) * Dt(i_f_r) ~ V_I_r - V_C_r - Rf*i_f_r + Lf*i_f_i
        (Lf/ω0) * Dt(i_f_i) ~ V_I_i - V_C_i - Rf*i_f_i - Lf*i_f_r
        # Capacitor: (C/ω0)*dV_C/dt = i_f - i_g + C*W*V_C
        (C/ω0) * Dt(V_C_r) ~ i_f_r - i_g_r + C*V_C_i
        (C/ω0) * Dt(V_C_i) ~ i_f_i - i_g_i - C*V_C_r
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
        Lf  # Inverter-side reactance [pu] (frequency-normalized inductance)
        C   # Filter susceptance [pu] (frequency-normalized capacitance)
        Rg
        Lg  # Grid-side reactance [pu] (frequency-normalized inductance)
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
        (Lf/ω0) * Dt(i_f_r) ~ V_I_r - V_C_r - Rf*i_f_r + Lf*i_f_i
        (Lf/ω0) * Dt(i_f_i) ~ V_I_i - V_C_i - Rf*i_f_i - Lf*i_f_r
        # Capacitor
        (C/ω0) * Dt(V_C_r) ~ i_f_r - i_g_r + C*V_C_i
        (C/ω0) * Dt(V_C_i) ~ i_f_i - i_g_i - C*V_C_r
        # Grid-side inductor
        (Lg/ω0) * Dt(i_g_r) ~ V_C_r - connected*terminal.u_r - Rg*i_g_r + Lg*i_g_i
        (Lg/ω0) * Dt(i_g_i) ~ V_C_i - connected*terminal.u_i - Rg*i_g_i - Lg*i_g_r
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

@component function CC1(; name, Lf, F, Fcoupl=1, KP, KI)
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
        V_I .~ -Fcoupl*Lf*W*i_f + KP*(i_f_ref - i_f) + KI*γ + F*V_C
    )
    System(eqs, t; name)
end

@component function VC(; name, C, F, Fcoupl=1, KP, KI)
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
        i_f_ref .~ -Fcoupl*C*W*V_C + KP*(V_C_ref - V_C) + KI*γ + F*i_g
    )
    System(eqs, t; name)
end

@component function CC2(; name, Lg, F, Fcoupl=1, KP, KI)
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
        V_C_ref .~ -Fcoupl*Lg*W*i_g + KP*(i_g_ref - i_g) + KI*γ + F*V_g
    )
    System(eqs, t; name)
end


"""
    VoltageSource(; name, Vset_input=false, filter_type=:LC, defaults...)

Grid-forming voltage source inverter with cascaded voltage and current control.

Implements two-loop control with voltage controller (VC) commanding current references
to inner current controller (CC1). Operates in a fixed dq-frame (no PLL). Suitable for:
- Grid-forming inverters establishing voltage and frequency
- Islanded or weak-grid operation
- Droop-controlled systems (via DroopInverter wrapper)

# Parameters
- `Vset_input`: If true, voltage setpoint comes from RealInput ports. If false, uses internal Vset/δset parameters.
- `filter_type`: `:LC` for LC filter or `:LCL` for LCL filter.
- `Lf`: Inverter-side filter reactance [pu]. Related to physical inductance by Lf = ω0 * Lf_actual.
- `C`: Filter susceptance [pu]. Related to physical capacitance by C = ω0 * C_actual.
- `Lg`: Grid-side filter reactance [pu] (LCL only). Related to physical inductance by Lg = ω0 * Lg_actual.
- `ω0`: Frame angular frequency [rad/s]. Default: 2π*50 rad/s.
- Various PI controller gains (CC1_KP, CC1_KI, VC_KP, VC_KI)
- `defaults...`: Additional parameter/variable defaults (e.g., `Lf=0.01, CC1_KP=0.1`)
"""
@component function VoltageSource(; name, Vset_input=false, filter_type=:LC, defaults...)
    @parameters begin
        # Filter parameters
        Rf  = 0.01
        Lf = 0.007, [description="Inverter-side filter reactance [pu] (frequency-normalized inductance)"]
        C = 0.5,  [description="Filter susceptance [pu] (frequency-normalized capacitance)"]
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
            Lg = 0.007, [description="Grid-side filter reactance [pu] (frequency-normalized inductance)"]
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
        @named filter = LCFilter(; ω0, Rf, Lf, C)
    else  # :LCL
        @named filter = LCLFilter(; ω0, Rf, Lf, C, Rg, Lg)
    end

    @named cc1 = CC1(; Lf=Lf, F=CC1_F, Fcoupl=CC1_Fcoupl, KP=CC1_KP, KI=CC1_KI)
    systems = @named begin
        terminal = Terminal()
        vc = VC(; C=C, F=VC_F, Fcoupl=VC_Fcoupl, KP=VC_KP, KI=VC_KI)
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
    sys = System(eqs, t; name, systems)
    set_mtk_defaults!(sys, defaults)
    return sys
end

"""
    CurrentSource(; name, iset_input=false, defaults...)

Grid-following current source inverter with triple-loop cascaded control and PLL.

Implements three-loop control: outer current controller (CC2) → voltage controller (VC)
→ inner current controller (CC1). Uses PLL for grid synchronization. Suitable for:
- Grid-following inverters injecting controlled current
- Renewable energy sources (solar, wind) in grid-connected mode
- Active/reactive power control applications

# Parameters
- `iset_input`: If true, current setpoint comes from RealInput ports. If false, uses internal iset_d/iset_q parameters.
- `Lf`: Inverter-side filter reactance [pu]. Related to physical inductance by Lf = ω0 * Lf_actual.
- `C`: Filter susceptance [pu]. Related to physical capacitance by C = ω0 * C_actual.
- `Lg`: Grid-side filter reactance [pu]. Related to physical inductance by Lg = ω0 * Lg_actual.
- `ω0`: Frame angular frequency [rad/s]. Default: 2π*50 rad/s.
- PLL and controller gains (PLL_Kp, PLL_Ki, CC1_*, VC_*, CC2_*)
- `defaults...`: Additional parameter/variable defaults (e.g., `Lf=0.01, PLL_Kp=100`)
"""
@component function CurrentSource(; name, iset_input=false, defaults...)
    @parameters begin
        # LCL Filter parameters
        Rf  = 0.01
        Lf = 0.007, [description="Inverter-side filter reactance [pu] (frequency-normalized inductance)"]
        C = 0.5,  [description="Filter susceptance [pu] (frequency-normalized capacitance)"]
        ω0  = 2π*50
        Rg  = 0.01
        Lg = 0.007, [description="Grid-side filter reactance [pu] (frequency-normalized inductance)"]

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
    @named filter = LCLFilter(; ω0, Rf, Lf, C, Rg, Lg)

    # Create PLL subsystem
    @named pll = SimplePLL(; Kp=PLL_Kp, Ki=PLL_Ki)

    # Create controller subsystems
    @named cc1 = CC1(; Lf=Lf, F=CC1_F, Fcoupl=CC1_Fcoupl, KP=CC1_KP, KI=CC1_KI)

    systems = @named begin
        terminal = Terminal()
        vc = VC(; C=C, F=VC_F, Fcoupl=VC_Fcoupl, KP=VC_KP, KI=VC_KI)
        cc2 = CC2(; Lg=Lg, F=CC2_F, Fcoupl=CC2_Fcoupl, KP=CC2_KP, KI=CC2_KI)
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
    sys = System(eqs, t; name, systems)
    set_mtk_defaults!(sys, defaults)
    return sys
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

"""
    DroopInverter(; name=:droop_inv, filter_type)

Grid-forming inverter with droop control for frequency and voltage regulation.

Wraps VoltageSource with DroopOuter controller implementing P-f and Q-V droop characteristics.
Suitable for:
- Grid-forming inverters in microgrids
- Virtual synchronous machine (VSM) implementations
- Parallel inverter operation with power sharing

# Parameters
- `filter_type`: `:LC` or `:LCL` filter topology (passed to VoltageSource)
- Droop controller parameters: `Kp` (P-f droop), `Kq` (Q-V droop), `τ_p`, `τ_q` (power filter time constants)
- Filter parameters inherited from VoltageSource: `Lf`, `C`, `Lg`, `ω0`
"""
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

"""
    SimplePLL(; name, Kp, Ki)

Basic PLL that tracks the grid voltage angle by driving the q-axis voltage to zero.
Uses a PI controller on the angle error signal `u_q`.

# Parameters
- `Kp`: Proportional gain
- `Ki`: Integral gain
"""
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

"""
    PLL_LPF(; name, Kp, Ki, τ_lpf)

PLL with an additional low-pass filter on the frequency output.
Uses a PI controller on the angle error, followed by a first-order LPF on the
estimated frequency `ω`. Matches the SimplusGT MATLAB PLL structure.

# Parameters
- `Kp`: Proportional gain
- `Ki`: Integral gain
- `τ_lpf`: Low-pass filter time constant [s]
"""
@component function PLL_LPF(; name, Kp, Ki, τ_lpf)
    vars = @variables begin
        ω_pll_i(t), [guess=0]
        ω(t), [guess=2π*50]
        θ(t), [guess=0]
        e_ang(t)
        u_r(t)
        u_i(t)
    end

    eqs = [
        # Error signal: v_q in MATLAB d-axis aligned convention = -v_d in our q-axis convention
        e_ang ~ -sin(θ)*u_r + cos(θ)*u_i
        Dt(ω_pll_i) ~ e_ang * Ki
        Dt(ω) ~ (ω_pll_i + e_ang*Kp - ω) / τ_lpf
        Dt(θ) ~ ω
    ]
    System(eqs, t; name)
end

"""
    SimpleGFL(; name, iset_input=false, defaults...)

Simplified grid-following inverter with L-filter, current control, and PLL.

Implements single-loop current control (CC1) with PLL synchronization and simple L-filter.
Suitable for:
- Basic grid-following inverter models
- Simplified renewable energy source representation
- Studies where detailed filter dynamics are not critical

# Parameters
- `iset_input`: If true, current setpoint comes from RealInput ports. If false, uses internal iset_d/iset_q parameters.
- `Lf`: Filter reactance [pu]. Related to physical inductance by Lf = ω0 * Lf_actual.
- `ω0`: Frame angular frequency [rad/s]. Default: 2π*50 rad/s.
- PLL and CC1 controller gains
- `defaults...`: Additional parameter/variable defaults (e.g., `Lf=0.05, PLL_Kp=50`)
"""
@component function SimpleGFL(; name, iset_input=false, defaults...)
    @parameters begin
        Rf  = 0.01
        Lf = 0.05, [description="Filter reactance [pu] (frequency-normalized inductance)"]
        ω0  = 2π*50
        # PLL
        PLL_Kp = 2π*10
        PLL_Ki = (2π*10)^2/4
        PLL_τ_lpf = 1/(2π*300)
        # CC1
        CC1_KP = 0.6
        CC1_KI = 565.5
        CC1_F  = 0
        CC1_Fcoupl = 0
    end

    if !iset_input
        @parameters begin
            iset_d, [guess=0, description="dq-frame current setpoint (d-axis)"]
            iset_q, [guess=0, description="dq-frame current setpoint (q-axis)"]
        end
    end

    @named filter = LFilter(; ω0, Rf, Lf)
    @named pll = PLL_LPF(; Kp=PLL_Kp, Ki=PLL_Ki, τ_lpf=PLL_τ_lpf)
    @named cc1 = CC1(; Lf=Lf, F=CC1_F, Fcoupl=CC1_Fcoupl, KP=CC1_KP, KI=CC1_KI)

    systems = @named begin
        terminal = Terminal()
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
        # Terminal connection
        connect(terminal, filter.terminal)
        # PLL measures terminal voltage
        pll.u_r ~ terminal.u_r
        pll.u_i ~ terminal.u_i
        # CC1 measurements in dq frame (using PLL angle)
        [cc1.i_f_d, cc1.i_f_q] .~ _ri_to_dq(filter.i_f_r, filter.i_f_i, pll.θ)
        [cc1.V_C_d, cc1.V_C_q] .~ _ri_to_dq(terminal.u_r, terminal.u_i, pll.θ)
        # CC1 output → filter voltage input
        [filter.V_I_r, filter.V_I_i] .~ _dq_to_ri(cc1.V_I_d, cc1.V_I_q, pll.θ)
        # Current setpoint
        cc1.i_f_ref_d ~ _iset_d
        cc1.i_f_ref_q ~ _iset_q
    ]
    sys = System(eqs, t; name, systems)
    set_mtk_defaults!(sys, defaults)
    return sys
end

"""
    SimpleGFLDC(; name, defaults...)

Grid-following inverter with DC-link dynamics and active power control.

Extends SimpleGFL with DC-link capacitor model and PI controller for DC voltage regulation.
The DC voltage controller generates active current reference (q-axis). Suitable for:
- Inverters with significant DC-link capacitance
- Renewable sources with DC power input (solar PV, battery storage)
- Studies requiring DC-side transient analysis

# Parameters
- `Lf`: Filter reactance [pu]. Related to physical inductance by Lf = ω0 * Lf_actual.
- `ω0`: Frame angular frequency [rad/s]. Default: 2π*50 rad/s.
- `C_dc`: DC-link capacitance [pu]
- `V_dc`: DC voltage reference [pu]
- `kp_v_dc`, `ki_v_dc`: DC voltage PI controller gains
- `P_dc`: External DC power draw [pu] (solved at initialization)
- PLL and CC1 controller gains
- `defaults...`: Additional parameter/variable defaults (e.g., `Lf=0.03, V_dc=2.0`)
"""
@component function SimpleGFLDC(; name, defaults...)
    @parameters begin
        Rf  = 0.01
        Lf = 0.03, [description="Filter reactance [pu] (frequency-normalized inductance)"]
        ω0  = 2π*50
        # PLL
        PLL_Kp = 2π*5
        PLL_Ki = (2π*5)^2/4
        PLL_τ_lpf = 1/(2π*300)
        # CC1
        CC1_KP = 0.36
        CC1_KI = 135.0
        CC1_F  = 0
        CC1_Fcoupl = 0
        # DC link
        C_dc = 1.25
        V_dc = 2.5, [description="DC voltage reference"]
        kp_v_dc, [description="DC voltage PI proportional gain"]
        ki_v_dc, [description="DC voltage PI integral gain"]
        # Reactive current setpoint (constant, from PF)
        iset_d, [guess=0, description="dq-frame reactive current setpoint"]
        # External DC power (solved at initialization)
        P_dc, [guess=0, description="External DC power draw"]
    end

    @variables begin
        v_dc_state(t), [guess=2.5, description="DC capacitor voltage"]
        v_dc_i(t), [guess=0, description="DC voltage PI integrator"]
    end

    @named filter = LFilter(; ω0, Rf, Lf)
    @named pll = PLL_LPF(; Kp=PLL_Kp, Ki=PLL_Ki, τ_lpf=PLL_τ_lpf)
    @named cc1 = CC1(; Lf=Lf, F=CC1_F, Fcoupl=CC1_Fcoupl, KP=CC1_KP, KI=CC1_KI)

    systems = @named begin
        terminal = Terminal()
    end
    push!(systems, filter)
    push!(systems, pll)
    push!(systems, cc1)

    # DC PI output → active current reference (q-axis in Julia convention)
    iset_q_dc = (V_dc - v_dc_state)*kp_v_dc + v_dc_i

    # P_ac: power from converter to AC filter (V_I · i_f in dq frame)
    P_ac = cc1.V_I_d * cc1.i_f_d + cc1.V_I_q * cc1.i_f_q

    eqs = [
        # Terminal connection
        connect(terminal, filter.terminal)
        # PLL measures terminal voltage
        pll.u_r ~ terminal.u_r
        pll.u_i ~ terminal.u_i
        # CC1 measurements in dq frame (using PLL angle)
        [cc1.i_f_d, cc1.i_f_q] .~ _ri_to_dq(filter.i_f_r, filter.i_f_i, pll.θ)
        [cc1.V_C_d, cc1.V_C_q] .~ _ri_to_dq(terminal.u_r, terminal.u_i, pll.θ)
        # CC1 output → filter voltage input
        [filter.V_I_r, filter.V_I_i] .~ _dq_to_ri(cc1.V_I_d, cc1.V_I_q, pll.θ)
        # Current setpoint: d from parameter (reactive), q from DC PI (active)
        cc1.i_f_ref_d ~ iset_d
        cc1.i_f_ref_q ~ iset_q_dc
        # DC link dynamics
        C_dc * Dt(v_dc_state) ~ (P_ac - P_dc) / v_dc_state
        Dt(v_dc_i) ~ (V_dc - v_dc_state) * ki_v_dc
    ]
    sys = System(eqs, t; name, systems)
    set_mtk_defaults!(sys, defaults)
    return sys
end

end
