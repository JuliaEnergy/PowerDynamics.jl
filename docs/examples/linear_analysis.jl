# # Linear Analysis Example with SimplusGT-style Synchronous Machine
#
# This model implements the SimplusGT SynchronousMachine Type 0 model:
# "constant field flux, rotor motion of torque"
#
# This implementation stays as close as possible to the MATLAB original,
# using the same parameter names, equations, and d-q reference frame convention.

using PowerDynamics
using PowerDynamics.Library
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as Dt
using OrdinaryDiffEqRosenbrock
using NetworkDynamics
using OrdinaryDiffEqNonlinearSolve

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
    end
    @components begin
        src = Terminal()
        dst = Terminal()
    end
    @variables begin
        i_mag(t)
    end
    @equations begin
        Dt(dst.i_r) ~ ωbase / L * (src.u_r - dst.u_r) - R / L * ωbase * dst.i_r + ωbase  * dst.i_i
        Dt(dst.i_i) ~ ωbase / L * (src.u_i - dst.u_i) - R / L * ωbase * dst.i_i - ωbase  * dst.i_r

        src.i_r ~ -dst.i_r
        src.i_i ~ -dst.i_i

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
        terminal.i_r ~ terminal.u_r / R
        terminal.i_i ~ terminal.u_i / R
    end
end

####
#### Testcase 1: sg connected to infinite bus through RL line
####

sm_bus = let
    @named sm = SyncMachineStatorDynamics(J=3.5, D=10, wL=0.3, R=0.003, w0=2*pi*60)
    @named shunt = ConstantRShunt(R=95.3062963614) # <- value obtained from simulink
    bus = compile_bus(MTKBus([sm, shunt]; name=:sm_bus), assume_io_coupling=false)

    @named pq = Library.PQConstraint(P=0.9, Q=0.3)
    pfmod = compile_bus(MTKBus(pq); name=:sm_bus_pfmod)

    set_pfmodel!(bus, pfmod)
    bus
end

slack_bus = pfSlack(V=0.995, name=:slack_bus)

line = let
    @named branch = DynRLLine(R=0, L=0.65, ωbase=2*pi*60)
    lm = compile_line(MTKLine(branch); name=:l12, src=:slack_bus, dst=:sm_bus)
    @named branch = PiLine(R=0, X=0.65)
    pfmod = compile_line(MTKLine(branch); name=:l12_pfmod)
    set_pfmodel!(lm, pfmod)
    lm
end


nw = Network([slack_bus, sm_bus], [line])
show_powerflow(solve_powerflow(nw))
s0 = initialize_from_pf!(nw; subverbose=true, tol=1e-7, nwtol=1e-7)
dump_initial_state(sm_bus)
jacobian_eigenvals(nw, s0) ./ (2*pi)
#=
Verdict:
We needed to include the shunt resistor to the model to make it solvable (DAE Index).
Otherwise, we cannot couple a current source with an EMT line.
This should be less of a problem in the other models because the have C-shunts!

I think this is the reason why our marginally stable mode from matlab (pure L line) is
very damped in this julia model.

   -27859.17401010486 - 60.00000030491102im   <- matlab has RE 0!
   -27859.17401010486 + 60.00000030491102im
 -0.18951881416957983 - 60.000000121045794im  <- Match Maltlab pair
 -0.18951881416957983 + 60.000000121045794im  <- Match Maltlab pair
 -0.11370497703857839 - 1.0110057816314497im  <- Match Maltlab pair
 -0.11370497703857839 + 1.0110057816314497im  <- Match Maltlab pair
=#
