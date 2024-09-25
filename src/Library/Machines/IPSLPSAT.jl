@mtkmodel IPSLPSATOrder4 begin
    @components begin
        terminal = Terminal()
        # inputs
        pm = RealInput() # mechanical power [pu]
        vf = RealInput(guess=1) # field voltage input [pu]
        # outputs
        δout = RealOutput() # rotor angle
        ωout = RealOutput() # rotor speed [pu]
        v_mag_out = RealOutput() # terminal voltage [pu]
        Pout = RealOutput() # active power [pu]
        Qout = RealOutput() # reactive power [pu]
    end
    @parameters begin
        # base P
        ra, [description="Armature resistance"]
        x1d, [description="d-axis transient reactance"]
        M, [description="Mechanical starting time, 2H [Ws/VA]"]
        D, [description="Damping coefficient"]
        S_b, [description="System power basis in MVA"]
        V_b, [description="System voltage basis in kV"]
        ω_b, [description="System base frequency in rad/s"]
        Sn, [description="Machine power rating in MVA"]
        Vn, [description="Machine voltage rating in kV"]
        # 4th order params
        xd=1.9, [description="d-axis synchronous reactance"]
        xq=1.7, [description="q-axis synchronous reactance"]
        x1q=0.5, [description="q-axis transient reactance"]
        T1d0=8, [description="d-axis open circuit transient time constant"]
        T1q0=0.8, [description="q-axis open circuit transient time constant"]
    end
    @variables begin
        # base vars
        v_arg(t), [description="Generator terminal angle"]
        vd(t), [description="d-axis voltage"]
        vq(t), [description="q-axis voltage"]
        id(t), [guess=0, description="d-axis current"]
        iq(t), [guess=0, description="q-axis current"]
        pe(t), [description="Electrical power transmitted through air gap"]

        δ(t), [guess=0, description="rotor angle"]
        ω(t), [guess=1, description="rotor speed [pu]"]
        v_mag(t), [description="terminal voltage [pu]"]
        P(t), [description="active power [pu]"]
        Q(t), [description="reactive power [pu]"]

        # 4th order vars
        e1d(t), [guess=0, description="d-axis transient voltage"]
        e1q(t), [guess=1, description="q-axis transient voltage"]
    end
    @equations begin
        # Park's transformations
        terminal.u_r ~ ( sin(δ) * vd + cos(δ) * vq) * V_b/Vn
        terminal.u_i ~ (-cos(δ) * vd + sin(δ) * vq) * V_b/Vn
        -terminal.i_r ~ (-sin(δ) * id - cos(δ) * iq) * Ibase(S_b, V_b)/Ibase(Sn, Vn)
        -terminal.i_i ~ ( cos(δ) * id - sin(δ) * iq) * Ibase(S_b, V_b)/Ibase(Sn, Vn)

        # observables
        v_mag ~ sqrt(terminal.u_r^2 + terminal.u_i^2)
        v_arg ~ atan(terminal.u_i, terminal.u_r)
        P ~ terminal.u_r*terminal.i_r + terminal.u_i*terminal.i_i
        Q ~ terminal.u_i*terminal.i_r - terminal.u_r*terminal.i_i

        #outputs
        Pout.u ~ P
        Qout.u ~ Q
        v_mag_out.u ~ v_mag
        δout.u ~ δ
        ωout.u ~ ω

        # swing equation
        pe ~ (vq + ra*iq)*iq + (vd + ra*id)*id
        Dt(δ) ~ ω_b*(ω - 1)
        Dt(ω) ~ (pm.u/S_b*Sn - pe - D*(ω - 1))/M

        # internal transients
        Dt(e1q) ~ ((-e1q) - (xd - x1d)*id + vf.u*V_b/Vn)/T1d0;
        Dt(e1d) ~ ifelse(abs(xd - x1q) < 1e-16,
                         ((-e1d) + (xq - x1q)*iq)/T1q0,
                         (-e1d)/T1q0)
        e1q ~ vq + ra*iq + x1d*id
        e1d ~ vd + ra*id - x1q*iq
    end
end
