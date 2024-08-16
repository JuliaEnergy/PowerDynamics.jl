export DQBus, DQSwing, DQLine, DQPiLine

@mtkmodel DQBus begin
    @variables begin
        u_r(t), [description = "d-voltage", output=true]
        u_i(t), [description = "q-voltage", output=true]
        i_r(t), [description = "d-current", input=true]
        i_i(t), [description = "d-current", input=true]
    end
end

rotm(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]
@mtkmodel DQSwing begin
    @extend DQBus()
    @variables begin
        ω(t) = 0.0, [description = "Rotor frequency"]
        θ(t) = 0.0, [description = "Rotor angle"]
        Pel(t), [description = "Electrical Power"]
    end
    @parameters begin
        M = 1, [description = "Inertia"]
        D = 0.1, [description = "Damping"]
        Pmech, [description = "Mechanical Power"]
        V = 1.0, [description = "Voltage magnitude"]
    end
    @equations begin
        Dt(θ) ~ ω
        Dt(ω) ~ 1/M * (Pmech - D*ω + Pel)
        Pel ~ real(Complex(u_r, u_i) * conj(Complex(i_r, i_i)))
        [u_r, u_i] ~ rotm(θ) * [V; 0]
    end
end

@mtkmodel DQLine begin
    @variables begin
        src_u_r(t), [description = "src d-voltage", input=false]
        src_u_i(t), [description = "src q-voltage", input=false]
        dst_u_r(t), [description = "dst d-voltage", input=false]
        dst_u_i(t), [description = "dst q-voltage", input=false]
        src_i_r(t), [description = "src d-current", output=true]
        src_i_i(t), [description = "src d-current", output=true]
        dst_i_r(t), [description = "dst d-current", output=true]
        dst_i_i(t), [description = "dst d-current", output=true]
    end
end

@mtkmodel DQPiLine begin
    @extend DQLine()
    @parameters begin
        R = 1.0, [description = "Resistance"]
        X = 1.0, [description = "Reactance"]
        src_B = 0.0, [description = "Shunt susceptance at src end"]
        dst_B = 0.0, [description = "Shunt susceptance at dst end"]
    end
    @equations begin
        src_i_r ~ -real(-(src_u_r + im*src_u_i - dst_u_r - im*dst_u_i)/(R + im*X) - (im*src_B)*(src_u_r + im*src_u_i))
        src_i_i ~ -imag(-(src_u_r + im*src_u_i - dst_u_r - im*dst_u_i)/(R + im*X) - (im*src_B)*(src_u_r + im*src_u_i))
        dst_i_r ~ -real(+(src_u_r + im*src_u_i - dst_u_r - im*dst_u_i)/(R + im*X) - (im*dst_B)*(dst_u_r + im*dst_u_i))
        dst_i_i ~ -imag(+(src_u_r + im*src_u_i - dst_u_r - im*dst_u_i)/(R + im*X) - (im*dst_B)*(dst_u_r + im*dst_u_i))
    end
end
