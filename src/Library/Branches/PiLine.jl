"""
    PiLine

Pi-section transmission branch with series impedance, shunt admittances, and optional transformer ratios.

$(PowerDynamics.ref_source_file(@__FILE__, @__LINE__))
"""
@cached_component function PiLine(; name)
    @named src = Terminal()
    @named dst = Terminal()
    params = @parameters begin
        R=0, [description="Resistance of branch in pu"]
        X=0.1, [description="Reactance of branch in pu"]
        G_src=0, [description="Conductance of src shunt"]
        B_src=0, [description="Susceptance of src shunt"]
        G_dst=0, [description="Conductance of dst shunt"]
        B_dst=0, [description="Susceptance of dst shunt"]
        r_src=1, [description="src end transformation ratio"]
        r_dst=1, [description="dst end transformation ratio"]
        active=1, [description="Line active or at fault"]
    end
    Z = R + im*X
    Ysrc = G_src + im*B_src
    Ydst = G_dst + im*B_dst
    Vsrc = src.u_r + im*src.u_i
    Vdst = dst.u_r + im*dst.u_i
    V₁ = r_src * Vsrc
    V₂ = r_dst * Vdst
    i₁ = Ysrc * V₁
    i₂ = Ydst * V₂
    iₘ = 1/Z * (V₁ - V₂)
    isrc = (-iₘ - i₁)*r_src
    idst = ( iₘ - i₂)*r_dst
    eqs = [
        src.i_r ~ active * simplify(real(isrc))
        src.i_i ~ active * simplify(imag(isrc))
        dst.i_r ~ active * simplify(real(idst))
        dst.i_i ~ active * simplify(imag(idst))
    ]
    System(eqs, t, [], params; name, systems=[src, dst])
end
