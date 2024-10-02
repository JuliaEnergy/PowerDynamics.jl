@mtkmodel PiLine begin
    @parameters begin
        R, [description="Resistance of branch in pu (base unclear?)"]
        X, [description="Reactance of branch in pu (base unclear?)"]
        G_src, [description="Conductance of src shunt (base unclear?)"]
        B_src, [description="Susceptance of src shunt (base unclear?)"]
        G_dst, [description="Conductance of dst shunt (base unclear?)"]
        B_dst, [description="Susceptance of dst shunt (base unclear?)"]
        r_src=1, [description="Src end transformation ratio"]
        r_dst=1, [description="Src end transformation ratio"]
        active=1, [description="Line active or at fault"]
    end
    @components begin
        src = Terminal()
        dst = Terminal()
    end
    begin
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
        isrc = -iₘ - i₁
        idst =  iₘ - i₂
    end
    @equations begin
        src.i_r ~ active * simplify(real(isrc))
        src.i_i ~ active * simplify(imag(isrc))
        dst.i_r ~ active * simplify(real(idst))
        dst.i_i ~ active * simplify(imag(idst))
    end
end
