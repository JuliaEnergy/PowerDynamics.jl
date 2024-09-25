#=
* Copyright (c) 2015-2019, RTE (http://www.rte-france.com)
* See AUTHORS.txt
* All rights reserved.
* This Source Code Form is subject to the terms of the Mozilla Public
* License, v. 2.0. If a copy of the MPL was not distributed with this
* file, you can obtain one at http://mozilla.org/MPL/2.0/.
* SPDX-License-Identifier: MPL-2.0
*
* This file is based on Dynawo, an hybrid C++/Modelica open source time domain simulation tool for power systems.
* It was ported to the Julia Programming language by the authors of PowerDynamics.jl.
=#

@mtkmodel DynawoPiLine begin
    @components begin
        src = Terminal()
        dst = Terminal()
    end
    @parameters begin
        RPu=0, [description="Resistance in pu (base SnRef)"]
        XPu=0.022522, [description="Reactance in pu (base SnRef)"]
        GPu=0, [description="Half-conductance in pu (base SnRef)"]
        BPu=0, [description="Half-susceptance in pu (base SnRef)"]
    end
    @equations begin
        # (dst.i_r + BPu*dst.u_i - GPu*dst.u_r)*RPu - (dst.i_i - BPu*dst.u_r - GPu*dst.u_i)*XPu ~ dst.u_r - src.u_r
        # (dst.i_r + BPu*dst.u_i - GPu*dst.u_r)*XPu + (dst.i_i - BPu*dst.u_r - GPu*dst.u_i)*RPu ~ dst.u_i - src.u_i
        # (src.i_r + BPu*src.u_i - GPu*src.u_r)*RPu - (src.i_i - BPu*src.u_r - GPu*src.u_i)*XPu ~ -dst.u_r + src.u_r
        # (src.i_r + BPu*src.u_i - GPu*src.u_r)*XPu + (src.i_i - BPu*src.u_r - GPu*src.u_i)*RPu ~ -dst.u_i + src.u_i

        # (src.i_r + im*src.i_i) ~ ((src.u_r + im*src.u_i) - (dst.u_r + im*dst.u_i))/(RPu+im*XPu) + (GPu+im*BPu) * (src.u_r + im*src.u_i)
        # (dst.i_r + im*dst.i_i) ~ ((dst.u_r + im*dst.u_i) - (src.u_r + im*src.u_i))/(RPu+im*XPu) + (GPu+im*BPu) * (dst.u_r + im*dst.u_i)

        simplify(-src.i_r ~ real(((src.u_r + im*src.u_i) - (dst.u_r + im*dst.u_i))/(RPu+im*XPu) + (GPu+im*BPu) * (src.u_r + im*src.u_i)))
        simplify(-src.i_i ~ imag(((src.u_r + im*src.u_i) - (dst.u_r + im*dst.u_i))/(RPu+im*XPu) + (GPu+im*BPu) * (src.u_r + im*src.u_i)))
        simplify(-dst.i_r ~ real(((dst.u_r + im*dst.u_i) - (src.u_r + im*src.u_i))/(RPu+im*XPu) + (GPu+im*BPu) * (dst.u_r + im*dst.u_i)))
        simplify(-dst.i_i ~ imag(((dst.u_r + im*dst.u_i) - (src.u_r + im*src.u_i))/(RPu+im*XPu) + (GPu+im*BPu) * (dst.u_r + im*dst.u_i)))
    end
end

#=
using ModelingToolkit
using ModelingToolkit: t_nounits as t
@variables begin
    t1u(t)::Complex
    t2u(t)::Complex
    t1i(t)::Complex
    t2i(t)::Complex
end
@parameters begin
    RPu, [description="Resistance of the generator transformer in pu (base U2Nom, SnRef)"]
    XPu, [description="Reactance of the generator transformer in pu (base U2Nom, SnRef)"]
    GPu, [description="Conductance of the generator transformer in pu (base U2Nom, SnRef)"]
    BPu, [description="Susceptance of the generator transformer in pu (base U2Nom, SnRef)"]
end

simplify((RPu+im*XPu) * (t2i - (GPu+im*BPu) * t2u) ~ t2u - t1u)
simplify((RPu+im*XPu) * (t1i - (GPu+im*BPu) * t1u) ~ t1u - t2u)

simplify(rTfoPu * rTfoPu * t1u ~ rTfoPu * t2u + (RPu + im*XPu) * t1i)
simplify(t1i ~ rTfoPu * ((GPu+im*BPu) * t2u - t2i))
=#
