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
        terminal1 = Terminal()
        terminal2 = Terminal()
    end
    @parameters begin
        RPu=0, [description="Resistance in pu (base SnRef)"]
        XPu=0.022522, [description="Reactance in pu (base SnRef)"]
        GPu=0, [description="Half-conductance in pu (base SnRef)"]
        BPu=0, [description="Half-susceptance in pu (base SnRef)"]
    end
    @equations begin
        # (terminal2.i_r + BPu*terminal2.u_i - GPu*terminal2.u_r)*RPu - (terminal2.i_i - BPu*terminal2.u_r - GPu*terminal2.u_i)*XPu ~ terminal2.u_r - terminal1.u_r
        # (terminal2.i_r + BPu*terminal2.u_i - GPu*terminal2.u_r)*XPu + (terminal2.i_i - BPu*terminal2.u_r - GPu*terminal2.u_i)*RPu ~ terminal2.u_i - terminal1.u_i
        # (terminal1.i_r + BPu*terminal1.u_i - GPu*terminal1.u_r)*RPu - (terminal1.i_i - BPu*terminal1.u_r - GPu*terminal1.u_i)*XPu ~ -terminal2.u_r + terminal1.u_r
        # (terminal1.i_r + BPu*terminal1.u_i - GPu*terminal1.u_r)*XPu + (terminal1.i_i - BPu*terminal1.u_r - GPu*terminal1.u_i)*RPu ~ -terminal2.u_i + terminal1.u_i

        # (terminal1.i_r + im*terminal1.i_i) ~ ((terminal1.u_r + im*terminal1.u_i) - (terminal2.u_r + im*terminal2.u_i))/(RPu+im*XPu) + (GPu+im*BPu) * (terminal1.u_r + im*terminal1.u_i)
        # (terminal2.i_r + im*terminal2.i_i) ~ ((terminal2.u_r + im*terminal2.u_i) - (terminal1.u_r + im*terminal1.u_i))/(RPu+im*XPu) + (GPu+im*BPu) * (terminal2.u_r + im*terminal2.u_i)

        simplify(-terminal1.i_r ~ real(((terminal1.u_r + im*terminal1.u_i) - (terminal2.u_r + im*terminal2.u_i))/(RPu+im*XPu) + (GPu+im*BPu) * (terminal1.u_r + im*terminal1.u_i)))
        simplify(-terminal1.i_i ~ imag(((terminal1.u_r + im*terminal1.u_i) - (terminal2.u_r + im*terminal2.u_i))/(RPu+im*XPu) + (GPu+im*BPu) * (terminal1.u_r + im*terminal1.u_i)))
        simplify(-terminal2.i_r ~ real(((terminal2.u_r + im*terminal2.u_i) - (terminal1.u_r + im*terminal1.u_i))/(RPu+im*XPu) + (GPu+im*BPu) * (terminal2.u_r + im*terminal2.u_i)))
        simplify(-terminal2.i_i ~ imag(((terminal2.u_r + im*terminal2.u_i) - (terminal1.u_r + im*terminal1.u_i))/(RPu+im*XPu) + (GPu+im*BPu) * (terminal2.u_r + im*terminal2.u_i)))
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
