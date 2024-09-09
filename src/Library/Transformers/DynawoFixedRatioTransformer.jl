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

@mtkmodel DynawoTransformerParameters begin
    @parameters begin
        RPu=0, [description="Resistance of the generator transformer in pu (base U2Nom, SnRef)"]
        XPu=0.00675, [description="Reactance of the generator transformer in pu (base U2Nom, SnRef)"]
        GPu=0, [description="Conductance of the generator transformer in pu (base U2Nom, SnRef)"]
        BPu=0, [description="Susceptance of the generator transformer in pu (base U2Nom, SnRef)"]
    end
end

@mtkmodel DynawoFixedRatioTransformer begin
    @extend DynawoTransformerParameters()
    @components begin
        src = Terminal()
        dst = Terminal()
    end
    @parameters begin
        rTfoPu=1, [description="Transformation ratio in pu: U2/U1 in no load conditions"]
    end
    @equations begin
        (rTfoPu^2)*src.u_r ~ RPu*src.i_r - XPu*src.i_i + rTfoPu*dst.u_r
        (rTfoPu^2)*src.u_i ~ RPu*src.i_i + XPu*src.i_r + rTfoPu*dst.u_i
        src.i_r ~ (-dst.i_r - BPu*dst.u_i + GPu*dst.u_r)*rTfoPu
        src.i_i ~ (-dst.i_i + BPu*dst.u_r + GPu*dst.u_i)*rTfoPu
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
    rTfoPu, [description="Transformation ratio in pu: U2/U1 in no load conditions"]
end

simplify(rTfoPu * rTfoPu * t1u ~ rTfoPu * t2u + (RPu + im*XPu) * t1i)
simplify(t1i ~ rTfoPu * ((GPu+im*BPu) * t2u - t2i))
=#
