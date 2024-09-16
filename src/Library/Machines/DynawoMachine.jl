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

@mtkmodel DynawoMachineParameterBase begin
    @parameters begin
        # General parameters of the synchronous machine
        UNom=24, [description="Nominal voltage in kV"]
        SNom=2220, [description="Nominal apparent power in MVA"]
        PNomTurb=2220, [description="Nominal active (turbine) power in MW"]
        PNomAlt=2200, [description="Nominal active (alternator) power in MW"]
        QNomAlt=sqrt(SNom * SNom - PNomAlt * PNomAlt), [description="Nominal reactive (alternator) power in Mvar"]
        ExcitationPu, [description="Choice of excitation base voltage"]
        H=3.5, [description="Kinetic constant = kinetic energy / rated power"]
        DPu=0, [description="Damping coefficient of the swing equation in pu"]

        # Transformer input parameters
        SnTfo=2220, [description="Nominal apparent power of the generator transformer in MVA"]
        UNomHV=24, [description="Nominal voltage on the network side of the transformer in kV"]
        UNomLV=24, [description="Nominal voltage on the generator side of the transformer in kV"]
        UBaseHV=24, [description="Base voltage on the network side of the transformer in kV"]
        UBaseLV=24, [description="Base voltage on the generator side of the transformer in kV"]
        RTfPu=0, [description="Resistance of the generator transformer in pu (base UBaseHV, SnTfo)"]
        XTfPu=0, [description="Reactance of the generator transformer in pu (base UBaseHV, SnTfo)"]

        # Mutual inductances saturation parameters, Shackshaft modelisation
        md=0.031, [description="Parameter for direct axis mutual inductance saturation modelling"]
        mq=0.031, [description="Parameter for quadrature axis mutual inductance saturation modelling"]
        nd=6.93, [description="Parameter for direct axis mutual inductance saturation modelling"]
        nq=6.93, [description="Parameter for quadrature axis mutual inductance saturation modelling"]

        # Transformer internal parameters
        RTfoPu = RTfPu * (UNomHV / UBaseHV) ^ 2 * (SNom / SnTfo), [description="Resistance of the generator transformer in pu (base SNom, UNom)"]
        XTfoPu = XTfPu * (UNomHV / UBaseHV) ^ 2 * (SNom / SnTfo), [description="Reactance of the generator transformer in pu (base SNom, UNom)"]
        # HACK: hardcoed one case
        # rTfoPu = ifelse((RTfPu > 0.0 || XTfPu > 0.0), UNomHV / UBaseHV / (UNomLV / UBaseLV), 1.0)#, [description="Ratio of the generator transformer in pu (base UBaseHV, UBaseLV)"]
        rTfoPu = ifelse((RTfPu > 0) | (XTfPu > 0), UNomHV / UBaseHV / (UNomLV / UBaseLV), 1.0)#, [description="Ratio of the generator transformer in pu (base UBaseHV, UBaseLV)"]
        # rTfoPu = UNomHV / UBaseHV / (UNomLV / UBaseLV), [description="Ratio of the generator transformer in pu (base UBaseHV, UBaseLV)"]
    end
end

@mtkmodel DynawoMachine begin
    @extend DynawoMachineParameterBase()
    @components begin
        terminal = Terminal()
        systembase = SystemBase()
        # ωRefPu = RealInput(), [input=true, description="Reference frequency in pu"]
        # PmPu = RealInput(), [input=true, description="Mechanical power in pu (base PNomTurb)"]
        # efdPu = RealInput(), [input=true, description="Input voltage of exciter winding in pu (user-selected base voltage)"]
        ωRefPu = RealInput()
        PmPu = RealInput()
        efdPu = RealInput()
    end
    @parameters begin
        # internal parameters for generator
        # Notation: Ra (resistance) + P ("'" or "Prim") + Pu (Per unit)
        Ra′Pu=0.003, [description="Armature resistance in pu"]
        Ld′Pu=0.15, [description="Direct axis stator leakage in pu"]
        Md′Pu=1.66, [description="Direct axis mutual inductance in pu"]
        LD′Pu=0.16634, [description="Direct axis damper leakage in pu"]
        RD′Pu=0.03339, [description="Direct axis damper resistance in pu"]
        Mrc′Pu=0, [description="Canay's mutual inductance in pu"]
        Lf′Pu=0.1699, [description="Excitation winding leakage in pu"]
        Rf′Pu=0.00074, [description="Excitation winding resistance in pu"]
        Lq′Pu=0.15, [description="Quadrature axis stator leakage in pu"]
        Mq′Pu=1.61, [description="Quadrature axis mutual inductance in pu"]
        LQ1′Pu=0.92815, [description="Quadrature axis 1st damper leakage in pu"]
        RQ1′Pu=0.00924, [description="Quadrature axis 1st damper resistance in pu"]
        LQ2′Pu=0.12046, [description="Quadrature axis 2nd damper leakage in pu"]
        RQ2′Pu=0.02821, [description="Quadrature axis 2nd damper resistance in pu"]
        MsalPu=0.05, [description="Constant difference between direct and quadrature axis saturated mutual inductances in pu"]
        # pu factor for excitation voltage
        Md′PuEfd, [description="Direct axis mutual inductance used to determine the excitation voltage in pu"]
        Md′PuEfdNom, [description="Direct axis mutual inductance used to determine the excitation voltage in nominal conditions in pu"]
    end
    @variables begin
        # # Input variables
        # ωRefPu(t), [input=true, description="Reference frequency in pu"]
        # PmPu(t), [input=true, description="Mechanical power in pu (base PNomTurb)"]
        # efdPu(t), [input=true, description="Input voltage of exciter winding in pu (user-selected base voltage)"]

        # Output variables
        ωPu(t), [output=true, description="Angular frequency in pu"]

        # d-q axis pu variables (base UNom, SNom)
        udPu(t), [description="Voltage of direct axis in pu"]
        uqPu(t), [description="Voltage of quadrature axis in pu"]
        idPu(t), [description="Current of direct axis in pu"]
        iqPu(t), [description="Current of quadrature axis in pu"]
        iDPu(t), [description="Current of direct axis damper in pu"]
        iQ1Pu(t), [description="Current of quadrature axis first damper sn pu"]
        iQ2Pu(t), [description="Current of quadrature axis second damper in pu"]
        ifPu(t), [description="Current of excitation winding in pu"]
        ufPu(t), [description="Voltage of exciter winding in pu (base voltage as per Kundur)"]
        λ_dPu(t), [description="Flux of direct axis in pu"]
        λ_qPu(t), [description="Flux of quadrature axis in pu"]
        λ_DPu(t), [description="Flux of direct axis damper in pu"]
        λ_fPu(t), [description="Flux of excitation winding in pu"]
        λ_Q1Pu(t), [description="Flux of quadrature axis 1st damper in pu"]
        λ_Q2Pu(t), [description="Flux of quadrature axis 2nd damper in pu"]

        # Other variables
        θ(t), [description="Rotor angle: angle between machine rotor frame and port phasor frame"]
        cmPu(t), [description="Mechanical torque in pu (base PNomTurb/ωNom)"]
        cePu(t), [description="Electrical torque in pu (base SNom/ωNom)"]
        PePu(t), [description="Electrical active power in pu (base SNom)"]

        # observables for inrospection
        PGenPu(t), [description="Active power generated in pu"]
        QGenPu(t), [description="Reactive power generated in pu"]
        PGen(t), [description="Active power generated in MW"]
        QGen(t), [description="Reactive power generated in Mvar"]
        UStatorPu(t), [description="Stator voltage magnitude in pu"]
        IStatorPu(t), [description="Stator current magnitude in pu"]

        # Saturated mutual inductances and related variables
        MdSat′Pu(t), [description="Direct axis saturated mutual inductance in pu"]
        MqSat′Pu(t), [description="Quadrature axis saturated mutual inductance in pu"]
        λ_AirGapPu(t), [description="Total air gap flux in pu"]
        λ_ADPu(t), [description="Common flux of direct axis in pu"]
        λ_AQPu(t), [description="Common flux of quadrature axis in pu"]
        mdsPu(t), [description="Direct axis saturated mutual inductance in the case when the total air gap flux is aligned on the direct axis in pu"]
        mqsPu(t), [description="Quadrature axis saturated mutual inductance in the case when the total air gap flux is aligned on the quadrature axis in pu"]
        cos2Eta(t), [description="Common flux of direct axis contribution to the total air gap flux in pu"]
        sin2Eta(t), [description="Common flux of quadrature axis contribution to the total air gap flux in pu"]
        miPu(t), [description="Intermediate axis saturated mutual inductance in pu"]
    end
    @equations begin
        # some observables
        PGenPu ~ terminal.u_r * terminal.i_r + terminal.u_i * terminal.i_i
        QGenPu ~ terminal.u_i * terminal.i_r - terminal.u_r * terminal.i_i
        PGen ~ PGenPu*systembase.SnRef;
        QGen ~ QGenPu*systembase.SnRef;
        UStatorPu ~ abs(1 / rTfoPu * (terminal.u_r + im*terminal.u_i - (terminal.i_r+im*terminal.i_i) * Complex(RTfoPu, XTfoPu) * systembase.SnRef / SNom))
        IStatorPu ~ abs(rTfoPu * (terminal.i_r+im*terminal.i_i));

        # Park's transformations
        terminal.u_r ~ sin(θ) * udPu + cos(θ) * uqPu;
        terminal.u_i ~ (-cos(θ) * udPu) + sin(θ) * uqPu;
        - terminal.i_r * systembase.SnRef / SNom ~ sin(θ) * idPu + cos(θ) * iqPu;
        - terminal.i_i * systembase.SnRef / SNom ~ (-cos(θ) * idPu) + sin(θ) * iqPu;
        # Flux linkages
        λ_dPu ~ (MdSat′Pu + Ld′Pu + XTfoPu) * idPu + MdSat′Pu * ifPu + MdSat′Pu * iDPu;
        λ_fPu ~ MdSat′Pu * idPu + (MdSat′Pu + Lf′Pu + Mrc′Pu) * ifPu + (MdSat′Pu + Mrc′Pu) * iDPu;
        λ_DPu ~ MdSat′Pu * idPu + (MdSat′Pu + Mrc′Pu) * ifPu + (MdSat′Pu + LD′Pu + Mrc′Pu) * iDPu;
        λ_qPu ~ (MqSat′Pu + Lq′Pu + XTfoPu) * iqPu + MqSat′Pu * iQ1Pu + MqSat′Pu * iQ2Pu;
        λ_Q1Pu ~ MqSat′Pu * iqPu + (MqSat′Pu + LQ1′Pu) * iQ1Pu + MqSat′Pu * iQ2Pu;
        λ_Q2Pu ~ MqSat′Pu * iqPu + MqSat′Pu * iQ1Pu + (MqSat′Pu + LQ2′Pu) * iQ2Pu;
        # Equivalent circuit equations in Park's coordinates
        udPu ~ (Ra′Pu + RTfoPu) * idPu - ωPu * λ_qPu;
        uqPu ~ (Ra′Pu + RTfoPu) * iqPu + ωPu * λ_dPu;
        ufPu ~ Rf′Pu * ifPu + Dt(λ_fPu) / systembase.ωNom;
        0 ~ RD′Pu * iDPu + Dt(λ_DPu) / systembase.ωNom;
        0 ~ RQ1′Pu * iQ1Pu + Dt(λ_Q1Pu) / systembase.ωNom;
        0 ~ RQ2′Pu * iQ2Pu + Dt(λ_Q2Pu) / systembase.ωNom;
        # Mechanical equations
        Dt(θ) ~ (ωPu - ωRefPu.u) * systembase.ωNom;
        2 * H * Dt(ωPu) ~ cmPu * PNomTurb / SNom - cePu - DPu * (ωPu - ωRefPu.u);
        cePu ~ λ_qPu * idPu - λ_dPu * iqPu;
        PePu ~ cePu * ωPu;

        # PmPu.u ~ cmPu * ωPu;
        cmPu ~ PmPu.u / ωPu;

        # Excitation voltage pu conversion
        # ufPu ~ efdPu.u * (Kuf * rTfoPu);
        # HACK: fixed to excitation type "noload"
        #=
        Kuf=if ExcitationPu == ExcitationPuType.Kundur
            1
        elseif ExcitationPu == ExcitationPuType.UserBase
            Rf′Pu / Md′PuEfd
        elseif ExcitationPu == ExcitationPuType.NoLoad
            Rf′Pu / Md′Pu
        elseif ExcitationPu == ExcitationPuType.NoLoadSaturated
            Rf′Pu * (1 + md) / Md′Pu
        else
            Rf′Pu / Md′PuEfdNom
        end;
        =#
        ufPu ~ efdPu.u * (Rf′Pu / Md′Pu * rTfoPu);
        # Mutual inductances saturation
        λ_ADPu ~ MdSat′Pu * (idPu + ifPu + iDPu);
        λ_AQPu ~ MqSat′Pu * (iqPu + iQ1Pu + iQ2Pu);
        λ_AirGapPu ~ sqrt(λ_ADPu ^ 2 + λ_AQPu ^ 2);
        mdsPu ~ Md′Pu / (1 + md * λ_AirGapPu ^ nd);
        mqsPu ~ Mq′Pu / (1 + mq * λ_AirGapPu ^ nq);
        cos2Eta ~ λ_ADPu ^ 2 / λ_AirGapPu ^ 2;
        sin2Eta ~ λ_AQPu ^ 2 / λ_AirGapPu ^ 2;
        miPu ~ mdsPu * cos2Eta + mqsPu * sin2Eta;
        MdSat′Pu ~ miPu + MsalPu * sin2Eta;
        MqSat′Pu ~ miPu - MsalPu * cos2Eta;
    end
end
