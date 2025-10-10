# PSSE IEEEST Power System Stabilizer
#
# Original work Copyright (c) 2016-2022 Luigi Vanfretti, ALSETLab, and contributors
# Original work licensed under BSD 3-Clause License
# Original source: https://github.com/OpenIPSL/OpenIPSL
#
# This Julia/PowerDynamics port maintains the same mathematical formulation
# while adapting to PowerDynamics/ModelingToolkit framework conventions.

@mtkmodel PSSE_IEEEST begin
    @structural_parameters begin
        warn = true
    end
    @parameters begin
        A1=0, [description="Power system stabilizer high frequency filter coefficient"]
        A2=0, [description="Power system stabilizer high frequency filter coefficient"]
        A3=0, [description="Power system stabilizer high frequency filter coefficient"]
        A4=0, [description="Power system stabilizer high frequency filter coefficient"]
        A5=0, [description="Power system stabilizer high frequency filter coefficient"]
        A6=0, [description="Power system stabilizer high frequency filter coefficient"]
        T1=0, [description="PSS first numerator (lead) time constant"]
        T2=0, [description="PSS first denominator (lag) time constant"]
        T3=0, [description="PSS second numerator (lead) time constant"]
        T4=0, [description="PSS second denominator (lag) time constant"]
        T5=1.65, [description="Stabilizer washout numerator time constant"]
        T6=1.65, [description="Stabilizer washout denominator time constant"]
        Ks=6.2, [description="Stabilizer Gain"]
        Lsmax=0.26, [description="Maximum output for stabilizer washout filter"]
        Lsmin=-0.1, [description="Minimum output for stabilizer washout filter"]
        Vcu=Inf, [description="Maximum power system stabilizer output"]
        Vcl=-Inf, [description="Minimum power system stabilizer output"]
        _IEEEST_active = 1, [description="Internal flag to disable/enable the stabilizer"]
    end
    begin
        _vcu = ModelingToolkit.getdefault(Vcu)
        _vcl = ModelingToolkit.getdefault(Vcl)
        if warn && _vcu != Inf && _vcl != -Inf
            @warn "You set explicit VCU/VCL limits for IEEEST. The output limiting is not active \
                   per default! You need to add a callback, which callback which watches vct(t) state \
                   in comparison to VCU/VCL and sets _IEEEST_active to 0/1 accordingly. Diable this wanrning \
                   by passing warn=false structural parmeter o the IEEEST constructor."
        end
    end
    begin
        _A1 = ModelingToolkit.getdefault(A1)
        _A2 = ModelingToolkit.getdefault(A2)
        _A3 = ModelingToolkit.getdefault(A3)
        _A4 = ModelingToolkit.getdefault(A4)
        _A5 = ModelingToolkit.getdefault(A5)
        _A6 = ModelingToolkit.getdefault(A6)

        # define coefficients (both defualts and symbolic)
        #             A6 s² + A5 s + 1
        #   ------------------------------------
        #   (A2 s² + A1 s + 1)(A4 s² + A3 s + 1)
        # or equally
        #            A6 s² + A5 s + 1
        #   --------------------------------
        #   (A2 A4) s⁴ + (A1 A4 + A2 A3) s³ + (A1 A3 + A2 + A4) s² + (A1 + A3) s + 1
        #
        # What we want to do here is to filter out all leading zeros in
        # the nom/denom vectors (otherwise the system is singular)
        # and ignore thos variables
        # I.e. the filter order is actually determined by the parameters
        # at build time!

        nom = [
            (A6) => _A6,
            (A5) => _A5,
            (1)  => 1
        ]

        denom = [
            (A2*A4)           => _A2*_A4,
            (A1*A4 + A2*A3)   => _A1*_A4 + _A2*_A3,
            (A1*A3 + A2 + A4) => _A1*_A3 + _A2 + _A4,
            (A1 + A3)         => _A1 + _A3,
            (1)               => 1
        ]

        leading_denom_nom = findfirst(x -> !isapprox(x[2], 0, atol=1e-10), nom)
        _nom = [x[2] for x in nom[leading_denom_nom:end]]
        leading_denom_denom = findfirst(x -> !isapprox(x[2], 0, atol=1e-10), denom)
        _denom = [x[2] for x in denom[leading_denom_denom:end]]

        A,B,C,D = siso_tf_to_ss(_nom, _denom)
        guesses = vcat(zeros(size(A,1)-1), 1)
    end

    @components begin
        VOTHSG_out = RealOutput()
        INPUT_in = RealInput()
        V_CT_in = RealInput() #terminal voltage

        filter12 = ss_to_mtkmodel(; A, B, C, D, guesses)
        filter2 = SimpleGain(K=1)  # Bypass second notch filter if A3=A4=0
        filter3 = LeadLag(K=1, T1=T1, T2=T2, guess=1)
        filter4 = LeadLag(K=1, T1=T3, T2=T4, guess=1)
        diff = Derivative(K=T5, T=T6)
    end
    @variables begin
        clamped_sig(t), [description="Clamped stabilizer signal"]
        vct(t), [description="Voltage measurement (disables stabilizer when outsise [Vcl,Vcu])"]
    end
    @equations begin
        filter12.in ~ INPUT_in.u
        filter3.in ~ filter12.out
        filter4.in ~ filter3.out
        diff.in ~ filter4.out * Ks
        clamped_sig  ~ clamp(diff.out, Lsmin, Lsmax)

        # output limiting is not implemented currently and should be done
        # via callback since its actually discontinuous (i.e. disables the output
        # once the operational reagion is left)
        vct ~ V_CT_in.u
        VOTHSG_out.u ~ clamped_sig * _IEEEST_active
    end
end
