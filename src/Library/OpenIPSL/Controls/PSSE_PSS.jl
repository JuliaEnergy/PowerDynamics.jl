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

ModelingToolkit.@component ss_to_mtkmodel(; A, B, C, D, kwargs...) = ss_to_mtkmodel(A, B, C, D; kwargs...)
function ss_to_mtkmodel(A, B, C, D; name=nothing, guesses=zeros(size(A,1)))
    t = ModelingToolkit.t_nounits
    Dt = ModelingToolkit.D_nounits

    n = size(A, 1)
    @assert size(D) == (1, 1) "Only SISO systems supported"

    # Symbolic system
    @variables in(t) out(t)
    _xs_names = [ Symbol("x", NetworkDynamics.subscript(i)) for i in 1:n]
    # dont use Symbolics.variables as it does not create all necessary metadata?
    # also needs to set guess in @variables not with MTK.setguess
    x = map(zip(_xs_names, guesses)) do (_name, _guess)
        only(@variables $(_name)(t) [guess=_guess])
    end

    ∂x = Dt.(x)
    eqs = vcat(
        ∂x .~ A*x .+ B*[in],
        [out] .~ (length(C)>0 ? C*x : 0) .+ D*[in]
    )
    eqs = Symbolics.simplify.(eqs)
    allp = mapreduce(Symbolics.get_variables, ∪, Iterators.flatten((A,B,C,D)))

    return System(eqs, t, vcat(x, [in, out]), allp; name=name)
end

#=
Taken and adapted from SymbolicControlSystems.jl

Copyright (c) 2020 Fredrik Bagge Carlson, MIT License
=#
function siso_tf_to_ss(num0, den0)
    T = Base.promote_type(eltype(num0), eltype(den0))

    # truncate leading zeros
    num0 = num0[findfirst(!iszero, num0):end]
    den0 = den0[findfirst(!iszero, den0):end]

    # check if it is proper
    denorder = length(den0) - 1
    numorder = length(num0) - 1
    if numorder > denorder
        error("Numerator degree > denominator degree not allowed (non-proper).")
    end


    # Normalize the numerator and denominator to allow realization of transfer functions
    # that are proper, but not strictly proper
    num = num0 ./ den0[1]
    den = den0 ./ den0[1]

    N = length(den) - 1 # The order of the rational function f

    # Get numerator coefficient of the same order as the denominator
    bN = length(num) == N+1 ? num[1] : zero(eltype(num))

    @views if N == 0 #|| num == zero(Polynomial{T})
        A = zeros(T, 0, 0)
        B = zeros(T, 0, 1)
        C = zeros(T, 1, 0)
    else
        A = LinearAlgebra.diagm(1 => ones(T, N-1))
        A[end, :] .= .-reverse(den)[1:end-1]

        B = zeros(T, N, 1)
        B[end] = one(T)

        C = zeros(T, 1, N)
        C[1:min(N, length(num))] = reverse(num)[1:min(N, length(num))]
        C[:] .-= bN .* reverse(den)[1:end-1] # Can index into polynomials at greater inddices than their length
    end
    D = fill(bN, 1, 1)

    return A, B, C, D
end
