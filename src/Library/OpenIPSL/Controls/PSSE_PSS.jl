@mtkmodel NotchFilter_Canonical begin
    @structural_parameters begin
        A  # numerator s term
        B  # numerator s^2 term
        C  # denominator s term
        D  # denominator s^2 term
    end
    @variables begin
        in(t), [input=true]
        out(t), [output=true]
        x1(t), [guess=1]  # x1 = y
        x2(t), [guess=0]  # x2 = y'
    end
    @equations begin
        Dt(x1) ~ x2
        Dt(x2) ~ (in + A*x2 + B*((in - x1 - C*x2)/D) - C*x2 - x1)/D
        out ~ x1
    end
end

@mtkmodel PSSE_IEEEST begin
    @structural_parameters begin
        A1     # Notch filter parameters
        A2     # Notch filter parameters
        A3     # Notch filter parameters
        A4     # Notch filter parameters
        A5     # Notch filter parameters
        A6     # Notch filter parameters
        T1     # Lead/lag time constant, sec
        T2     # Lead/lag time constant, sec
        T3     # Lead/lag time constant, sec
        T4     # Lead/lag time constant, sec
        T5     # Wahsout numerator time constant, sec
        T6     # Washout denominator time constant, sec
        Ks     # Stabilizer gains
        Lsmax  # Maximum stabilizer output, pu
        Lsmin  # Minimum stabilizer output, pu
        Vcu    # Stabilizer input cutoff threshold, pu
        Vcl    # Stabilizer input cutoff threshold, pu
    end
    begin
        @warn "PSSE_IEEEST is not fully implemented! Filter bypass logic is broken"
    end
    @components begin
        VOTHSG_out = RealOutput()
        INPUT_in = RealInput()
        V_CT_in = RealInput() #terminal voltage

        # filter1 = SimpleGain(K=1)
        filter1 = NotchFilter_Canonical(A=A5, B=A6, C=A1, D=A2)
        # filter2 = NotchFilter_Canonical(A=A5, B=A6, C=A3, D=A4)
        filter2 = SimpleGain(K=1)  # Bypass second notch filter if A3=A4=0
        filter3 = LeadLag(K=1, T1=T1, T2=T2, guess=1)
        filter4 = LeadLag(K=1, T1=T3, T2=T4, guess=1)
        diff = Derivative(K=T5, T=T6)
    end
    @variables begin
        clamped_sig(t), [description="Clamped stabilizer signal"]
        active(t), [description="Stabilizer active flag"]
        vct(t), [description="Voltage measurement for stabilizer input"]
    end
    @equations begin
        filter1.in ~ INPUT_in.u
        filter2.in ~ filter1.out
        filter3.in ~ filter2.out
        filter4.in ~ filter3.out
        diff.in ~ filter4.out * Ks
        clamped_sig  ~ clamp(diff.out, Lsmin, Lsmax)
        vct ~ V_CT_in.u
        VOTHSG_out.u ~ clamped_sig*0
        # VOTHSG_out.u ~ clamped_sig
        # active ~ ifelse(Vcl < vct, 1, 0)*ifelse(vct < Vcu, 1, 0)
        # VOTHSG_out.u ~ clamped_sig * active
    end
end


function RationalTransferFunction(num::AbstractVector, den::AbstractVector; name=nothing)
    m = length(num)-1
    n = length(den)-1
    if m > n
        error("Numerator degree > denominator degree not allowed (non-proper).")
    end

    # Feedthrough and remainder
    r, D = split_proper(num, den)
    if length(r) < n
        r = vcat(zeros(n - length(r)), r)
    end
    C = reverse(r)'

    # Leading denominator coefficient
    leading = den[1]
    @assert !iszero(leading)

    # Companion canonical form
    A = convert(Matrix{Any}, zeros(n,n))
    if n > 1
        A[1:end-1, 2:end] .= LinearAlgebra.I(n-1)
    end
    # scale last row by 1/leading
    A[end, :] .= -reverse(den[2:end]) / leading
    B = convert(Matrix{Any}, zeros(n,1))
    B[end] = 1 / leading

    # Symbolic system
    @variables in(t) out(t)
    _xs_names = [ Symbol("x", NetworkDynamics.subscript(i)) for i in 1:n ]
    x = map(_xs_names) do nm
        Symbolics.variable(nm; T=Symbolics.FnType)(t)
    end
    ∂x = Dt.(x)
    eqs = vcat(
        ∂x .~ A*x .+ B*[in],
        [out] .~ (length(C)>0 ? C*x : 0) .+ D*[in]
    )
    eqs = Symbolics.simplify.(eqs)
    allp = reduce(∪, Symbolics.get_variables.(vcat(num, den)))

    return System(eqs, t, vcat(x, [in, out]), allp; name=name)
end

# Symbolic-safe split of numerator into feedthrough D and strictly proper remainder r(s)
function split_proper(num::AbstractVector, den::AbstractVector)
    n = length(den) - 1
    m = length(num) - 1
    if m < n
        # strictly proper: feedthrough D = 0
        D = zero(eltype(num))
        r = num
        # pad with zeros on left to match denominator order later
    elseif m == n
        # proper: feedthrough D = leading coefficient ratio
        D = num[1] / den[1]
        # remainder r(s) = num(s) - D * den(s)
        r_full = [ num[i] - D * den[i] for i in 1:length(num) ]
        # drop leading coefficient (coefficient of s^n)
        r = r_full[2:end]
    else
        error("Numerator degree > denominator degree not allowed")
    end
    return r, D
end
