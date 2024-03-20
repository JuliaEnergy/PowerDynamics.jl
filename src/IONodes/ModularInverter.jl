using Unitful
using BlockSystems
using .MIComponents
using PowerDynamics.IOComponents

export MIParameters
export Inverter, InverterCS
export VSwithLoad, DroopControl, Synchronverter, FixedVoltage, FixedCurrent, PLLCurrent
export PerfectSource, PerfectCSource, PT1Source, VoltageControlFOM, CurrentControlFOM, VoltageControlCSFOM, ConstantPower

MIBase = let ω0base = 2π*50u"rad/s",
             Sbase = 10000u"W",
             Vbase = sqrt(3)*230u"V",
             Ibase = Sbase/(Vbase) |> u"A",
             Cbase = Ibase/Vbase,
             Lbase = Vbase/Ibase,
             Rbase = (Vbase^2)/Sbase
    (;ω0base, Sbase, Vbase, Ibase, Cbase, Lbase, Rbase)
end

MIParameters = Dict(
    # electrical parameters
    :Rf => uconvert(NoUnits, 0.01u"Ω"/MIBase.Rbase),
    :Rg => uconvert(NoUnits, 0.01u"Ω"/MIBase.Rbase),
    :Lf => ustrip(u"s", 350e-6u"H"/MIBase.Lbase),
    :Lg => ustrip(u"s", 350e-6u"H"/MIBase.Lbase),
    :C  => ustrip(u"s", 100e-6u"F"/MIBase.Cbase),
    :ω0 => ustrip(u"rad/s", MIBase.ω0base),
    # CC1 control parameters
    :CC1₊KP => 1    * ustrip(u"A/V", MIBase.Ibase/MIBase.Vbase),
    :CC1₊KI => 1000 * ustrip(u"A/V", MIBase.Ibase/MIBase.Vbase),
    # PR configs
    :CC1₊PR1₊K  => 50*ustrip(u"A/V", MIBase.Ibase/MIBase.Vbase),
    :CC1₊PR1₊ω  => 3*2π*50,
    :CC1₊PR1₊ωc => 5,
    :CC1₊PR2₊K  => 25*ustrip(u"A/V", MIBase.Ibase/MIBase.Vbase),
    :CC1₊PR2₊ω  => 6*2π*50,
    :CC1₊PR2₊ωc => 5,
    :CC1₊PR3₊K  => 10*ustrip(u"A/V", MIBase.Ibase/MIBase.Vbase),
    :CC1₊PR3₊ω  => 9*2π*50,
    :CC1₊PR3₊ωc => 5,
    :CC1₊PR4₊K  => 5*ustrip(u"A/V", MIBase.Ibase/MIBase.Vbase),
    :CC1₊PR4₊ω  => 12*2π*50,
    :CC1₊PR4₊ωc => 5,
    # VC control parameters
    :VC₊KP => 0.06 * ustrip(u"V/A", MIBase.Vbase/MIBase.Ibase),
    :VC₊KI => 20  * ustrip(u"V/A", MIBase.Vbase/MIBase.Ibase),
    # CC2 control parameters
    :CC2₊KP => 3  * ustrip(u"A/V", MIBase.Ibase/MIBase.Vbase),
    :CC2₊KI => 10 * ustrip(u"A/V", MIBase.Ibase/MIBase.Vbase)
)

function Inverter(inner, outer; name=Symbol(string(outer.name)*"_"*string(inner.name)))
    @assert BlockSpec([:i_i, :i_r, :u_ref_r, :u_ref_i], [:u_r, :u_i])(inner) "Inner ctrl loop :$(inner.name)  does not meet expectation."
    @assert BlockSpec([], [:u_ref_r, :u_ref_i]; in_strict=false)(outer)  "Outer ctrl loop :$(outer.name) does not meet expectation."

    outerinputs = ModelingToolkit.getname.(outer.inputs)
    :u_r ∈ outerinputs && @error "Outer shoud use :u_meas_r instead of :u_r"
    :u_i ∈ outerinputs && @error "Outer shoud use :u_meas_i instead of :u_i"
    :i_r ∈ outerinputs && @error "Outer shoud use :i_meas_r instead of :i_r"
    :i_i ∈ outerinputs && @error "Outer shoud use :i_meas_i instead of :i_i"

    meas = MIComponents.UIMeas()

    sys = IOSystem(:autocon, [outer, inner, meas];
                   outputs=[inner.u_r, inner.u_i],
                   globalp=[:i_r, :i_i],
                   name)
    closed = connect_system(sys)

    @assert BlockSpec([:i_i, :i_r], [:u_r, :u_i])(closed) "Closed loop does not match expectation! $closed"

    return closed
end

function InverterCS(inner, outer; name=Symbol(string(outer.name)*"_"*string(inner.name)))
    # @assert BlockSpec([:u_i, :u_r, :i_ref_r, :i_ref_i], [:i_r, :i_i])(inner) "Inner ctrl loop :$(inner.name)  does not meet expectation."
    # @assert BlockSpec([], [:i_ref_r, :i_ref_i]; in_strict=false)(outer)  "Outer ctrl loop :$(outer.name) does not meet expectation."

    outerinputs = ModelingToolkit.getname.(outer.inputs)
    :u_r ∈ outerinputs && @error "Outer shoud use :u_meas_r instead of :u_r"
    :u_i ∈ outerinputs && @error "Outer shoud use :u_meas_i instead of :u_i"
    :i_r ∈ outerinputs && @error "Outer shoud use :i_meas_r instead of :i_r"
    :i_i ∈ outerinputs && @error "Outer shoud use :i_meas_i instead of :i_i"

    meas = MIComponents.UIMeas()
    sys = IOSystem(:autocon, [outer, inner, meas];
                   outputs=[inner.i_r, inner.i_i],
                   globalp=[:u_r, :u_i],
                   name)
    closed = connect_system(sys)

    @assert BlockSpec([:u_i, :u_r], [:i_r, :i_i])(closed) "Closed loop does not match expectation! $closed"

    return closed
end

####
#### Outer Control Loops
####
"""
    DroopControl(; params...)

Return block for droop control outer ocntrol.
"""
function DroopControl(; params...)
    @named Pfil = IOComponents.LowPassFilter(; τ=:τ_P, input=:P_meas, output=:P_fil)
    @named Pdroop = IOComponents.DroopControl(; x_ref=:P_ref, K=:K_P, u_ref=:ω_ref, u=:ω, x=:P_fil)
    Psys = @connect Pfil.P_fil => Pdroop.P_fil outputs=:remaining name=:Psys

    @named Qfil = IOComponents.LowPassFilter(; τ=:τ_Q, input=:Q_meas, output=:Q_fil)
    @named Qdroop = IOComponents.DroopControl(; x_ref=:Q_ref, K=:K_Q, u_ref=:V_ref, u=:V, x=:Q_fil)
    Qsys = @connect Qfil.Q_fil => Qdroop.Q_fil outputs=:remaining name=:Qsys

    @named PQmeas = IOComponents.Power(; u_r=:u_meas_r, u_i=:u_meas_i,
                                     i_r=:i_meas_r, i_i=:i_meas_i,
                                     P=:P_meas, Q=:Q_meas)

    refgen = MIComponents.VRefGen()

    @named droop = IOSystem(:autocon, [PQmeas, Psys, Qsys, refgen], outputs=:remaining)
    con = connect_system(droop)

    return replace_vars(con, params)
end

function Synchronverter(; params...)
    # Frequency Loop
    @variables t ΔT(t) ω(t) θ(t)
    @parameters Te(t) P_ref Dp ω0 J
    dt = Differential(t)
    @named floop = IOBlock([ΔT ~ P_ref/ω0 - Dp*ω - Te,
                            dt(ω) ~ 1/J * ΔT,
                            dt(θ) ~ ω],
                           [Te],
                           [ω, θ])

    # Voltage Loop
    @variables t ΔQ(t) MfIf(t) V(t)
    @parameters u_meas_r(t) u_meas_i(t) Q(t) Q_ref Dq V_ref Kv
    @named vloop = IOBlock([V ~ sqrt(u_meas_r^2 + u_meas_i^2),
                            ΔQ ~ Q_ref - Q + Dq*(V_ref - V),
                            dt(MfIf) ~ ΔQ/Kv],
                           [Q, u_meas_r, u_meas_i],
                           [MfIf])
    BlockSystems.WARN[] = true

    # main model
    @variables t Te(t) u_ref_r(t) u_ref_i(t) Q(t)
    @parameters i_meas_r(t) i_meas_i(t) MfIf(t) ω(t) θ(t) ω0
    @named machine = IOBlock([Te      ~  sqrt(3/2) * MfIf * (cos(θ)*i_meas_r + sin(θ)*i_meas_i),
                              u_ref_r ~  sqrt(3/2) * MfIf * (ω0 + ω) * cos(θ),
                              u_ref_i ~  sqrt(3/2) * MfIf * (ω0 + ω) * sin(θ),
                              Q       ~ -sqrt(3/2) * MfIf * (ω0 + ω) * (-sin(θ)*i_meas_r + cos(θ)*i_meas_i)],
                             [ω, θ, i_meas_r, i_meas_i, MfIf],
                             [u_ref_r, u_ref_i, Q, Te])

    @named syncvert = IOSystem(:autocon, [floop, vloop, machine], outputs=:remaining, globalp=[:ω0])
    con = connect_system(syncvert)
    return replace_vars(con, params)
end

function FixedVoltage(; params...)
    @variables t u_ref_r(t) u_ref_i(t)
    @parameters u_fix_r u_fix_i
    blk = IOBlock([u_ref_r ~ u_fix_r,
             u_ref_i ~ u_fix_i],
            [], [u_ref_r, u_ref_i],
            name=:FixedU)
    replace_vars(blk, params)
end

function PLLCurrent(; params...)
    pll = MIComponents.JuanPLL(; Kp=250, Ki=1000, u_r=:u_meas_r, u_i=:u_meas_i)

    pol2cart = IOComponents.Polar2Cart(;arg=:δ_pll, mag=:i_mag_ref, x=:i_ref_r, y=:i_ref_i)

    con = IOSystem([pll.δ_pll => pol2cart.δ_pll],
                   [pll, pol2cart];
                   outputs = [pol2cart.i_ref_r, pol2cart.i_ref_i],
                   name=:PLL_current) |> connect_system
    con = make_iparam(con, :i_mag_ref)
end

function ConstantPower(; params...)
    @variables t i_ref_r(t) i_ref_i(t)
    @parameters u_meas_r(t) u_meas_i(t) P_ref Q_ref
    blk = IOBlock([i_ref_r ~ ( P_ref*u_meas_r + Q_ref*u_meas_i)/(u_meas_r^2+u_meas_i^2),
                   i_ref_i ~ ( P_ref*u_meas_i - Q_ref*u_meas_r)/(u_meas_r^2+u_meas_i^2)],
        [u_meas_r, u_meas_i], [i_ref_r, i_ref_i],
        name = :ConstantPower)
end

export FiltConstantPower
function FiltConstantPower(; params...)
    lpf_r = IOComponents.LowPassFilter(input=:u_meas_r, output=:u_filt_r)
    lpf_i = IOComponents.LowPassFilter(input=:u_meas_i, output=:u_filt_i)

    @variables t i_ref_r(t) i_ref_i(t)
    @parameters u_filt_r(t) u_filt_i(t) P_ref Q_ref
    blk = IOBlock([i_ref_r ~ ( P_ref*u_filt_r + Q_ref*u_filt_i)/(u_filt_r^2+u_filt_i^2),
                   i_ref_i ~ ( P_ref*u_filt_i - Q_ref*u_filt_r)/(u_filt_r^2+u_filt_i^2)],
        [u_filt_r, u_filt_i], [i_ref_r, i_ref_i],
        name = :PowerCalc)
    sys = IOSystem(:autocon, [lpf_r, lpf_i, blk]; globalp=[:τ], outputs=:remaining)
    con = connect_system(sys)
    if !isempty(params)
        con = replace_vars(con, params)
    end
    con
end

####
#### Inner Control Loops
####
"""
    PT1Source(; params...)

Create Schiffer Voltage source which follows angle directly but
amplitude with lag.
"""
function PT1Source(; params...)
    @variables t A(t) u_r(t) u_i(t)
    @parameters τ u_ref_r(t) u_ref_i(t) i_r(t) i_i(t)
    dt = Differential(t)
    blk = IOBlock([dt(A) ~ 1/τ*(√(u_ref_r^2 + u_ref_i^2) - A),
                   u_r ~ A/√(u_ref_r^2 + u_ref_i^2) * u_ref_r,
                   u_i ~ A/√(u_ref_r^2 + u_ref_i^2) * u_ref_i],
                  [u_ref_r, u_ref_i, i_r, i_i],
                  [u_r, u_i],
                  name=:PT1Src,
                  warn=false)

    if !isempty(params)
        blk = replace_vars(blk, params)
    end
    return blk
end

"""
    PerfectSource()

Perfect Voltage source which follows the reference directly.
"""
function PerfectSource(; params...)
    @variables t u_r(t) u_i(t)
    @parameters u_ref_r(t) u_ref_i(t) i_r(t) i_i(t)
    Vsource = IOBlock([u_r ~ u_ref_r,
                       u_i ~ u_ref_i],
                      [u_ref_r, u_ref_i, i_r, i_i],
                      [u_r, u_i],
                      name=:VSrc,
                      warn=false)
end

"""
    PerfectCSource()

Perfect current source which follows the reference directly.
"""
function PerfectCSource(; params...)
    @variables t i_r(t) i_i(t)
    @parameters i_ref_r(t) i_ref_i(t) u_r(t) u_i(t)
    Vsource = IOBlock([i_r ~ i_ref_r,
                       i_i ~ i_ref_i],
                      [i_ref_r, i_ref_i, u_r, u_i],
                      [i_r, i_i],
                      name=:CSrc,
                      warn=false)
end


"""
    VoltageControlFOM(;params..)

"""
function VoltageControlFOM(;pr=true, ff=true)
    LC  = MIComponents.LC()
    CC1 = if pr
        MIComponents.CC1_PR()
    else
        MIComponents.CC1()
    end
    VC  = MIComponents.VC()
    vcinner = IOSystem(:autocon,
                       [LC,CC1,VC];
                       globalp=[:ω0, :Rf, :Rg, :Lf, :Lg, :C, :i_g_r, :i_g_i],
                       name=:VCFOM,
                       outputs=[:LC₊V_C_r, :LC₊V_C_i],
                       autopromote=false) |> connect_system

    vcinner = replace_vars(vcinner;
                           :LC₊V_C_r => :u_r,
                           :LC₊V_C_i => :u_i,
                           :VC₊V_C_ref_r => :u_ref_r,
                           :VC₊V_C_ref_i => :u_ref_i,
                           :i_g_r => :i_r,
                           :i_g_i => :i_i)
    params = copy(MIParameters)
    params[:CC1₊F] = false
    params[:VC₊F] = ff
    vcinner = replace_vars(vcinner, params; warn=false)
end

"""
    VoltageControlCSFOM(;params..)

"""
function VoltageControlCSFOM(;pr=true, ff=true, params...)
    LCL = MIComponents.LCL()
    CC1 = if pr
        MIComponents.CC1_PR()
    else
        MIComponents.CC1()
    end
    VC  = MIComponents.VC()
    vcinner = IOSystem(:autocon,
                       [LCL,CC1,VC];
                       globalp=[:ω0, :Rf, :Rg, :Lf, :Lg, :C, :V_g_r, :V_g_i],
                       name=:VCFOM,
                       outputs=[:LCL₊i_g_r, :LCL₊i_g_i],
                       autopromote=false) |> connect_system

    vcinner = replace_vars(vcinner;
                           :V_g_r => :u_r,
                           :V_g_i => :u_i,
                           :VC₊V_C_ref_r => :u_ref_r,
                           :VC₊V_C_ref_i => :u_ref_i,
                           :LCL₊i_g_r => :i_r,
                           :LCL₊i_g_i => :i_i)
    params = copy(MIParameters)
    params[:CC1₊F] = false
    params[:VC₊F] = ff
    vcinner = replace_vars(vcinner, params; warn=false)
end

"""
    CurrentControlFOM(;params..)

"""
function CurrentControlFOM(; pr=true)
    LCL = MIComponents.LCL()
    CC1 = if pr
        MIComponents.CC1_PR()
    else
        MIComponents.CC1()
    end
    VC  = MIComponents.VC()
    CC2 = MIComponents.CC2()
    ccinner = IOSystem(:autocon,
                       [LCL,CC1,VC,CC2];
                       globalp=[:ω0, :Rf, :Rg, :Lf, :Lg, :C, :V_g_r, :V_g_i],
                       name=:CCFOM,
                       outputs=[:LCL₊i_g_r, :LCL₊i_g_i],
                       autopromote=false) |> connect_system

    ccinner = replace_vars(ccinner;
                           :V_g_r => :u_r,
                           :V_g_i => :u_i,
                           :CC2₊i_g_ref_r => :i_ref_r,
                           :CC2₊i_g_ref_i => :i_ref_i,
                           :LCL₊i_g_r => :i_r,
                           :LCL₊i_g_i => :i_i)

    ccparams = copy(MIParameters)
    ccparams[:CC1₊F] = false
    ccparams[:VC₊F] = false
    ccparams[:CC2₊F] = true

    replace_vars(ccinner, ccparams; warn=false)
end
