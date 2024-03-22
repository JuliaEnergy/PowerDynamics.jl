module MIComponents

using BlockSystems
using LinearAlgebra

"""
    VRefGen(; name=:Vrefgen, renamings...)

Create dq-reference from (ω, V) reference.
"""
function VRefGen(; name=:Vrefgen, renamings...)
    @variables t δ(t) u_ref_r(t) u_ref_i(t)
    @parameters V(t) ω(t)
    dt = Differential(t)
    block = IOBlock([dt(δ) ~ ω,
                     u_ref_r ~ V*cos(δ),
                     u_ref_i ~ V*sin(δ)],
                    [V, ω],
                    [u_ref_r, u_ref_i];
                    name)
    replace_vars(block; renamings...)
end

"""
    UIMeas(; name=:uimeas, renamings...)

Creates a very simple block which is mainly there for renaming.

TODO: is this necessary? Maybe make i_i and i_r "globalp" inputs
"""
function UIMeas(; name=:ui_meas, renamings...)
    @variables t u_meas_r(t) u_meas_i(t) i_meas_r(t) i_meas_i(t)
    @parameters u_r(t) u_i(t) i_i(t) i_r(t)
    block = IOBlock([u_meas_r ~ u_r,
                     u_meas_i ~ u_i,
                     i_meas_r ~ i_r,
                     i_meas_i ~ i_i],
                    [u_r, u_i, i_r, i_i],
                    [u_meas_r, u_meas_i, i_meas_r, i_meas_i];
                    name)

    replace_vars(block; renamings...)
end

function JuanPLL(; name=:pll, renamings...)
    @variables t δ_pll(t) ω_pll(t) u_q(t)
    @parameters u_r(t) u_i(t) Kp Ki
    dt = Differential(t)
    block = IOBlock([u_q ~ 1/sqrt(u_r^2+u_i^2)*(u_r*sin(-δ_pll) + u_i*cos(-δ_pll)),
                     dt(δ_pll) ~ ω_pll + Kp*u_q,
                     dt(ω_pll) ~ Ki*u_q
                     ],
                    [u_r, u_i],
                    [δ_pll, ω_pll];
                    name)
    block = substitute_algebraic_states(block)
    replace_vars(block; renamings...)
end

function L(; kwargs...)
    @parameters ω0 Rf Rg Lf Lg
    @variables t i_f_r(t) i_f_i(t)
    @parameters V_I_r(t) V_I_i(t) V_C_r(t) V_C_i(t)
    dt = Differential(t)

    W = [0 1; -1 0]
    i_f = [i_f_r, i_f_i]

    x = [i_f_r, i_f_i]
    u = [V_I_r, V_I_i, V_C_r, V_C_i]
    A = W*ω0 - Rf/Lf * I
    B = [1/Lf*I(2) -1/Lf*I]

    L = IOBlock(dt.(x) .~ A*x + B*u,
                 [V_I_r, V_I_i, V_C_r, V_C_i],
                 [i_f_r, i_f_i],
                 name = :L)
    replace_vars(L; kwargs...)
end

function LC(; name=:LC, kwargs...)
    @parameters ω0 Rf Rg Lf Lg C
    @variables t i_f_r(t) i_f_i(t) V_C_r(t) V_C_i(t)
    @parameters  V_I_r(t) V_I_i(t) i_g_r(t) i_g_i(t)
    dt = Differential(t)

    W = [0 1; -1 0]
    i_f = [i_f_r, i_f_i]
    V_C = [V_C_r, V_C_i]

    x = [i_f_r, i_f_i, V_C_r, V_C_i]
    u = [V_I_r, V_I_i, i_g_r, i_g_i]
    A = [W*ω0 - Rf/Lf*I   -1/Lf*I ;
         1/C*I              W*ω0]
    B = [1/Lf*I     zeros(2,2);
         zeros(2,2)  -1/C*I]

    LC = IOBlock(dt.(x) .~ A*x + B*u,
                 [i_g_r, i_g_i, V_I_r, V_I_i],
                 [i_f_r, i_f_i, V_C_r, V_C_i];
                 name)

    replace_vars(LC; kwargs...)
end

function LCL(; kwargs...)
    @parameters ω0 Rf Rg Lf Lg C
    @variables t i_f_r(t) i_f_i(t) i_g_r(t) i_g_i(t) V_C_r(t) V_C_i(t)
    @parameters  V_I_r(t) V_I_i(t) V_g_r(t) V_g_i(t)
    dt = Differential(t)
    W = [0 1; -1 0]

    i_f = [i_f_r, i_f_i]
    V_C = [V_C_r, V_C_i]
    i_g = [i_g_r, i_g_i]

    x = [i_f_r, i_f_i, V_C_r, V_C_i, i_g_r, i_g_i]
    u = [V_I_r, V_I_i, V_g_r, V_g_i]
    A = [W*ω0 - Rf/Lf*I   -1/Lf*I   zeros(2,2)   ;
         1/C*I              W*ω0    -1/C*I      ;
         zeros(2,2)       1/Lg*I  W*ω0 - Rg/Lg*I ]
    B = [1/Lf*I     zeros(2,2);
         zeros(2,2) zeros(2,2);
         zeros(2,2)  -1/Lg*I  ]

    LCL = IOBlock(dt.(x) .~ A*x + B*u,
                  [V_g_r, V_g_i, V_I_r, V_I_i],
                  [i_f_r, i_f_i, V_C_r, V_C_i, i_g_r, i_g_i],
                  name = :LCL)
    replace_vars(LCL; kwargs...)
end

function CC1(; kwargs...)
    @parameters ω0 Rf Rg Lf Lg C
    @variables t γ_r(t) γ_i(t) V_I_r(t) V_I_i(t)
    @parameters KP KI i_f_r(t) i_f_i(t) i_f_ref_r(t) i_f_ref_i(t) V_C_r(t) V_C_i(t) F
    dt = Differential(t)
    W = [0 1; -1 0]

    γ = [γ_r, γ_i]
    i_f = [i_f_r, i_f_i]
    i_f_ref = [i_f_ref_r, i_f_ref_i]
    V_I = [V_I_r, V_I_i]
    V_C = [V_C_r, V_C_i]

    CC1 = IOBlock(vcat(dt.(γ) .~ i_f_ref - i_f,
                      V_I .~ -Lf*ω0*W*i_f + KP*(i_f_ref - i_f) + KI*γ + F*V_C),
                 [i_f_r, i_f_i, i_f_ref_r, i_f_ref_i, V_C_r, V_C_i],
                 [V_I_r, V_I_i],
                 name=:CC1)

    replace_vars(CC1; kwargs...)
end

function VC(; kwargs...)
    @parameters ω0 Rf Rg Lf Lg C
    @variables t γ_r(t) γ_i(t) i_f_ref_r(t)  i_f_ref_i(t)
    @parameters KP KI F V_C_r(t) V_C_i(t) V_C_ref_r(t) V_C_ref_i(t) i_g_r(t) i_g_i(t) ω0 C
    dt = Differential(t)
    W = [0 1; -1 0]

    γ = [γ_r, γ_i]
    V_C = [V_C_r, V_C_i]
    i_g = [i_g_r, i_g_i]
    V_C_ref = [V_C_ref_r, V_C_ref_i]
    i_f_ref = [i_f_ref_r, i_f_ref_i]

    VC = IOBlock(vcat(dt.(γ) .~ V_C_ref - V_C,
                      i_f_ref .~ - C*ω0*W*V_C + KP*(V_C_ref - V_C) + KI*γ + F*i_g),
                 [i_g_r, i_g_i, V_C_r, V_C_i, V_C_ref_r, V_C_ref_i],
                 [i_f_ref_r, i_f_ref_i],
                 name=:VC)
    replace_vars(VC; kwargs...)
end

function CC2(;kwargs...)
    @parameters ω0 Rf Rg Lf Lg C
    @variables t γ_r(t) γ_i(t) V_C_ref_r(t) V_C_ref_i(t)
    @parameters KP KI i_g_r(t) i_g_i(t) i_g_ref_r(t) i_g_ref_i(t) V_g_r(t) V_g_i(t) F
    dt = Differential(t)
    W = [0 1; -1 0]

    γ = [γ_r, γ_i]
    i_g = [i_g_r, i_g_i]
    i_g_ref = [i_g_ref_r, i_g_ref_i]
    V_C_ref = [V_C_ref_r, V_C_ref_i]
    V_g = [V_g_r, V_g_i]

    CC2 = IOBlock(vcat(dt.(γ) .~ i_g_ref - i_g,
                      V_C_ref .~ -Lf*ω0*W*i_g + KP*(i_g_ref - i_g) + KI*γ + F*V_g),
                 [i_g_r, i_g_i, i_g_ref_r, i_g_ref_i, V_g_r, V_g_i],
                 [V_C_ref_r, V_C_ref_i],
                 name=:CC2)
    replace_vars(CC2; kwargs...)
end

function CC1_PR()
    dq = CC1()

    renamings = (; u_r=:i_f_r, u_ref_r=:i_f_ref_r, u_i=:i_f_i, u_ref_i=:i_f_ref_i)
    pr1 = PR(;name=:PR1, renamings...)
    pr2 = PR(;name=:PR2, renamings...)
    pr3 = PR(;name=:PR3, renamings...)
    pr4 = PR(;name=:PR4, renamings...)

    @variables t
    @parameters dq_V_I_r(t) dq_V_I_i(t) pr_V_I_r(t) pr_V_I_i(t)
    @variables V_I_r(t) V_I_i(t)

    add = IOBlock([V_I_r ~ dq_V_I_r + pr_V_I_r,
                   V_I_i ~ dq_V_I_i + pr_V_I_i],
                  [dq_V_I_r, dq_V_I_i, pr_V_I_r, pr_V_I_i],
                  [V_I_r, V_I_i],
                  name=:DQPR_add)

    CC1pr = IOSystem([dq.V_I_r => add.dq_V_I_r,
                      dq.V_I_i => add.dq_V_I_i,
                      pr1.y_r + pr2.y_r + pr3.y_r + pr4.y_r => add.pr_V_I_r,
                      pr1.y_i + pr2.y_i + pr3.y_i + pr4.y_i => add.pr_V_I_i],
                     [add, dq, pr1, pr2, pr3, pr4],
                     globalp=[:i_f_r, :i_f_i, :i_f_ref_r, :i_f_ref_i],
                     outputs=[add.V_I_r, add.V_I_i],
                     namespace_map=[add.V_I_r=>:V_I_r,
                                    add.V_I_i=>:V_I_i], name=:CC1_pr)
    connect_system(CC1pr; name=:CC1)
end

function PR(;name=:PR, kwargs...)
    @variables t x₁_r(t) x₂_r(t) x₁_i(t) x₂_i(t) y_r(t) y_i(t)
    @parameters K ω ωc u_r(t) u_i(t) u_ref_r(t) u_ref_i(t)
    dt = Differential(t)
    blk = IOBlock([dt(x₁_r) ~ x₂_r,
                   dt(x₂_r) ~ -ω^2*x₁_r - 2*ωc*x₂_r + u_ref_r - u_r,
                   y_r ~ 2*K*ωc*x₂_r,
                   dt(x₁_i) ~ x₂_i,
                   dt(x₂_i) ~ -ω^2*x₁_i - 2*ωc*x₂_i + u_ref_r - u_i,
                   y_i ~ 2*K*ωc*x₂_i],
                  [u_r, u_i, u_ref_r, u_ref_i], [y_r, y_i]; name)

    replace_vars(blk; kwargs...)
end

end # module
