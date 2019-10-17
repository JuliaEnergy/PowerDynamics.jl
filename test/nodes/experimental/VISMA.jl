using Test: @testset, @test
using SymPy: @syms
using PowerDynamics: VISMA, construct_vertex, symbolsof
using LinearAlgebra: I

include("../NodeTestBase.jl")

@testset "VISMA Tests" begin
    @syms Ω R_e L_e R_d R_q L_d L_q R_D R_Q L_D L_Q M_Dd M_Qq M_ed M_eD Z_p H positive=true
    @syms P_m u_e real=true
    @syms θ dθ ω dω  real=true
    @syms Ψ_d dΨ_d Ψ_q dΨ_q Ψ_D dΨ_D Ψ_Q dΨ_Q Ψ_e dΨ_e complex=true
    @syms H_INVALID R_D_invalid negative=true

    @test_throws AssertionError construct_vertex(VISMA(P_m=P_m,Ω=Ω,u_e=u_e,R_e=R_e,L_e=L_e,R_d=R_d,R_q=R_q,L_d=L_d,L_q=L_q,R_D=R_D_invalid,R_Q=R_Q,L_D=L_D,L_Q=L_Q,M_Dd=M_Dd,M_Qq=M_Qq,M_ed=M_ed,M_eD=M_eD,Z_p=Z_p,H=H))
    @test_throws AssertionError construct_vertex(VISMA(P_m=P_m,Ω=Ω,u_e=u_e,R_e=R_e,L_e=L_e,R_d=R_d,R_q=R_q,L_d=L_d,L_q=L_q,R_D=R_D,R_Q=R_Q,L_D=L_D,L_Q=L_Q,M_Dd=M_Dd,M_Qq=M_Qq,M_ed=M_ed,M_eD=M_eD,Z_p=Z_p,H=H_INVALID))

    visma = VISMA(P_m=P_m,Ω=Ω,u_e=u_e,R_e=R_e,L_e=L_e,R_d=R_d,R_q=R_q,L_d=L_d,L_q=L_q,R_D=R_D,R_Q=R_Q,L_D=L_D,L_Q=L_D,M_Dd=M_Dd,M_Qq=M_Qq,M_ed=M_ed,M_eD=M_eD,Z_p=Z_p,H=H)
    visma_vertex = construct_vertex(visma)

    @test symbolsof(visma) == [:u_r, :u_i,:θ,:Ψ_d,:Ψ_q,:Ψ_D,:Ψ_Q,:Ψ_e,:ω]
    @test visma_vertex.mass_matrix == [0,0,1,1,1,1,1,1,1]

    smoketest_rhs(visma_vertex, int_x=[θ,Ψ_d,Ψ_q,Ψ_D,Ψ_Q,Ψ_e,ω], int_dx=[dθ,dΨ_d,dΨ_q,dΨ_D,dΨ_Q,dΨ_e,dω])
end
