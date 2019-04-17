using PowerDynBase
using NetworkDynamics
using Test

@testset "construct_node_dynamics should return ODEVertex" begin
        sw_eq = SwingEq(H=10.0, P=2.0, D=0.5, Ω=0.1)
        swing_dyn = construct_node_dynamics(sw_eq)
        @test isa(swing_dyn, ODEVertex)
        @test swing_dyn.massmatrix == nothing

        pq_par = PQAlgebraic(S=S=10+3*im)
        pq_dyn = construct_node_dynamics(pq_par)
        @test pq_dyn.massmatrix == [0,0]

        
        #dx = [0.1 0.05 0.0]
        #x = [0.2 0.5 0.01]
        #e_s = [[2 0.3]]
        #e_d = [[1.5 0.4]]
        #p = []
        #t = 0.1
        #res = swing_dyn.f!(dx, x, e_s, e_d, p, t)
end
