using Test: @testset, @test
using PowerDynamics: GFI, construct_vertex, symbolsof, VSIVoltagePT1
using LinearAlgebra: I

@testset "Compare GFI_MTK vs. VSIVoltagePT1 Node" begin
    for i in 1:10
        para = Dict(
            :τ_v => rand(), # time constant voltage control delay
            :τ_P => rand(), # time constant active power measurement
            :τ_Q => rand(), # time constant reactive power measurement
            :K_P => rand(), # droop constant frequency droop
            :K_Q => rand(), # droop constant voltage droop
            :V_r => rand(), # reference/ desired voltage
            :P   => rand(), # active (real) power infeed
            :Q   => rand(), # reactive (imag) power infeed                .
            :ω_r => 0.0)    # refrence/ desired frequency

        ## create IONode
        nt = NamedTuple{Tuple(keys(para))}(values(para))
        node_bs = GFI(; nt...)
        f_bs = construct_vertex(node_bs).f

        ## create PDNode
        ## the PD node does not have the explicit ω_r parameter
        para_pd = delete!(copy(para), :ω_r)
        nt = NamedTuple{Tuple(keys(para_pd))}(values(para_pd))
        node_pd = VSIVoltagePT1(; nt...)

        f_pd = construct_vertex(node_pd).f

        ## create fake "edge data", 4 incoming, 4 outgooing with 4 states each
        edges = [randn(4) for i in 1:4]

        ## select random time
        t = rand()

        ## chose random initial state and account for initial ω in PD node
        x_bs = randn(4)
        x_pd = copy(x_bs)
        x_pd[3] = - para[:K_P] * (x_bs[3] - para[:P])

        ## create arrays for the results
        dx_bs = similar(x_bs)
        dx_pd = similar(x_pd)

        ## call both functions
        f_bs(dx_bs, x_bs, edges, nothing, t)
        f_pd(dx_pd, x_pd, edges, nothing, t)

        ## compare results
        ## we have to correct variable 3 of bs implementation to match dω
        @test dx_bs[1] ≈ dx_pd[1]                # u_r
        @test dx_bs[2] ≈ dx_pd[2]                # u_i
        @test - para[:K_P] * dx_bs[3] ≈ dx_pd[3] # ω
        @test dx_bs[4] ≈ dx_pd[4]                # q_filtered
    end

    τ_v, τ_P, τ_Q, K_P, K_Q = rand(5)
    P,Q,V_r = rand(3)
    q_m, dq_m, omega, domega = rand(4)
    gfi = GFI(τ_v=τ_v,τ_P=τ_P,τ_Q=τ_Q,K_P=K_P,K_Q=K_Q,V_r=V_r,P=P,Q=Q,ω_r=0.0)
    gfi_vertex = construct_vertex(gfi)
    dint = [domega,dq_m]; int = [omega,q_m]; int_test = copy(int)

    @test symbolsof(gfi) == [:u_r, :u_i, :p_filter₊output, :q_filter₊output]
    @test gfi.mass_matrix == I

    smoketest_rhs(gfi_vertex, int_x=[q_m, omega], int_dx=[dq_m, domega])
end
