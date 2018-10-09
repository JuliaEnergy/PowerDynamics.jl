begin
        using PowerDynBase
        # using SymEngine
        using SymPy
        Core.eval(PowerDynBase, :(using SymPy))
        Core.eval(PowerDynBase, :(pi = PI))
        using MacroTools
        using Test
        using Random
        @show random_seed = 1234
        Random.seed!(random_seed)
        using LinearAlgebra
end

println("# Definition of helper functions #")
begin
    # within construct_network_rhs, the array of real values (necessary for DifferentialEqautions.jl)
    # is converted to an array of complex values, mostly for the complex voltage.
    # As this is a small hack, we need to prepare here that SymPy Arrays can do the same.
    # In the following, we define something that emulates a reinterpreted vector.
    # Then we overload complex_view to use this emulated version instead
    # of the real reinterpret as done in GridDynamics.jl.
    using Base.Iterators: flatten, partition
    import Base: setindex!, size, getindex, IndexStyle
    import PowerDynBase: complexview

    struct ComplexTestVector{T} <: AbstractVector{T}
        vec::AbstractVector{T}
        function ComplexTestVector(vec::AbstractVector{T}) where T
            @assert length(vec) % 2 == 0 "need an even length to interpret it as an array of complex numbers"
            new{T}(vec)
        end
    end
    realindex(i::Integer) = 2*(i-1)+1
    imagindex(i::Integer) = realindex(i)+1
    function Base.setindex!(A::ComplexTestVector{T}, v, I::Vararg{Int}) where {T}
        (i, ) = I
        A.vec[[realindex(i), imagindex(i)]] = [real(v), imag(v)]
     end
    function getindex(v::ComplexTestVector{T}, i::Integer) where T
        i_real = 2*(i-1)+1
        i_imag = i_real + 1
        v.vec[i_real]+ v.vec[i_imag]*im
    end
    size(v::ComplexTestVector{T}) where T = (div(length(v.vec), 2),)
    Base.IndexStyle(::Type{<:ComplexTestVector{T}}) where T = IndexLinear()
    function complexview(vec::Array{Sym,1}, begin_index, num_complex_numbers)
        begin_index_original = 2*(begin_index-1)+1
        vec_view = view(vec, begin_index_original:begin_index_original + 2*(num_complex_numbers-1)+1)
        ComplexTestVector(vec_view)
    end

    # and two helper functions
    "Split an array of (complex) symbols into their real and complex parts."
    splitcomplexsymbols(vec::AbstractVector{Sym}) =
        collect(Base.Iterators.flatten([func(el) for func in (real, imag), el in vec]))
    "Merge an array of (real) symbols into (half the number of) complex symbols."
    mergecomplexsymbols(vec::AbstractVector{Sym}) =
        [vec[i]+im*vec[i+1] for i in 1:2:length(vec)]
    "Create a matrix of SymPy symbols."
    SymbolMatrix(var::AbstractString, n, m; kwargs...) = [symbols("$(var)_$i$j"; kwargs...) for i in 1:n, j in 1:m]
    SymbolMatrix(var::AbstractString, n; kwargs...) = SymbolMatrix(var, n, n; kwargs...)
end

println("#### Single Node Tests ####")
@time @testset "Single Node Tests" begin
    ################################################################################
    # define variables as SymPy symbols for all node tests
    begin
        @syms u i
        @syms t real=true
        # derivaed dynamic variables
        v = abs(u)
        s = u*conj(i)
        p, q = real(s), imag(s)
    end
    ################################################################################

    # PowerDynBase.jl is only assembling and exporting stuff, nothing to be tested
    # NodeDynamicsBase.jl is test with all files in the NodeDynamics directory (see below)

    let
        println("## NodeDynamics/PQAlgebraic.jl ##")
        @syms S
        pq_par = PQAlgebraic(S=S)
        pq_dyn = construct_node_dynamics(pq_par)
        @test pq_par === parametersof(pq_dyn)
        dint = []; int = []
        @test pq_dyn.ode_dynamics.rhs(dint, u, i, int, t) == S - s
        @test internalsymbolsof(pq_dyn) == []
        @test internaldsymbolsof(pq_dyn) == []

        println("## Functor of OrdinaryNodeDynamicsWithMass ##")
        @syms du
        ints = PowerDynBase.ODEVariable(val = int)
        us = PowerDynBase.ODEVariable(val = [u])
        pq_dyn(1, us, i, ints, t)
        @test us.ddt[1] == S - s
    end

    let
        println("## NodeDynamics/PVAlgebraic.jl ##")
        @syms P real=true
        @syms V positive=true
        pv_dyn = construct_node_dynamics(PVAlgebraic(P=P, V=V))
        dint = []; int = []
        @test pv_dyn.ode_dynamics.rhs(dint, u, i, int, t) == (v-V) + im*(p-P)
        @test internalsymbolsof(pv_dyn) == []
        @test internaldsymbolsof(pv_dyn) == []
    end

    let
        println("## NodeDynamics/SlackAlgebraic.jl ##")
        @syms U
        slack_dyn = construct_node_dynamics(SlackAlgebraic(U=U))
        dint = []; int = []
        @test slack_dyn.ode_dynamics.rhs(dint, u, i, int, t) == u - U
        @test internalsymbolsof(slack_dyn) == []
        @test internaldsymbolsof(slack_dyn) == []
    end

    let
        println("## NodeDynamics/SwingEquation.jl ##")
        println("# SwingEq #")
        @syms H D positive=true
        @syms P Ω real=true
        @syms omega domega real=true
        swing_par = SwingEq(H=H, P=P, D=D, Ω=Ω)
        swing_dyn = construct_node_dynamics(swing_par)
        @test swing_par === parametersof(swing_dyn)
        dint = [domega]; int = [omega]; int_test = copy(int)
        @test swing_dyn.rhs(dint, u, i, int, t) == u*im*omega
        @test expand.(dint) == expand.([(P - D*omega - p)*2PI*Ω/H])
        @test internalsymbolsof(swing_dyn) == [:ω]
        @test internaldsymbolsof(swing_dyn) == [:dω]

        println("## Functor of OrdinaryNodeDynamics ##")
        @syms du
        ints = PowerDynBase.ODEVariable(val = int, ddt = dint)
        us = PowerDynBase.ODEVariable(val = [u], ddt = [du])
        swing_dyn(1, us, i, ints, t)
        @test us.ddt[1] == u*im*omega
        @test expand.(ints.ddt) == expand.([(P - D*omega - p)*2PI*Ω/H])

        println("# SwingEqLVS #")
        @syms V Γ positive=true
        swing_lvs_dyn = construct_node_dynamics(SwingEqLVS(H=H, P=P, D=D, Ω=Ω, Γ=Γ, V=V))
        dint = [domega]; int = [omega]; int_test = copy(int)
        @test swing_lvs_dyn.rhs(dint, u, i, int, t) == u*im*omega - u/v * Γ * (v-V)
        @test expand.(dint) == expand.([(P - D*omega - p)*2PI*Ω/H])
        @test internalsymbolsof(swing_lvs_dyn) == [:ω]
        @test internaldsymbolsof(swing_lvs_dyn) == [:dω]

        println("# SwingEq as PowerDynBase.AlgebraicNodeDynamics #")
        algebraic_swing_dyn = convert(PowerDynBase.AlgebraicNodeDynamics, swing_dyn)
        @test swing_par === parametersof(algebraic_swing_dyn)
        @syms out_omega real=true
        @syms du
        int_out = [out_omega]; dint = [domega]; int = [omega]
        @test algebraic_swing_dyn.root(int_out, du,  dint, u, i, int, t) == u*im*omega - du
        @test expand.(int_out) == expand.([(P - D*omega - p)*2PI*Ω/H]) .- dint
        @test internalsymbolsof(algebraic_swing_dyn) == [:ω]
        @test internaldsymbolsof(algebraic_swing_dyn) == [:dω]
        @test internaloutsymbolsof(algebraic_swing_dyn) == [:outω] # automatically generated symbol
    end

    let
        @syms H  D  Ω  T_d_dash T_q_dash X_q_dash X_d_dash X_d X_q positive=true
        @syms P  E_f real=true
        @syms omega domega real=true
        fourth_dyn = construct_node_dynamics(FourthEq(H=H, P=P, D=D, Ω=Ω, E_f=E_f, T_d_dash=T_d_dash ,T_q_dash=T_q_dash ,X_q_dash=X_q_dash ,X_d_dash=X_d_dash,X_d=X_d, X_q=X_q))
        dint = [domega]; int = [omega]; int_test = copy(int)
        du = fourth_dyn.rhs(dint, u, i, int, t)
        @test real(du) == (1 / T_q_dash)* (- real(u) + Ω * 2PI / H*(X_q - X_q_dash)* imag(i))
        @test imag(du) == (1 / T_d_dash)* (- imag(u) - Ω * 2PI / H*(X_d - X_d_dash) * real(i) + E_f)
        @test expand.(dint) == expand.([(P - D*omega - p)*2PI*Ω/H])

    end
end


println("#### Grid Construction Tests ####")
println("## GridDynamics.jl ##")
println("# OrdinaryGridDynamics #")
@time @testset "OrdinaryGridDynamics" begin
        @syms t real=true
        num = 2
        pars = [SwingEq(H=symbols("H_$i", positive=true), P=symbols("P_$i", real=true), D=symbols("D_$i", positive=true), Ω=symbols("Omega_$i", real=true)) for i=1:num ]
        dyns = construct_node_dynamics.(pars)
        # LY = [symbols("LY_$i$j") for i in 1:num, j in 1:num]
        LY = SymbolMatrix("LY", num)
        grid_dyn = GridDynamics(pars, LY, skip_LY_check=true)
        @test grid_dyn isa PowerDynBase.OrdinaryGridDynamics
        @test pars == grid_dyn |> Nodes .|> parametersof

        # evalute the grid rhs symbolically
        us = [symbols("u_$i") for i in 1:num ]
        omegas = [symbols("omega_$i",real=true) for i=1:num ]
        xs = [splitcomplexsymbols(us); omegas]
        dxs = [
                [symbols("d$u")[1] for u in us]|> splitcomplexsymbols;
                [symbols("d$omega") for omega in omegas]
                ]
        dxs_test = copy(dxs)
        grid_dyn(dxs, xs, nothing, t)
        dus = dxs[1:2*num] |> mergecomplexsymbols
        domegas = dxs[2*num+1:end]

        # Rebuild the grid rhs with the assumption that the node dynamics is working
        # already (as it was tested separately before).
        us = [symbols("u_$i") for i in 1:num ] .|> complex
        us_var = PowerDynBase.ODEVariable(us)
        omegas = [symbols("omega_$i",real=true) for i=1:num ]
        omegas_var = PowerDynBase.ODEVariable(omegas)
        t_i = LY*us
        map(i -> dyns[i](i, us_var, t_i, view(omegas_var, i:i), nothing), 1:num)
        @test all(complex.(dus) .== complex.(us_var.ddt)) # check whether the voltage differentials match
        @test all(complex.(domegas) .== complex.(omegas_var.ddt)) # check whether the internal differentials match
end


println("# OrdinaryGridDynamicsWithMass #")
@time @testset "OrdinaryGridDynamicsWithMass" begin
        @syms t real=true
        swing_num = 2
        pq_num = 2
        num = pq_num + swing_num
        pars = [
        [SwingEq(H=symbols("H_$i", positive=true), P=symbols("P_$i", real=true), D=symbols("D_$i", positive=true), Ω=symbols("Omega_$i", real=true)) for i=1:swing_num ];
        [PQAlgebraic(S=symbols("S_$i")) for i in 1:pq_num]
        ]
        dyns = construct_node_dynamics.(pars)
        LY = SymbolMatrix("LY", num)
        grid_dyn = GridDynamics(pars, LY, skip_LY_check=true)
        @test typeof(grid_dyn) === PowerDynBase.OrdinaryGridDynamicsWithMass
        @test PowerDynBase.masses(grid_dyn) == [true, true, true, true, false, false, false, false, true, true]
        @test pars == grid_dyn |> Nodes .|> parametersof

        # evalute the grid rhs symbolically
        dyns = convert(Array{OrdinaryNodeDynamicsWithMass}, dyns)
        us = [symbols("u_$i") for i in 1:num ]
        omegas = [symbols("omega_$i",real=true) for i=1:swing_num ]
        xs = [splitcomplexsymbols(us); omegas]
        dxs = [
        [symbols("d$u")[1] for u in us]|> splitcomplexsymbols;
        [symbols("d$omega") for omega in omegas]
        ]
        dxs_test = copy(dxs)
        grid_dyn(dxs, xs, nothing, t)
        dus = dxs[1:2*num] |> mergecomplexsymbols
        domegas = dxs[2*num+1:end]

        # Rebuild the grid rhs with the assumption that the node dynamics is working
        # already (as it was tested separately before).
        us = [symbols("u_$i") for i in 1:num ] .|> complex
        us_var = PowerDynBase.ODEVariable(us)
        omegas = [symbols("omega_$i",real=true) for i=1:swing_num ]
        omegas_var = PowerDynBase.ODEVariable(omegas)
        t_i = LY*us
        map(i -> dyns[i](i, us_var, t_i, view(omegas_var, i:i), nothing), 1:swing_num);
        map(i -> dyns[i](i, us_var, t_i, view(omegas_var,1:0), nothing), swing_num + 1:num);
        @test all(complex.(dus) .== complex.(us_var.ddt)) # check whether the voltage differentials match
        @test all(complex.(domegas) .== complex.(omegas_var.ddt)) # check whether the internal differentials match
end

println("# complexview #")
@time @testset "complexview" begin
        function cmp_complex_real(c_arr, r_arr)
                cmp_arr = Array{Bool}(undef, length(c_arr))
                for i in 1:length(c_arr)
                        cmp_arr[i] = (c_arr[i] == r_arr[2*i-1] + im*r_arr[2*i])
                end
                all(cmp_arr)
        end
        num_real = 21 # should be an odd number
        @assert isodd(num_real) "test for odd numbers, too"
        num_complex = 3
        @assert 2*num_complex <= num_real "need to get less complex number than half the real numbers"
        real_array = rand(num_real)
        complex_array = PowerDynBase.complexview(real_array, 1, 3)
        @test length(complex_array) == 3
        @test cmp_complex_real(complex_array, real_array)
        complex_array[:] = rand(3) .+ im.*rand(3)
        @test cmp_complex_real(complex_array, real_array)
        real_array[:] = rand(num_real)
        @test cmp_complex_real(complex_array, real_array)
end


println("# OrdinaryGridDynamicsWithMass (non-symbolically) #")
@time @testset "OrdinaryGridDynamicsWithMass (non-symbolically)" begin
        nodes = [SwingEq(H=1, P=1, D=1, Ω=50), SwingEq(H=1, P=-1, D=1, Ω=50)]
        LY = [im -im; -im im]
        grid = GridDynamics(nodes, LY)
        x = rand(SystemSize(grid))
        dx = similar(x)
        grid(dx, x, nothing, 0)
        @test true # if the code runs until here, everything succeeded
end


println("#### Testing Output and Errors ####")
struct DummyNodeDynamics{N <: PowerDynBase.AbstractNodeParameters} <: PowerDynBase.AbstractNodeDynamics{N} end
@time @testset "Output & Errors" begin
    println("## show ##")
    @test_nowarn println(PQAlgebraic)
    @test_nowarn println(SwingEq)
    @test_nowarn println(PQAlgebraic(S=3+4im))
    @test_nowarn println(SwingEq(H=1, P=2, D=3, Ω=4))
    println("# Error Output #")
    Base.showerror(stdout, NodeDynamicsError("My message."))
    println()
    Base.showerror(stdout, GridDynamicsError("My message."))
    println()
    println("# check_admittance_matrix #")
    LY = [1 -1; 0 0]
    @test_throws GridDynamicsError PowerDynBase.checkLY(LY)
    @test_throws GridDynamicsError PowerDynBase.checkLY(transpose(LY))

    LY = [1 -1; -1 1] # now a valid admittance laplacian
    @test_throws GridDynamicsError GridDynamics([DummyNodeDynamics{PowerDynBase.SwingEq}()], LY)

    println("# UndefKeywordError #")
    @test_throws UndefKeywordError SlackAlgebraic()
end

println("#### Testing DynamicNode Macro ####")
@time @testset "@DynamicNode Macro" begin
        begin
            _removelinereferences(b) = b
            function _removelinereferences(b::Expr)
                if ~(b.head in (:block, :macrocall))
                    return b
                end
                b.args = [ el for el in b.args if ~(el isa LineNumberNode) ]
                b
            end
            removelinereferences(q::Expr) = MacroTools.postwalk( _removelinereferences , q)
            rlr = removelinereferences
            _removesingleblocks(b) = b
            function _removesingleblocks(b::Expr)
                if ~(b.head in (:block, )) || length(b.args) != 1
                    return b
                end
                b.args[1]
            end
            removesingleblocks(q::Expr) = MacroTools.postwalk( _removesingleblocks , q)
            rsb = removesingleblocks
        end
        struct_def = PowerDynBase.buildparameterstruct(:MyParams, [:p1, :p2])
        struct_def_tests = :(struct MyParams <: AbstractNodeParameters
                p1
                p2
                MyParams(; p1, p2) = new(p1, p2)
                end)
        @test struct_def |> rlr == struct_def_tests |> rlr |> rsb
        internals = PowerDynBase.getinternalvars(Val{:OrdinaryNodeDynamics}, :([[a, da], [b, db]]))
        @test internals == (vars=[:a, :b], dvars=[:da, :db])
        @test PowerDynBase.getinternalvars(Val{:OrdinaryNodeDynamicsWithMass}, :([[a, da], [b, db]])) == internals
        @test_throws NodeDynamicsError PowerDynBase.getinternalvars(Val{:unknown}, :([[a, da], [b, db]]))
        dynamicscall = :(OrdinaryNodeDynamicsWithMass(a=1, b=3))
        funcbody = quote
            funcbodyline1
            funcbodyline2
            end
        prep = quote
                prepline1
                prepline2
            end
        cndfunction = quote
                function construct_node_dynamics(par::MyParams)
                        p1 = par.p1
                        p2 = par.p2
                        prepline1
                        prepline2
                end
        end |> rlr |> rsb
        PowerDynBase.cndfunction_builder!(Val{:OrdinaryNodeDynamicsWithMass},
                internals,
                dynamicscall,
                funcbody,
                cndfunction)
        cnd_function_test = quote
                function construct_node_dynamics(par::MyParams)
                        p1 = par.p1
                        p2 = par.p2
                        prepline1
                        prepline2
                        function rhs!(dint::AbstractVector, u, i, int::AbstractVector, t)
                                a = int[1]
                                b = int[2]
                                funcbodyline1
                                funcbodyline2
                                try
                                        dint[1] = da
                                        dint[2] = db
                                        return du
                                catch e
                                        if typeof(e) === UndefVarError
                                                throw(NodeDynamicsError("you need to provide $(e.var)"))
                                        else
                                                throw(e)
                                        end
                                end
                        end
                        OrdinaryNodeDynamicsWithMass(a=1, b=3, rhs=rhs!, n_int=2, symbols=ODENodeSymbols([:a, :b], [:da, :db]), parameters=par)
                end
        end |> rlr |> rsb
        @test cndfunction |> rlr |> rsb == cnd_function_test
        full_macro_return = PowerDynBase.DynamicNode(:(MyParams(p1, p2) <: OrdinaryNodeDynamicsWithMass(a=1, b=3)), prep, :([[a, da], [b, db]]), funcbody ) |> rlr |> rsb
        full_macro_return_test = quote
                @__doc__ $(struct_def_tests)
                $(cnd_function_test)
        end |> rlr |> rsb
        @test full_macro_return == full_macro_return_test
end

println("#### States ####")
@time @testset "States" begin
    nodes = [SlackAlgebraic(U=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1), SwingEqLVS(H=1, P=-1, D=1, Ω=50, Γ=20, V=1)]
    numnodes = length(nodes)
    LY = SymbolMatrix("LY", numnodes)
    grid = GridDynamics(nodes, LY, skip_LY_check=true)
    @syms u_Sl u_Sw1 u_Sw2
    @syms omega1 omega2 real=true
    state = State(grid, [real(u_Sl), imag(u_Sl), real(u_Sw1), imag(u_Sw1), real(u_Sw2), imag(u_Sw2), omega1, omega2])
    println("## getindex ##")
    @test state[1, :u] == complex(u_Sl)
    @test state[2, :u] == complex(u_Sw1)
    @test state[1:2, :u] == [complex(u_Sl), complex(u_Sw1)]
    @test state[:, :u] == [complex(u_Sl), complex(u_Sw1), complex(u_Sw2)]
    @test state[2, :v] == complex(u_Sw1) |> abs
    @test state[3, :v] == complex(u_Sw2) |> abs
    @test state[2, :int, 1] == omega1
    @test state[2, :ω] == omega1
    @test state[3, :int, 1] == omega2
    @test state[3, :ω] == omega2
    @test state[2:3, :int, 1] == [omega1, omega2]
    @test state[2:3, :ω] == [omega1, omega2]
    @test state[2:3, :int, [1, 1]] == [omega1, omega2]
    @test state[1, :i] == (LY * complex.([u_Sl, u_Sw1, u_Sw2]))[1]
    @test state[:, :i] == LY * complex.([u_Sl, u_Sw1, u_Sw2])
    @test state[1, :iabs] == abs((LY * complex.([u_Sl, u_Sw1, u_Sw2]))[1])
    @test state[:, :iabs] == abs.((LY * complex.([u_Sl, u_Sw1, u_Sw2])))
    @test state[2, :s] ==  state[2, :u] * conj((LY * complex.([u_Sl, u_Sw1, u_Sw2]))[2])
    @test state[2, :p] ==  real(state[2, :u] * conj((LY * complex.([u_Sl, u_Sw1, u_Sw2]))[2]))
    @test state[2, :q] ==  imag(state[2, :u] * conj((LY * complex.([u_Sl, u_Sw1, u_Sw2]))[2]))
    @test state[:, :s] ==  state[:, :u] .* conj.((LY * complex.([u_Sl, u_Sw1, u_Sw2])))
    @test state[:, :p] ==  real(state[:, :u] .* conj.((LY * complex.([u_Sl, u_Sw1, u_Sw2]))))
    @test state[:, :q] ==  imag(state[:, :u] .* conj.((LY * complex.([u_Sl, u_Sw1, u_Sw2]))))

    @test_throws BoundsError state[numnodes+1, :u]
    @test_throws BoundsError state[2, :int, 2]
    @test_throws BoundsError state[1, :int, 1]
    @test_throws BoundsError state[numnodes+1, :int, 1]
    @test_throws StateError state[1, :ω]
    println("## setindex! ##")
    # exchange the u as test
    state[1, :u] = u_Sw1
    @test state[1, :u] == complex(u_Sw1)
    state[2, :u] = u_Sl
    @test state[2, :u] == complex(u_Sl)
    # change back
    state[:, :u] = [u_Sl, u_Sw1, u_Sw2]
    @test state[:, :u] == [complex(u_Sl), complex(u_Sw1), complex(u_Sw2)]
    @test_throws BoundsError state[numnodes+1, :u] = u_Sw1
    # exhange omega
    state[2, :int, 1] = omega2
    @test state[2:3, :int, 1] == [omega2, omega2]
    state[2:3, :int, 1] = [omega1, omega2]
    @test state[2:3, :int, 1] == [omega1, omega2]
    state[2, :ω] = omega2
    @test state[2:3, :ω] == [omega2, omega2]
    state[2:3, :ω] = [omega1, omega2]
    @test state[2:3, :ω] == [omega1, omega2]
    @test_throws BoundsError state[1, :int, 1] = omega2
    @test_throws StateError state[1, :ω] = omega2
    println("## getindex for angle (numerically) ##")
    v_Sl = rand()
    v_Sw1 = rand()
    v_Sw2 = rand()
    φ_Sl = rand()
    φ_Sw1 = rand()
    φ_Sw2 = rand()
    u_Sl = v_Sl*exp(im*φ_Sl)
    u_Sw1 = v_Sw1*exp(im*φ_Sw1)
    u_Sw2 = v_Sw2*exp(im*φ_Sw2)
    omega1 = 3.
    omega2 = 5.
    state = State(grid, [real(u_Sl), imag(u_Sl), real(u_Sw1), imag(u_Sw1), real(u_Sw2), imag(u_Sw2), omega1, omega2])
    @test state[1, :φ] ≈ φ_Sl
    @test state[2, :φ] ≈ φ_Sw1
    @test state[3, :φ] ≈ φ_Sw2
end
