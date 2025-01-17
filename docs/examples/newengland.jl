using OpPoDyn
using OpPoDyn.Library
using ModelingToolkit
using NetworkDynamics
using CairoMakie
using JSON
using OrderedCollections
using DataFrames
using OrdinaryDiffEqRosenbrock
using OrdinaryDiffEqNonlinearSolve
using DiffEqCallbacks
using CSV
using LinearAlgebra
using Test

json = JSON.parsefile(joinpath(pkgdir(OpPoDyn), "docs", "examples", "data", "ieee39.json"),dicttype=OrderedDict)

# define functions to load the "reference" timeseries from csv files
function load_data(i)
    CSV.read(joinpath(pkgdir(OpPoDyn), "docs", "examples", "data", "rmspowersims", "bus$i.csv"), DataFrame)
end
function load_ts(i, sym; f=x->x)
    df = load_data(i)
    collect(zip(df.t, f.(getproperty(df, sym))))
end

# helperfunciton to transform json dicts into dataframes
function dict_to_df(dict; cols...)
    df = DataFrame(; cols...)
    for (k, v) in dict
        push!(df, v; cols=:subset)
    end
    df
end

# define the main branch dataframe
branch_df = sort!(dict_to_df(json["branch"]; index=Int[], f_bus=Int[], t_bus=Int[], br_r=Float64[], br_x=Float64[], tap=Float64[], shift=Float64[], br_status=Int[], transformer=Int[], g_fr=Float64[], b_fr=Float64[], g_to=Float64[], b_to=Float64[]), :index)

# define commbinde bus dataframe which contains all the dynamical parameters
bus_df = let
    bus_df = sort!(dict_to_df(json["bus"]; index=Int[], bus_type=Int[], base_kv=Float64[], vm=Float64[], va=Float64[]), :index)

    # add loads data to bus_df
    _load_df = dict_to_df(json["load"]; load_bus=Int[], pd=Float64[], qd=Float64[])
    leftjoin!(bus_df, select(_load_df, :load_bus, :pd, :qd); on=:index=>:load_bus)

    # add gen data to bus_Df
    _gen_df = dict_to_df(json["gen"]; gen_bus=Int[], pg=Float64[], qg=Float64[], vg=Float64[], mbase=Float64[], vbase=Float64[])
    leftjoin!(bus_df, select(_gen_df, :gen_bus, :pg, :qg, :vg, :vbase, :mbase); on=:index=>:gen_bus)

    # add dynamic load data
    _dynload_df = DataFrame(;)
    for (i, loaddat) in enumerate(values(json["load"]))
        local p = loaddat["dynamic_model"]["parameters"]
        p = OrderedDict(p)
        p["bus"] =  loaddat["load_bus"]
        push!(_dynload_df, p; cols=:union)
    end
    leftjoin!(bus_df, _dynload_df, on=:index=>:bus)

    _machine_df = DataFrame(;)
    _avr_df = DataFrame(;)
    _gov_df = DataFrame(;)
    for (i, gendat) in enumerate(values(json["gen"]))
        # machine parameters
        local p = gendat["dynamic_model"]["parameters"]
        p = OrderedDict(p)
        p["bus"] =  gendat["gen_bus"]
        delete!(p, "index") # remove index key to avoid conflict on leftjoin
        push!(_machine_df, p; cols=:union)

        # avr parameters
        if haskey(gendat["dynamic_model"]["controllers"], "IEEET1")
            p = gendat["dynamic_model"]["controllers"]["IEEET1"]["parameters"]
            p["bus"] =  gendat["gen_bus"]
            push!(_avr_df, p; cols=:union)
        end
        # gov parameeters
        if haskey(gendat["dynamic_model"]["controllers"], "TGOV1")
            p = gendat["dynamic_model"]["controllers"]["TGOV1"]["parameters"]
            p["bus"] =  gendat["gen_bus"]
            push!(_gov_df, p; cols=:union)
        end
    end

    leftjoin!(bus_df, _machine_df; on=:index=>:bus)
    leftjoin!(bus_df, _avr_df; on=:index=>:bus)
    leftjoin!(bus_df, _gov_df; on=:index=>:bus)
    bus_df
end

# TODO: check vref/pset after init
_missing_to_zero(x) = ismissing(x) ? 0.0 : x
ptotal(i) = _missing_to_zero(bus_df[i, :pg]) - _missing_to_zero(bus_df[i, :pd])
qtotal(i) = _missing_to_zero(bus_df[i, :qg]) - _missing_to_zero(bus_df[i, :qd])
has_load(i) = !ismissing(bus_df[i, :pd])
has_gen(i) = !ismissing(bus_df[i, :pg])
has_avr(i) = !ismissing(bus_df[i, :Vref])
has_gov(i) = !ismissing(bus_df[i, :P_set])

@time branches = let
    branches = []
    for p in eachrow(branch_df)
        println("Branch: $(p.f_bus) -> $(p.t_bus)")

        tapkw = if p.transformer == true
            src_kv = bus_df[p.f_bus, :base_kv]
            dst_kv = bus_df[p.t_bus, :base_kv]
            # TODO: looks like the tap is allways at src bus no matter where the high voltage side is
            (; r_src = 1/p.tap)
            # if src_kv > dst_kv
            #     (; r_src = 1/p.tap)
            # elseif src_kv < dst_kv
            #     (; r_dst = 1/p.tap)
            # end
        else
            (;)
        end

        @named pibranch = PiLine(;
            R=p.br_r, X=p.br_x,
            B_src = p.b_fr, B_dst = p.b_to,
            G_src = p.g_fr, G_dst = p.g_to,
            tapkw...
        )
        line = Line(MTKLine(pibranch); src=p.f_bus, dst=p.t_bus)
        push!(branches, line)
    end
    branches
end;

@time busses = let
    busses = []
    for i in 1:39
        print("Bus $i:")

        local p = bus_df[i, :]
        components = [] # vector to collect dynamical components

        if has_load(i)
            @named load = ZIPLoad(;
                KpZ=p.Kpz,KpI=p.Kpi, KqZ=p.Kqz, KqI=p.Kqi,
                Pset=-p.pd, Qset=-p.qd
            )
            push!(components, load)
            print(" Load")
        end

        if has_gen(i)
            print(" Machine")
            controleqs = Equation[] # equations connecting avr and gov to machine
            @named machine = SauerPaiMachine(;
                R_s = p.Rs, X_ls = p.Xl,
                X_d = p.Xd, X_q = p.Xq,
                X′_d = p.Xd_d, X′_q = p.Xq_d,
                X″_d = p.Xd_dd, X″_q = p.Xq_dd,
                T′_d0 = p.Td_d, T′_q0 = p.Tq_d,
                T″_d0 = p.Td_dd, T″_q0 = p.Tq_dd,
                H = p.H,
                S_b = json["baseMVA"], Sn = p.mbase,
                V_b = p.base_kv, Vn = p.vbase,
                ω_b = json["dynamic_model_parameters"]["f_nom"]*2π,
            )

            # add dynamic or fixed avr
            if has_avr(i)
                print(" AVR")
                @named avr = AVRTypeI(
                    ceiling_function=:quadratic,
                    # vref = p.Vref # let this free for initialization
                    Ka = p.Ka, Ke = p.Ke, Kf = p.Kf,
                    Ta = p.Ta, Tf = p.Tf, Te = p.Te, Tr = p.Tr,
                    vr_min = p.Vrmin, vr_max = p.Vrmax,
                    E1 = p.E1, Se1 = p.Se1, E2 = p.E2, Se2 = p.Se2,
                )
                append!(controleqs, [connect(machine.v_mag_out, avr.vh), connect(avr.vf, machine.vf_in)])
            else
                @named avr = AVRFixed()
                append!(controleqs, [connect(avr.vf, machine.vf_in)])
            end

            # add dynamic or fixed gove
            if has_gov(i)
                print(" Gov")
                @named gov = TGOV1(;
                    # p_ref = p.P_set, # let this free for initialization
                    V_min = p.Vmin, V_max = p.Vmax, DT = 0,
                    R = p.Rd, T1 = p.T1, T2 = p.T2, T3 = p.T3,
                )
                append!(controleqs, [connect(gov.τ_m, machine.τ_m_in), connect(machine.ωout, gov.ω_meas)])
            else
                @named gov = GovFixed()
                append!(controleqs, [connect(gov.τ_m, machine.τ_m_in)])
            end

            comp = CompositeInjector(
                [machine, avr, gov],
                controleqs,
                name=:ctrld_gen
            )
            push!(components, comp)
        end

        # we also add the powerflow model based on the bus_type
        pfmodel = if p.bus_type == 1
            # PQ bus with S = S_gen - S_demand
            pfPQ(P=ptotal(i), Q=qtotal(i))
        elseif p.bus_type == 2
            # PV bus with P = P_gen - P_demand and V from (static) generator reference
            pfPV(P=ptotal(i), V=p.vg)
        elseif p.bus_type == 3
            # Slack bus with voltage from generator reference andangle 0
            pfSlack(;V=p.vg, δ=0)
        end

        bus = Bus(MTKBus(components...); vidx=i, pf=pfmodel, name=Symbol("bus$i"))
        push!(busses, bus)
        println()
    end
    busses
end;

nw = Network(copy.(busses), copy.(branches))
OpPoDyn.solve_powerflow!(nw)

# Buses 31 and 39 have a load attached, we need to manualy initialize the Vset for those
set_default!(nw.im.vertexm[31], :load₊Vset, norm(get_initial_state.(Ref(nw.im.vertexm[31]), [:busbar₊u_r,:busbar₊u_i])))
set_default!(nw.im.vertexm[39], :load₊Vset, norm(get_initial_state.(Ref(nw.im.vertexm[39]), [:busbar₊u_r,:busbar₊u_i])))
OpPoDyn.initialize!(nw)

@testset "Test initialization" begin
    # test powerflow results agains reference
    df = OpPoDyn.show_powerflow(nw)
    @test df."vm [pu]" ≈ bus_df.vm
    @test df."varg [deg]" ≈ rad2deg.(bus_df.va)
    # test initialized controller references Vref and P_set against reference
    @test bus_df[30:38, :Vref] ≈ get_initial_state.(nw.im.vertexm[30:38], :ctrld_gen₊avr₊vref)
    @test bus_df[30:38, :P_set] ≈ get_initial_state.(nw.im.vertexm[30:38], :ctrld_gen₊gov₊p_ref)
end

####
#### solve dyn system
####
function affect(integrator)
    println("Change something at $(integrator.t)")
    s = NWState(integrator) # get indexable parameter object
    s.p.v[16, :load₊Pset] *= 1.2
    auto_dt_reset!(integrator); save_parameters!(integrator)
end
cb = PresetTimeCallback(1, affect)
s0 = NWState(nw)
prob = ODEProblem(nw, copy(uflat(s0)), (0,15), copy(pflat(s0)); callback=cb)
sol = solve(prob, Rodas5P());
break

# plot he desired
let
    fig = Figure()
    ax = Axis(fig[1,1]; title="Power demand at bus 16", xlabel="Time [s]", ylabel="P [pu]")
    lines!(ax, load_ts(16, :Pd;f=x->-x))
    lines!(ax, sol, idxs=vidxs(sol, 16, :busbar₊P), color=Cycled(2))
    fig
end

# plot of generator states
function plot_gen_states(title, rmssym, ndsym)
    fig = Figure(size=(1000,800))
    row = 1; col = 1
    for row in 1:3, col in 1:3
        i = 30 + 3*(row-1) + col
        ax = Axis(fig[row, col], title=title*" at bus $i")
        try lines!(ax, load_ts(i, rmssym)); catch; end
        idxs = ndsym isa Symbol ? VIndex(i, ndsym) : ndsym(i)
        try lines!(ax, sol; idxs, color=Cycled(2)); catch; end
    end
    fig
end
plot_gen_states("Voltage magnitude", :V, :busbar₊u_mag)
plot_gen_states("Machine frequency", :ω, :ctrld_gen₊machine₊ω)
plot_gen_states("Machine angle", :δ, i->@obsex(VIndex(i, :ctrld_gen₊machine₊δ)-VIndex(39, :ctrld_gen₊machine₊δ)+load_data(39).δ[1]))
plot_gen_states("Mechanical Torque", :Tm, :ctrld_gen₊machine₊τ_m)
plot_gen_states("Excitation voltage (contains ceiling)", :Efd, :ctrld_gen₊machine₊vf)
plot_gen_states("Regulator voltage (limited)", :Vr, :ctrld_gen₊avr₊vr)
plot_gen_states("Gov valve (limited)", :Pv, :ctrld_gen₊gov₊xg1)
