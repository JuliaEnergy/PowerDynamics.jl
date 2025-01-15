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

function load_data(i)
    CSV.read(joinpath(pkgdir(OpPoDyn), "docs", "examples", "data", "rmspowersims", "bus$i.csv"), DataFrame)
end
function load_ts(i, sym; f=x->x)
    df = load_data(i)
    collect(zip(df.t, f.(getproperty(df, sym))))
end

function load_model(dict)
    @assert dict["dynamic_model"]["model_type"] == "ZIPLoad"
    p = dict["dynamic_model"]["parameters"]
    ZIPLoad(; KpZ=p["Kpz"],KpI=p["Kpi"], KqZ=p["Kqz"], KqI=p["Kqi"],
            Pset=-dict["pd"], Qset=-dict["qd"], name=:load)
end

function gen_model(dict; S_b, V_b, ω_b)
    @assert dict["dynamic_model"]["model_type"] == "SixthOrderModel"
    controllers = dict["dynamic_model"]["controllers"]

    p = dict["dynamic_model"]["parameters"]
    @named machine = SauerPaiMachine(;
        R_s = p["Rs"],
        X_d = p["Xd"],
        X_q = p["Xq"],
        X′_d = p["Xd_d"],
        X′_q = p["Xq_d"],
        X″_d = p["Xd_dd"],
        X″_q = p["Xq_dd"],
        X_ls = p["Xl"],
        T′_d0 = p["Td_d"],
        T″_d0 = p["Td_dd"],
        T′_q0 = p["Tq_d"],
        T″_q0 = p["Tq_dd"],
        H = p["H"],
        S_b,
        Sn = dict["mbase"],
        V_b,
        ω_b,
    )
    controleqs = Equation[]
    if haskey(controllers, "IEEET1")
        avrp = controllers["IEEET1"]["parameters"]
        @named avr = AVRTypeI(
            ceiling_function=:quadratic,
            Ka = avrp["Ka"],
            Ke = avrp["Ke"],
            Kf = avrp["Kf"],
            Ta = avrp["Ta"],
            Tf = avrp["Tf"],
            Te = avrp["Te"],
            Tr = avrp["Tr"],
            vr_min = avrp["Vrmin"],
            vr_max = avrp["Vrmax"],
            E1 = avrp["E1"],
            Se1 = avrp["Se1"],
            E2 = avrp["E2"],
            Se2 = avrp["Se2"],
        )
        append!(controleqs, [connect(machine.v_mag_out, avr.vh), connect(avr.vf, machine.vf_in)])
    else
        @assert isempty(controllers)
        @named avr = AVRFixed()
        append!(controleqs, [connect(avr.vf, machine.vf_in)])
    end

    if haskey(controllers, "TGOV1")
        govp = controllers["TGOV1"]["parameters"]
        @named gov = TGOV1(;
            # p_ref = govp["P_set"],
            V_min = govp["Vmin"],
            V_max = govp["Vmax"],
            R = govp["Rd"],
            T1 = govp["T1"],
            T2 = govp["T2"],
            T3 = govp["T3"],
            DT = 0,
        )
        append!(controleqs, [connect(gov.τ_m, machine.τ_m_in), connect(machine.ωout, gov.ω_meas)])
    else
        @assert isempty(controllers)
        @named gov = GovFixed()
        append!(controleqs, [connect(gov.τ_m, machine.τ_m_in)])
    end

    CompositeInjector(
        [machine, avr, gov],
        controleqs,
        name=:ctrld_gen
    )
end

json = JSON.parsefile(joinpath(pkgdir(OpPoDyn), "docs", "examples", "data", "ieee39.json"),dicttype=OrderedDict)

function dict_to_df(dict; cols...)
    df = DataFrame(; cols...)
    for (k, v) in dict
        push!(df, v; cols=:union)
    end
    df
end
bus_df = sort!(dict_to_df(json["bus"]; index=Int[], name=String[], bus_type=Int[], base_kv=Float64[], vm=Float64[], va=Float64[]), :index)
gen_df = sort!(dict_to_df(json["gen"]; index=Int[], gen_bus=Int[], pg=Float64[], qg=Float64[], vg=Float64[], mbase=Float64[]), :index)
load_df = sort!(dict_to_df(json["load"]; index=Int[], load_bus=Int[], pd=Float64[], qd=Float64[], status=Int[]), :index)
branch_df = sort!(dict_to_df(json["branch"]; index=Int[], f_bus=Int[], t_bus=Int[], br_r=Float64[], br_x=Float64[], tap=Float64[], shift=Float64[], br_status=Int[], transformer=Int[], g_fr=Float64[], b_fr=Float64[], g_to=Float64[], b_to=Float64[]), :index)

machine_df = let
    machine_df = DataFrame()
    for (i, gendat) in enumerate(values(json["gen"]))
        p = gendat["dynamic_model"]["parameters"]
        p = OrderedDict(p)
        p["bus"] =  gendat["gen_bus"]
        push!(machine_df, p; cols=:union)
    end
    sort!(machine_df, "bus")
end

# set the correct powerflow models
leftjoin!(bus_df, select(load_df, :load_bus, :pd, :qd); on=:index=>:load_bus)
leftjoin!(bus_df, select(gen_df, :gen_bus, :pg, :qg, :vg); on=:index=>:gen_bus)
bus_df.ptotal = replace(bus_df.pg, missing=>0.0) - replace(bus_df.pd, missing=>0.0)
bus_df.qtotal = replace(bus_df.qg, missing=>0.0) - replace(bus_df.qd, missing=>0.0)

branches = let
    branches = []
    for (i, bdict) in enumerate(values(json["branch"]))
        println("Branch $i: $(bdict["f_bus"]) -> $(bdict["t_bus"])")

        tapkw = if bdict["transformer"] == true
            src_kv = bus_df[bdict["f_bus"], :base_kv]
            dst_kv = bus_df[bdict["t_bus"], :base_kv]
            if src_kv > dst_kv
                (; r_src = 1/bdict["tap"])
            elseif src_kv < dst_kv
                # TODO: looks like the tap is allways at src bus no matter where the high voltage side is
                (; r_src = 1/bdict["tap"])
                # (; r_dst = 1/bdict["tap"])
            end
        else
            (;)
        end

        @assert bdict["br_status"] == 1
        @assert bdict["shift"] == 0.0
        @named pibranch = PiLine(;
            R=bdict["br_r"], X=bdict["br_x"],
            B_src = bdict["b_fr"], B_dst = bdict["b_to"],
            G_src = bdict["g_fr"], G_dst = bdict["g_to"],
            tapkw...
        )
        line = Line(MTKLine(pibranch); src=bdict["f_bus"], dst=bdict["t_bus"])
        push!(branches, line)

    end
    branches
end;

@time busses = let
    busses = []
    for i in 1:39
        busdata = json["bus"][string(i)]
        gen_keys = findall(json["gen"]) do gendat
            gendat["gen_bus"] == i
        end
        load_keys = findall(json["load"]) do gendat
            gendat["load_bus"] == i
        end
        printstyled("Bus $(lpad(i,2)):  "; bold=true)
        print("$(length(gen_keys)) generators\t")


        println("$(length(load_keys)) loads")
        name = Symbol(replace(busdata["name"], " "=>""))

        S_b = json["baseMVA"]
        ω_b = json["dynamic_model_parameters"]["f_nom"]*2π
        V_b = busdata["base_kv"]

        if isempty(gen_keys) && isempty(load_keys)
            # empty bus
            bus = Bus(MTKBus(); vidx=i, name)
        elseif isempty(gen_keys) && length(load_keys) == 1
            # load bus
            load = load_model(json["load"][only(load_keys)])
            bus = Bus(MTKBus(load); vidx=i, name)
        elseif length(gen_keys) == 1 && isempty(load_keys)
            # pure generator
            gen = gen_model(json["gen"][only(gen_keys)]; S_b, ω_b, V_b)
            bus = Bus(MTKBus(gen); vidx=i, name)
        elseif length(gen_keys) == 1 && length(load_keys) == 1
            # generator and load
            load = load_model(json["load"][only(load_keys)])
            gen = gen_model(json["gen"][only(gen_keys)]; S_b, ω_b, V_b)
            bus = Bus(MTKBus(gen, load); vidx=i, name)
        else
            error()
        end
        push!(busses, bus)
    end
    busses
end;

nw = Network(copy.(busses), copy.(branches))
for row in eachrow(bus_df)
    busm = nw.im.vertexm[row.index]
    if row.bus_type == 1
        set_metadata!(busm, :pfmodel, pfPQ(P=row.ptotal, Q=row.qtotal))
    elseif row.bus_type == 2
        set_metadata!(busm, :pfmodel, pfPV(P=row.ptotal, V=row.vg))
    elseif row.bus_type == 3
        set_metadata!(busm, :pfmodel, pfSlack(;V=row.vm, δ=row.va))
    else
        error()
    end
end
OpPoDyn.solve_powerflow!(nw)

# compare powerflow
let
    df = OpPoDyn.show_powerflow(nw)
    V = Float64[]
    θ = Float64[]
    for i in 1:39
        _V, _θ = load_data(i)[1,[:V,:θ]]
        push!(V, _V)
        push!(θ, rad2deg(_θ))
    end
    df = DataFrame(N=1:39, V_nd=df."vm [pu]", V=V, θ_nd=df."varg [deg]", θ=θ)
    df.ΔV = df.V - df.V_nd
    df.Δθ = df.θ - df.θ_nd
    df
end

# Buses 31 and 39 have a load attached, we need to manualy initialize the Vset for those
set_default!(nw.im.vertexm[31], :load₊Vset, norm(get_initial_state.(Ref(nw.im.vertexm[31]), [:busbar₊u_r,:busbar₊u_i])))
set_default!(nw.im.vertexm[39], :load₊Vset, norm(get_initial_state.(Ref(nw.im.vertexm[39]), [:busbar₊u_r,:busbar₊u_i])))
OpPoDyn.initialize!(nw)

####
#### solve dyn system
####
function affect(integrator)
    println("Change something at $(integrator.t)")
    s = NWState(integrator) # get indexable parameter object
    s.p.v[16, :load₊Pset] *= 1.2
    # s.p.v[3, :load₊Pset] = 1000
    # s.p.e[1, :pibranch₊active] = 0
    # s.v[30, :ctrld_gen₊machine₊δ] += 0.001
    auto_dt_reset!(integrator); save_parameters!(integrator)
end
cb = PresetTimeCallback(1, affect)
s0 = NWState(nw)
prob = ODEProblem(nw, copy(uflat(s0)), (0,15), copy(pflat(s0)); callback=cb)
sol = solve(prob, Rodas5P());

# plot of fault
let
    fig = Figure()
    ax = Axis(fig[1,1])
    lines!(ax, load_ts(16, :Pd;f=x->-x))
    lines!(ax, sol, idxs=vidxs(sol, 16, :busbar₊P), color=Cycled(2))
    fig
end

# plot of generator states
function plot_gen_states(title, rmssym, ndsym)
    fig = Figure(size=(1000,800))
    row = 1
    col = 1
    for row in 1:3, col in 1:3
        i = 30 + 3*(row-1) + col
        ax = Axis(fig[row, col], title=title*" at bus $i")
        try lines!(ax, load_ts(i, rmssym)); catch; end
        try lines!(ax, sol; idxs=VIndex(i, ndsym), color=Cycled(2)); catch; end
    end
    fig
end
plot_gen_states("Voltage magnitude", :V, :busbar₊u_mag)
plot_gen_states("Machine frequency", :ω, :ctrld_gen₊machine₊ω)
plot_gen_states("Machine angle", :δ, :ctrld_gen₊machine₊δ)
plot_gen_states("Gov power", :Tm, :ctrld_gen₊machine₊τ_m)
plot_gen_states("Excitation voltage", :Efd, :ctrld_gen₊machine₊vf)
plot_gen_states("Regulator voltage (limited)", :Vr, :ctrld_gen₊avr₊vr)
plot_gen_states("Gov valve (limited)", :Pv, :ctrld_gen₊gov₊xg1)
