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

function load_model(dict)
    @assert dict["dynamic_model"]["model_type"] == "ZIPLoad"
    p = dict["dynamic_model"]["parameters"]
    ZIPLoad(; KpZ=p["Kpz"],KpI=p["Kpi"], KqZ=p["Kqz"], KqI=p["Kqi"],
            Pset=-dict["pd"], Qset=-dict["qd"], name=:load)
    # PQLoad(; Pset=-dict["pd"], Qset=-dict["qd"], name=:load)
    # Library.PQFiltLoad(; Pset=-dict["pd"], Qset=-dict["qd"], τ=0.1, name=:load)
end

function gen_model(dict)
    @assert dict["dynamic_model"]["model_type"] == "SixthOrderModel"
    controllers = dict["dynamic_model"]["controllers"]

    p = dict["dynamic_model"]["parameters"]
    @named machine = SauerPaiMachine(
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
        S_b = 100,
        Sn = dict["mbase"],
        V_b = dict["vbase"],
        ω_b = p["ωs"],
    )
    controleqs = Equation[]

    if haskey(controllers, "IEEET1")
        avrp = controllers["IEEET1"]["parameters"]
        # Ae, Be = Library.solve_ceilf(avrp["E1"]=>avrp["Se1"], avrp["E2"]=>avrp["Se2"]; u0=[1.0, 1.0])
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

    # @named avr = AVRFixed()
    # append!(controleqs, [connect(avr.vf, machine.vf_in)])
    # @named gov = GovFixed()
    # append!(controleqs, [connect(gov.τ_m, machine.τ_m_in)])

    CompositeInjector(
        [machine, avr, gov],
        controleqs,
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
        # bdict = collect(values(json["branch"]))[1]

        @assert bdict["br_status"] == 1
        @assert bdict["shift"] == 0.0
        @named pibranch = PiLine(;
            R=bdict["br_r"], X=bdict["br_x"],
            B_src = bdict["b_fr"], B_dst = bdict["b_to"],
            G_src = bdict["g_fr"], G_dst = bdict["g_to"],
            r_dst = bdict["tap"]
        )
        line = Line(MTKLine(pibranch); src=bdict["f_bus"], dst=bdict["t_bus"])
        push!(branches, line)
    end
    branches
end;

# sort!(gen_df, :gen_bus)

@time busses = let
    # slacks = [34, 35, 37, 38]
    slacks = []
    # slacks = collect(31:39)
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

        if i ∈ slacks
            bus = Bus(MTKBus(Library.UrUiConstraint(;name=:uconstraint)); vidx=i, name)
        elseif isempty(gen_keys) && isempty(load_keys)
            # empty bus
            bus = Bus(MTKBus(); vidx=i, name)
        elseif isempty(gen_keys) && length(load_keys) == 1
            # load bus
            load = load_model(json["load"][only(load_keys)])
            bus = Bus(MTKBus(load); vidx=i, name)
            set_default!(bus, :load₊Vset, busdata["vm"])
        elseif length(gen_keys) == 1 && isempty(load_keys)
            # pure generator
            gen = gen_model(json["gen"][only(gen_keys)])
            bus = Bus(MTKBus(gen); vidx=i, name)
        elseif length(gen_keys) == 1 && length(load_keys) == 1
            # generator and load
            load = load_model(json["load"][only(load_keys)])
            gen = gen_model(json["gen"][only(gen_keys)])
            bus = Bus(MTKBus(gen, load); vidx=i, name)
            set_default!(bus, :load₊Vset, busdata["vm"])
        else
            error()
        end
        set_voltage!(bus; mag=busdata["vm"], arg=busdata["va"])
        push!(busses, bus)
    end
    busses
end;

nw = Network(copy.(busses), copy.(branches))
for row in eachrow(bus_df)
    busm = nw.im.vertexm[row.index]
    if row.bus_type == 1
        # set_metadata!(busm, :pfmodel, pfPQ(P=row.ptotal, Q=row.qtotal))
    elseif row.bus_type == 2
        set_metadata!(busm, :pfmodel, pfPV(P=row.ptotal, V=row.vg))
    elseif row.bus_type == 3
        set_metadata!(busm, :pfmodel, pfSlack(;V=row.vm, δ=row.va))
    else
        error()
    end
end
OpPoDyn.solve_powerflow!(nw)
OpPoDyn.initialize!(nw)

function reinitialize!(cf)
    while OpPoDyn.get_initial_state(cf, :machine_avr_gov₊machine₊vf) < 0
        println("Try reinit...")
        for sym in cf.sym
            set_guess!(cf, sym, randn())
        end
        try
            initialize_component!(cf; verbose=false)
        catch e
        end
    end
end
for i in 30:39
    :machine_avr_gov₊machine₊vf ∈ obssym(nw.im.vertexm[i]) || continue
    println("Check reinitialization of $i")
    reinitialize!(nw.im.vertexm[i])
end

####
#### solve dyn system
####
function affect(integrator)
    println("Change something at $(integrator.t)")
    s = NWState(integrator) # get indexable parameter object
    # s.p.v[3, :load₊Pset] = 1000
    # s.p.e[1, :pibranch₊active] = 0
    # s.v[30, :machine_avr_gov₊machine₊δ] += 0.001
    auto_dt_reset!(integrator); save_parameters!(integrator)
end
cb = PresetTimeCallback(0.1, affect)

s0 = NWState(nw)
# s0.p.v[30, :machine_avr_gov₊machine₊D] .= 1
# s0.p.v[32:39, :machine_avr_gov₊machine₊D] .= 1
prob = ODEProblem(nw, uflat(s0), (0,100), pflat(s0); callback=cb, dtmax=1e-1)
sol = solve(prob, Rodas5P());

let
    i = 33
    fig = Figure(size=(1000,800))
    ax = Axis(fig[1, 1])
    lines!(ax, sol, idxs=vidxs(sol, i, :busbar₊u_arg); color=Cycled(1), label="u_arg")
    axislegend(ax; position=:lb)
    ax = Axis(fig[2,1])
    lines!(ax, sol, idxs=vidxs(sol, i, r"δ$"); color=Cycled(1), label="δ")
    axislegend(ax; position=:lb)
    ax = Axis(fig[3, 1])
    lines!(ax, sol, idxs=vidxs(sol, i, r"ω$"); color=Cycled(1), label="ω")
    axislegend(ax; position=:lb)
    ax = Axis(fig[4, 1])
    lines!(ax, sol, idxs=vidxs(sol, i, r"τ_m$"); color=Cycled(1), label="τ_m")
    lines!(ax, sol, idxs=vidxs(sol, i, r"τ_e$"); color=Cycled(2), label="τ_e")
    axislegend(ax; position=:rb)
    ax = Axis(fig[5, 1])
    # lines!(ax, sol, idxs=vidxs(sol, i, r"machine₊vf$"); color=Cycled(1), label="vf")
    lines!(ax, sol, idxs=vidxs(sol, i, r"avr₊vr$"); color=Cycled(1), label="vr")
    axislegend(ax; position=:rb)
    ax = Axis(fig[6,1])
    lines!(ax, sol, idxs=vidxs(sol, i, r"machine₊v_mag$"); color=Cycled(1), label="v_mag")
    axislegend(ax; position=:rb)
    ax = Axis(fig[7, 1])
    lines!(ax, sol, idxs=vidxs(sol, i, r"busbar₊P$"); color=Cycled(2), label="Bus P")
    axislegend(ax; position=:rb)
    fig
end

let
    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, sol, idxs=vidxs(sol, :, :busbar₊u_arg))
    ax = Axis(fig[1, 2])
    lines!(ax, sol, idxs=vidxs(sol, :, :busbar₊u_mag))
    # axislegend(ax; position=:lb)
    fig
end


let
    fig, ax, p = lines(sol, idxs=vidxs(sol, :, r"machine₊vf$"))
    # axislegend(ax; position=:lb)
    fig
end
fig, ax, p = lines(sol, idxs=vidxs(sol, :, r"load₊P$"))
fig, ax, p = lines(sol, idxs=vidxs(sol, :, :busbar₊u_mag))
# fig, ax, p = lines(sol, idxs=vidxs(sol, :, :busbar₊u_arg))

fig, ax, p = lines(sol, idxs=vidxs(sol, 30, r"δ$"))

lines!(sol, idxs=vidxs(sol, 30, r"busbar₊u_arg"))
fig

lines(sol, idxs=vidxs(sol, 3, :busbar₊P))
lines(sol, idxs=vidxs(sol, 3, :busbar₊u_mag))


####
#### check for anomalies
####
for (i, bus) in enumerate(nw.im.vertexm)
    :machine_avr_gov₊machine₊δ ∈ bus.sym || continue
    δ = get_initial_state(bus, :machine_avr_gov₊machine₊δ) |> normalize_angle
    u_arg = get_initial_state(bus, :busbar₊u_arg) |> normalize_angle
    # println("Bus $i: u_arg = $u_arg, δ = $δ, diff = $(normalize_angle(u_arg - δ))" )
    println(i, " ", get_initial_state(bus, :machine_avr_gov₊machine₊vf))
    # println(i, " ", get_initial_state(bus, :machine_avr_gov₊machine₊Q))
end
dump_state(nw.im.vertexm[30])
dump_state(nw.im.vertexm[32])
gen_df

cf = nw.im.vertexm[34]
cf
cf32a = copy(cf)

dump_state(cf32a)
dump_state(cf32b)

cf30a.sym .=> default_or_init_state(cf30a) - default_or_init_state(cf30b)


normalize_angle(get_initial_state(cf, :machine_avr_gov₊machine₊δ))
normalize_angle(get_initial_state(cf32b, :machine_avr_gov₊machine₊δ))

let
    δs = Float64[]
    v_ds = Float64[]
    v_qs = Float64[]
    for i in 1:2
        for sym in cf.sym
            set_guess!(cf, sym, rand())
        end
        initialize_component!(cf; verbose=false)
        if init_residual(cf) > 1e-8
            @warn "Encountered high residual"
        end
        # @show init_residual(cf)
        push!(δs, get_initial_state(cf, :machine_avr_gov₊machine₊δ))
        push!(v_ds, get_initial_state(cf, :machine_avr_gov₊machine₊E′_d))
        push!(v_qs, get_initial_state(cf, :machine_avr_gov₊machine₊E′_q))
    end

    fig = Figure(); ax = Axis(fig[1, 1])
    perm = sortperm(δs)
    ax.aspect = DataAspect()
    xlims!(-1.1,1.1)
    ylims!(-1.1,1.1)
    lines!(ax, [(0,0), (get_initial_state(cf, :busbar₊u_r), get_initial_state(cf, :busbar₊u_i))])
    scatter!(ax, cos.(δs[perm]), sin.(δs[perm]))
    scatter!(ax, v_ds[perm], v_qs[perm])
    fig
    nothing
    cf
end


# WHY IS IT EVER SO SLIGHTLY UNSTABLE?!?!
