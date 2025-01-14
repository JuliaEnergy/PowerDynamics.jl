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
net = Dict{String,Any}(
    "bus" => Dict{String,Any}(
        "1" => Dict{String,Any}(
            "name" => "bus_1",
            "bus_type" => 3,
            "va" => 0.0,
            "vm" => 1.0,
            "vmin" => 0.0,
            "base_kv" => 110.0,
            "vmax" => 1.05,
            "index" => 1,
        ),
        "2" => Dict{String,Any}(
            "name" => "bus_2",
            "bus_type" => 1,
            "va" => 0.0,
            "vm" => 1.0,
            "vmin" => 0.0,
            "base_kv" => 110.0,
            "vmax" => 1.05,
            "index" => 2,
        ),
    ),
    "branch" => Dict{String,Any}(
        "1" => Dict{String,Any}(
            "f_bus" => 1,
            "t_bus" => 2,
            "br_r" => 0.0,
            "br_x" => 0.01,
            "b_fr" => 0.0,
            "b_to" => 0.0,
            "br_status" => 1,
            "shift" => 0.0,
            "index" => 1,
            "name" => "Line",
            "angmin" => -30.0,
            "angmax" => 30.0,
            "transformer" => false,
            "tap" => 1.0,
            "g_to" => 0.0,
            "g_fr" => 0.0,
        ),
    ),
    "gen" => Dict{String,Any}(
        "1" => Dict{String,Any}(
            "name" => "Synchronous Machine",
            "index" => 1,
            "vg" => 1.0,
            "vbase" => 120.0,
            "mbase" => 200.0,
            "qg" => 0.0,
            "gen_status" => 1,
            "pg" => 0.0,
            "gen_bus" => 1,
            "dynamic_model" => Dict{String,Any}(
                "model_type" => "SixthOrderModel",
                "parameters" => Dict{String,Any}(
                    "ωs" => 1,
                    "H" => 1.3,
                    "Xl" => 0.172,
                    "Rs" => 0.0,
                    "Xq" => 2.0,
                    "Xq_d" => 0.3,
                    "Xq_dd" => 0.2,
                    "Xd" => 2.0,
                    "Xd_d" => 0.3,
                    "Xd_dd" => 0.2,
                    "Tq_d" => 6.66667,
                    "Td_dd" => 0.075,
                    "Td_d" => 6.66667,
                    "Tq_dd" => 0.075,
                    "consider_ωr_variations" => true,
                ),
                "controllers" => Dict{String,Any}(
                    "IEEET1" => Dict{String,Any}(
                        "model_type" => "IEEET1",
                        "parameters" => Dict{String,Any}(
                            "Te" => 0.2,
                            "Ta" => 0.03,
                            "Tf" => 1.5,
                            "Tr" => 0.02,
                            "Ke" => 1.0,
                            "Ka" => 200,
                            "Kf" => 0.05,
                            "E1" => 3.036,
                            "Se1" => 0.66,
                            "E2" => 4.048,
                            "Se2" => 0.88,
                            "Vrmin" => -10.0,
                            "Vrmax" => 10.0,
                        )
                    ),
                    "TGOV1" => Dict{String,Any}(
                        "model_type" => "TGOV1",
                        "parameters" => Dict{String,Any}(
                            "T1" => 0.5,
                            "T2" => 2.1,
                            "T3" => 7.0,
                            "Rd" => 0.05,
                            "Vmin" => 0.0,
                            "Vmax" => 1.0,
                        ),
                    ),
                ),
            ),
        ),
    ),
    "load" => Dict{String,Any}(
        "1" => Dict{String,Any}(
            "name" => "General Load",
            "load_bus" => 2,
            "status" => 1,
            "vm" => 0.01,
            "qd" => 0.03,
            "pd" => 0.5,
            "index" => 1,
            "dynamic_model" => Dict{String,Any}(
                "model_type" => "ZIPLoad",
                "parameters" => Dict{String,Any}(
                    "Kpz" => 1.0,
                    "Kpi" => 0.0,
                    "Kqz" => 1.0,
                    "Kqi" => 0.0,
                ),
            ),
        ),
    ),
    "shunt" => Dict{String,Any}(),
    "switch" => Dict{String,Any}(),
    "dcline" => Dict{String,Any}(),
    "storage" => Dict{String,Any}(),
)

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
        V_b = 110,
        Vn = dict["vbase"],
        ω_b = 2*π*50,
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

bdict = net["branch"]["1"]
@named pibranch = PiLine(;
    R=bdict["br_r"], X=bdict["br_x"],
    B_src = bdict["b_fr"], B_dst = bdict["b_to"],
    G_src = bdict["g_fr"], G_dst = bdict["g_to"],
    r_dst = bdict["tap"]
)
branch = Line(MTKLine(pibranch); src=bdict["f_bus"], dst=bdict["t_bus"])

busdata = net["bus"]["1"]
name = Symbol(replace(busdata["name"], " "=>""))
gen = gen_model(net["gen"]["1"])
genbus = Bus(MTKBus(gen); vidx=1, name, pf=pfSlack(;V=busdata["vm"], δ=busdata["va"]))

busdata = net["bus"]["2"]
name = Symbol(replace(busdata["name"], " "=>""))
load = load_model(net["load"]["1"])
loadbus = Bus(MTKBus(load); vidx=2, name, pf=pfPQ(;P=-net["load"]["1"]["pd"], Q=-net["load"]["1"]["qd"]))

nw = Network([genbus, loadbus], [branch])
OpPoDyn.solve_powerflow!(nw)
# set_default!(genbus, :machine_avr_gov₊avr₊vf₊u, 1.0973046159580935)
OpPoDyn.initialize!(nw)
@assert get_initial_state(genbus, :machine_avr_gov₊machine₊vf) > 0
@assert get_initial_state(genbus, :machine_avr_gov₊machine₊τ_m) > 0


function load_data(i)
    CSV.read(joinpath(pkgdir(OpPoDyn), "docs", "examples", "data", "rmspowersims_single_machine", "bus$i.csv"), DataFrame)
end

function load_ts(i, sym; f=x->x)
    df = load_data(i)
    collect(zip(df.t, f.(getproperty(df, sym))))
end


function affect(integrator)
    println("Change load at $(integrator.t)")
    s = NWState(integrator) # get indexable parameter object
    s.p.v[2, :load₊Pset] = -0.6
    auto_dt_reset!(integrator); save_parameters!(integrator)
end
cb = PresetTimeCallback(1.0, affect)

s0 = NWState(nw)
prob = ODEProblem(nw, uflat(s0), (0,10), pflat(s0); callback=cb, dtmax=1e-4)
sol = solve(prob, Rodas5P());


let
    fig = Figure()
    ax = Axis(fig[1,1])
    lines!(ax, load_ts(2, :Pd;f=x->-x))
    lines!(ax, sol; idxs=VIndex(2, :load₊P), color=Cycled(2))
    fig
end
let
    fig = Figure()
    ax = Axis(fig[1,1])
    lines!(ax, load_ts(1, :V))
    lines!(ax, sol; idxs=VIndex(1, :busbar₊u_mag), color=Cycled(2))
    ax = Axis(fig[2,1])
    lines!(ax, load_ts(2, :V))
    lines!(ax, sol; idxs=VIndex(2, :busbar₊u_mag), color=Cycled(2))
    fig
end

let
    fig = Figure()
    ax = Axis(fig[1,1]; title="Machine frequency")
    lines!(ax, load_ts(1, :ω))
    lines!(ax, sol; idxs=VIndex(1, :machine_avr_gov₊machine₊ω), color=Cycled(2))
    ax = Axis(fig[1,2]; title="Machine Ed")
    lines!(ax, load_ts(1, :Ed))
    lines!(ax, sol; idxs=VIndex(1, :machine_avr_gov₊machine₊E′_d), color=Cycled(2))
    ax = Axis(fig[1,3]; title="Machine Eq")
    lines!(ax, load_ts(1, :Eq))
    lines!(ax, sol; idxs=VIndex(1, :machine_avr_gov₊machine₊E′_q), color=Cycled(2))
    fig
end

# AVR
let
    fig = Figure()
    ax = Axis(fig[1,1]; title="vout")
    lines!(ax, load_ts(1, :Efd))
    lines!(ax, sol; idxs=VIndex(1, :machine_avr_gov₊avr₊vf₊u), color=Cycled(2))
    ax = Axis(fig[2,1]; title="vr (saturated)")
    lines!(ax, load_ts(1, :Vr))
    lines!(ax, sol; idxs=VIndex(1, :machine_avr_gov₊avr₊vr), color=Cycled(2))
    ax = Axis(fig[3,1]; title="vfb")
    lines!(ax, load_ts(1, :Vf))
    lines!(ax, sol; idxs=VIndex(1, :machine_avr_gov₊avr₊v_fb), color=Cycled(2))
    fig
end

names(load_data(1))
# GOV
let
    fig = Figure()
    ax = Axis(fig[1,1]; title="τ_m output")
    lines!(ax, load_ts(1, :Tm))
    lines!(ax, sol; idxs=VIndex(1, :machine_avr_gov₊gov₊τ_m₊u), color=Cycled(2))
    ax = Axis(fig[2,1]; title="valve (saturated)")
    lines!(ax, load_ts(1, :Pv))
    lines!(ax, sol; idxs=VIndex(1, :machine_avr_gov₊gov₊xg1), color=Cycled(2))
    ax = Axis(fig[3,1]; title="Pm")
    lines!(ax, load_ts(1, :Pm))
    lines!(ax, sol; idxs=VIndex(1, :machine_avr_gov₊gov₊xg2), color=Cycled(2))
    fig
end
