using CSV
using DataFrames
using LightGraphs
using Plots
using PowerDynamics
using OrdinaryDiffEq
using HDF5
using Interpolations
using GLMakie, GraphMakie, ColorSchemes, NetworkLayout

## Load Network Structure and Parameters from Data

BusData = CSV.read("./examples/ieee14bus/Data/Bus.csv", DataFrame; types = [Int, Int, String, Int, Int])
LineData = CSV.read("./examples/ieee14bus/Data/Line.csv",DataFrame)
    LineData[!, :source] = LineData.source .|> Int
    LineData[!, :dest] = LineData.dest .|> Int
GeneratorData = CSV.read("./examples/ieee14bus/Data/Generator.csv",DataFrame)
LoadData = CSV.read("./examples/ieee14bus/Data/Load.csv",DataFrame)
FlowData = CSV.read("./examples/ieee14bus/Data/Flow.csv",DataFrame)

const N = nrow(BusData)
const L = nrow(LineData)
const G = nrow(GeneratorData)


g = SimpleGraph(N)
    for i = 1:L
        add_edge!(g, Int64(LineData[i, :source]), Int64(LineData[i, :dest]))
    end

node_df = outerjoin(BusData, GeneratorData, LoadData, FlowData, on=:ID, makeunique=true)

slack_idx = argmax(skipmissing(node_df.P_Gen))
#slack_idx = 13

## construct power grid

nodes = []
    for n in eachrow(node_df)
        if n.Number == slack_idx
            # in the data set, n.Vm is not exactly 1,
            # so our steady state will be slightly different
            push!(nodes, SlackAlgebraic(; U=complex(1.)) ) # n.Vm
        else
            if n.P_Gen |> ismissing
                #push!(nodes, PQAlgebraic(; P=-n.P_Load, Q=-n.Q_Load) )
                push!(nodes, PQAlgebraic(; P=-n.P_Load, Q=0.0) )
            else
                ### replace this with other node types ###
                push!(nodes, SwingEqLVS(; H=n.Inertia, P=n.P_Gen-n.P_Load, D=0.01, Ω=100π, Γ=10., V=1.) )
                #push!(nodes, SwingEq(; H=1., P=n.P_Gen-n.P_Load, D=0.1, Ω=100π) )
            end
        end
    end

lines = []
    for l in eachrow(LineData)
        if node_df[l.source, :Base_V] == node_df[l.dest, :Base_V] #normal line
            push!(lines, StaticLine(; from=l.source, to=l.dest, Y=inv(complex(l.r, l.x))))
        else
            #t_ratio = node_df[l.dest, :Base_V] / node_df[l.source, :Base_V]
            t_ratio = 1.0 # the tap is already included via the p.u. system
            push!(lines, Transformer(; from=l.source, to=l.dest, y=inv(complex(l.r, l.x)), t_ratio=t_ratio))
        end
    end

pg = PowerGrid(g, nodes, lines)

## power flow solution

data, result = power_flow(pg);
v = [result["solution"]["bus"][string(k)]["vm"] for k in 1:length(pg.nodes)]
va = [result["solution"]["bus"][string(k)]["va"] for k in 1:length(pg.nodes)]

## exchange slack with swing equation

Hs = node_df[slack_idx,:].Inertia
gen_idx = [gen for (gen,val) in data["gen"] if val["gen_bus"] == slack_idx][1]
P_Gen = result["solution"]["gen"][string(gen_idx)]["pg"]

SwEq = SwingEqLVS(; H=Hs, P=P_Gen, D=0.01, Ω=100π, Γ=10., V=1.)
nodes[slack_idx] = SwEq
pg = PowerGrid(g, nodes, lines)

## construct initial condition

ic_guess = PowerDynamics.initial_guess(pg, v .* exp.(1im .* va))
x0 = State(pg, ic_guess);

## simulate the system

tspan = (0.0, 100.0)
ω_idx = findall(:ω .∈ symbolsof.(pg.nodes));

ode = ODEProblem(rhs(pg), x0.vec, tspan);
dqsol = solve(ode, Rodas4());
sol = PowerGridSolution(dqsol, pg);

## plot results

Plots.plot(sol, :, :v, legend=false, xlim=(1e-3, 100.), xscale=:log10, ylabel="Vm")
Plots.plot(sol, :, :φ, legend=false, xlim=(1e-3, 100.), xscale=:log10, ylabel="Va")
Plots.plot(sol, ω_idx, :ω, legend=false, xlim=(1e-3, 100.), ylabel="\\omega")

## fluctuation time series

loads = [idx for (idx,val) in enumerate(pg.nodes) if typeof(val) == PQAlgebraic];
gens = [idx for (idx,val) in enumerate(pg.nodes) if typeof(val) == SwingEqLVS];

function time_interpolation(series,Δt)
  itp = interpolate(series, BSpline(Linear()))
  etp = extrapolate(itp, Interpolations.Flat()) #extrapolate first and last value to all times
  t -> etp(t/Δt .+1)
end

series = h5open("./examples/ieee14bus/fluctuations_time_series.h5","r") do file
    read(file, "sum_of_fluctuations")
end

Fluct = time_interpolation(series,0.01)

## simulate flucutation

t_idx = Float64[]

for bus in loads

    acp = nodes[bus].P
    rep = nodes[bus].Q

    nodes[bus] = FluctuationNode(t -> acp + 0.2*Fluct(t), t -> rep)
    pg_dynamic = PowerGrid(nodes, lines)

    timespan = (-5.,50.)
    #fault = ChangeInitialConditions(node="bus2", var=:ω, f=Inc(0.0)) # dummy fault
    #sol = simulate(fault, pg_dynamic, x0.vec, timespan)
    ode = ODEProblem(rhs(pg), x0.vec, tspan)
    dqsol = solve(ode, Rodas4())
    sol = PowerGridSolution(dqsol, pg);
    println("Simulate fluctuation at ",bus, ": ", sol.dqsol.retcode)
    
    Δt = 0.01; T = 50
    norm = sum(sol(0.0:Δt:T,gens,:ω).^2)*Δt / T / length(gens) # 1/N ∑(1/T ∫ ωᵢ² dt)
    push!(t_idx,norm)

    nodes[bus] = PQAlgebraic(P=acp,Q=rep)

end

## plot fluctuation

include("plotexample.jl")
dg, nlabels, arrow_size, node_color, node_shape = flow_plot(x0,t_idx)
fig = Figure()
ax = Axis(fig[1,1]);
graphplot!(dg, node_color=node_color,node_size=15,node_marker=node_shape,layout=Stress())
fig[1,2] = Colorbar(fig; vertical=true, width=20,
                    limits = (minimum(t_idx),maximum(t_idx)),
                    colormap=:copper, label="troublemaker index")
#save("rts96.png", fig)