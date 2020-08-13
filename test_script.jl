
using Pkg;
  Pkg.instantiate();
  cd(@__DIR__);
  using PowerDynamics;
  using OrderedCollections: OrderedDict

begin
  busses_dict = OrderedDict(
    "bus1"=>SlackAlgebraic(U=1.),
    "bus2"=>SwingEqLVS(H=1, P=1, D=1, Ω=2π*50,Γ=0.1,V=1.))

  busses =[SlackAlgebraic(U=1.),SwingEqLVS(H=1, P=1., D=1, Ω=2π*50,Γ=0.1,V=1.)]
  lines= [PiModelLine(;from=1, to=2, y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.)]
  #,"line2"=>  PiModelLine(;from=2, to=3, y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.))

  lines_dict = OrderedDict("line1"=> PiModelLine(;from="bus1", to="bus2", y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.))
  #,"line2"=>  PiModelLine(;from=2, to=3, y = 1/(0.1152 + im*0.0458), y_shunt_km = 0.,  y_shunt_mk = 0.))
end#

begin
  pg = PowerGrid(busses, lines);
  #pg_dict = PowerGrid(collect(values(busses_dict)), collect(values(lines_dict)));
  pg_dict = PowerGrid(busses_dict, lines_dict);

  #nsc = NodeShortCircuit(node = 2, R = 0.3, sc_timespan=(1.,2.));

  # a short circuit at the slack should have no effect
  nsc = NodeShortCircuit(;
      node_number = 2,
      Y = complex(160., 0.),
      tspan_fault = (0.5, 0.65),
  )

  ic_guess = ones(systemsize(pg_dict))# randn(systemsize(pg))#+ +

  pd = PowerPerturbation(node_number="bus2",fraction=0.1,tspan_fault=[0.1,0.2])
end 

begin # simulation
  #op = find_operationpoint(pg_dict)
  op_dict = find_operationpoint(pg_dict)
  #sol_2 = simulate2(nsc,pg,op_point,(0.,10.))
  sol_nc = simulate(pd,pg_dict,op_dict,(0.,5.))
end

begin # plotting
  include("plotting.jl");
  p2=plot(sol_nc,"bus1",:p)
  display(p2)
  #p=plot_res(sol_nc,pg_dict,"bus2")
  #display(p)
  #plot_res(sol_2,pg,2)
end
p2=plot(sol_nc,"bus1",:p)
display(p2)