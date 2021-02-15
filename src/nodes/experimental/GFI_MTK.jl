# First attempt at VSC in BlockSystems (nee IOSystems)

# ]add "https://github.com/hexaeder/IOSystems_prototype"
using BlockSystems
using ModelingToolkit


## Common time

@parameters t 
D = Differential(t)

## Low pass filter

@parameters τ input(t) 
@variables filtered(t)

lpf = IOBlock([D(filtered) ~ 1/τ * (- filtered + input)], [input], [filtered], name=:lp_filter)

## integrator

@parameters x(t)
@variables int(t)
    
integrator = IOBlock([D(int) ~ x], [x], [int], name=:integrator)

## Droop control

@parameters K u_ref x_ref x(t)
@variables u(t)

droop_control = IOBlock([
    u ~ - K * (x - x_ref) + u_ref # output is the droop voltage v
    ], [x], [u], name = :droop)

##

p_filter = IOBlock(lpf, name = :active_power_filter)
q_filter = IOBlock(lpf, name = :reactive_power_filter)
p_droop = IOBlock(droop_control, name = :active_power_droop)
q_droop = IOBlock(droop_control, name = :reactive_power_droop)
f_integrator = IOBlock(integrator, name = :frequency_integrator)

##
@variables ϕ(t) v(t)
@parameters p(t) q(t)

gfi = IOSystem([f_integrator.x => p_droop.u,
          p_droop.x => p_filter.filtered,
          q_droop.x => q_filter.filtered],
          [p_filter, q_filter, p_droop, q_droop, f_integrator],
          name = :GridForming,
          namespace_map = [p_filter.input => p, q_filter.input => q, f_integrator.int => ϕ, q_droop.u => v],
          outputs = [ϕ, v])

connected_gfi = connect_system(gfi)

##

gen = generate_io_function(connected_gfi, f_states=[v, ϕ], f_inputs=[p, q])

##


##
begin
  Base.@__doc__ struct GFI <: AbstractNode
          τ_P
          τ_Q
          P_ref
          V_refP
          K_P
          V_refQ 
          Q_ref 
          K_Q 
          Y_n
      end
  GFI(; τ_P,τ_Q,P_ref,V_refP,K_P,V_refQ,Q_ref,K_Q, Y_n = 0) = GFI(τ_P,τ_Q,P_ref,V_refP,K_P,V_refQ,Q_ref,K_Q,Y_n)
  function construct_vertex(par::GFI)
      τ_P = par.τ_P#par.active_power_filter₊τ
      τ_Q = par.τ_Q#reactive_power_filter₊τ
      P_ref = par.P_ref#active_power_droop₊x_ref
      V_refP = par.V_refP#active_power_droop₊u_ref
      K_P = par.K_P#active_power_droop₊K
      V_refQ = par.V_refQ#reactive_power_droop₊u_ref
      Q_ref = par.Q_ref#reactive_power_droop₊x_ref
      K_Q = par.K_Q#reactive_power_droop₊K
      Y_n = par.Y_n
      function rhs!(dx, x, e_s, e_d, p, t)
        # states sind in gen.states
        # TODO current sollte input sein, voltage sollte output sein
          i = total_current(e_s, e_d) + Y_n * (x[1] + x[2] * im)#Y_n * x[1]*exp(1im*x2)
          u = x[1] + x[2] * im #u = x[1]*exp(x[2] * im)
          active_power = real(u * conj(i))
          reactive_power = imag(u * conj(i))
          # TODO input statt reactive_power und active_power i_r and i_i
          odefun(dx, x, p, t) = gen.f_ip(dx, x, [active_power,reactive_power], p, t)
          try
              dx[1] = real(du)
              dx[2] = imag(du)
              return nothing
          catch e
              if typeof(e) === UndefVarError
                  throw(NodeDynamicsError("you need to provide $(e.var)"))
              else
                  throw(e)
              end
          end
      end
      ODEVertex(f! = rhs!, dim = length(gen.states), mass_matrix = gen.massm, sym = [:u_r, :u_i, :p_fil, :q_fil])
  end
  symbolsof(::GFI) = begin
          [:u_r, :u_i,:p_fil, :q_fil]
      end
  dimension(::GFI) = begin
          length(gen.states)
      end
end
##

export GFI