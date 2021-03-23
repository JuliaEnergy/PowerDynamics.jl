using BlockSystems
using ModelingToolkit
using PowerDynamics.IOComponents

export PQ_Load

"""
    PQ_Load(;P,Q)

Implementation of `PQAlgebraic` as an IONode.

         +-----+
i_r(t) --|  P  |-- u_r(t)
i_i(t) --|  Q  |-- u_i(t)
         +-----+

"""

function PQ_Load(;P,Q)

    para = Dict(:P => P,    # active power set point
                :Q => Q)    # reactive power set point

    pq_constraint = PowerConstraint(name=:pq_constraint, P=:P, Q=:Q)

    pql = IOSystem([],
                   [pq_constraint],
                   name=:PQLoad,
                   outputs=[pq_constraint.u_i, pq_constraint.u_r])

    connected = connect_system(pql)
    IONode(connected, para)
end