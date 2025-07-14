using BlockSystems
using MarinePowerDynamics.IOComponents

export PQ_Load

"""
    PQ_Load(;P, Q, name=gensym(:PQ_load))

Implementation of `PQAlgebraic` as an IONode.

             +-----+
    i_r(t) --|  P  |-- u_r(t)
    i_i(t) --|  Q  |-- u_i(t)
             +-----+
"""
function PQ_Load(;P, Q, name=gensym(:PQ_load))
    para = Dict(:P => P,    # active power set point
                :Q => Q)    # reactive power set point

    pq_constraint = PowerConstraint(;name=name)
    IONode(pq_constraint, para)
end
