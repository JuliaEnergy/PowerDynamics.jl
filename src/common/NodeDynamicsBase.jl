# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

#FIXME if removed we get MethodError in PowerGridSolution -> investigate this
import Base: convert, promote_rule

"A variable to be used when no internal masses are present for a node dynamics type."
const no_internal_masses = Vector{Bool}()

function MassMatrix(;m_u::Bool = false, m_int = no_internal_masses)
    mm = [m_u, m_u] # double mass for u, because it is a complex variable
    append!(mm, m_int)
    return mm .|> Int
end
