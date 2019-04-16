using PowerDynBase
using NetworkDynamics
using Test
using MacroTools
#include("../src/MassMatrix.jl")
@testset "Macro 2" begin
    exp = quote
        m_u=false
        m_int=[false, true]
    end
    #,:(m_int=no_internal_masses)
    @capture(exp, fields__)
    #println(eval(fields[1].args[1]))
    m_u = eval(fields[1].args[1])
    m_int = eval(fields[2].args[1])
    if length(fields) == 0
        mass_list = nothing
    else
        mass_list = [m_u]
        append!(mass_list,m_int)
    end
    println(mass_list)
end

mass_list =  nothing, [1, 0.] [false, true]
