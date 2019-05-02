using NetworkDynamics
using Parameters
import Base.convert

@with_kw struct StaticLine2!
    Y
end

function(sl::StaticLine2!)(e, v_s, v_d, p, t)
    source_voltage = v_s[1] + v_s[2] * im
    destination_voltage = v_d[1] + v_d[2] * im
    complex_current = sl.Y * (destination_voltage - source_voltage)
    push!(e, real(complex_current))
    push!(e, imag(complex_current))
    push!(e, real(complex_current))
    push!(e, imag(complex_current))
    nothing
end

convert(::Type{StaticEdge}, sl::StaticLine2!) = StaticEdge(f! = sl, dim = 4)
