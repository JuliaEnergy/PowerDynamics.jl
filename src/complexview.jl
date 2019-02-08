# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

"""
    complexview(vec::AbstractArray, i0, n)

Interpret (part of) an array of real values as an array with complex values.
`i0` is the index where to start.
`n` is the number of complex values that should be extracted.
"""
function complexview(vec::AbstractArray{T}, i0, n) where T
    vec_view = view(vec, i0:i0 + 2*n -1)
    reinterpret(Complex{T}, vec_view)
end
