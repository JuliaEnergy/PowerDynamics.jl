# (C) 2018 Potsdam Institute for Climate Impact Research, authors and contributors (see AUTHORS file)
# Licensed under GNU GPL v3 (see LICENSE file)

"Interpret (part of) an array of real values as an array with complex values."
function complexview(vec::AbstractArray{T}, i0, n) where T
    vec_view = view(vec, i0:i0 + 2*n -1)
    reinterpret(Complex{T}, vec_view)
end
