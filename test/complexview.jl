
include("testing_base.jl")

function cmp_complex_real(c_arr, r_arr)
        cmp_arr = Array{Bool}(undef, length(c_arr))
        for i in 1:length(c_arr)
                cmp_arr[i] = (c_arr[i] == r_arr[2*i-1] + im*r_arr[2*i])
        end
        all(cmp_arr)
end
num_real = 21 # should be an odd number
@assert isodd(num_real) "test for odd numbers, too"
num_complex = 3
@assert 2*num_complex <= num_real "need to get less complex number than half the real numbers"
real_array = rand(num_real)
complex_array = PowerDynBase.complexview(real_array, 1, 3)
@test length(complex_array) == 3
@test cmp_complex_real(complex_array, real_array)
complex_array[:] = rand(3) .+ im.*rand(3)
@test cmp_complex_real(complex_array, real_array)
real_array[:] = rand(num_real)
@test cmp_complex_real(complex_array, real_array)
