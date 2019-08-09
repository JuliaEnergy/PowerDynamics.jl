Pkg.instantiate()
cd(@__DIR__)
using PowerDynamics

begin
    node_list = []
    append!(node_list, [CSIMinimal(
        I_r=1.,Y_n=1.
                )])
end
