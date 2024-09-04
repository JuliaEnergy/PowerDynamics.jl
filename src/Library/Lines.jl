@mtkmodel LineEnd begin
    @components begin
        terminal = Terminal()
    end
    @variables begin
        u_r(t), [description="line end d-voltage", input=true]
        u_i(t), [description="line end q-voltage", input=true]
        i_r(t), [description="line end d-current", output=true]
        i_i(t), [description="line end d-current", output=true]
    end
    @equations begin
        u_r ~  terminal.u_r
        u_i ~  terminal.u_i
        i_r ~ -terminal.i_r
        i_i ~ -terminal.i_i
    end
end

function Line(_branches...; name=:line)
    branches = Tuple[]
    for branch in _branches
        if branch isa Tuple
            push!(branches, branch)
        else
            push!(branches, (branch, branch.terminal1, branch.terminal2))
        end
    end
    systems = @named begin
        src = LineEnd()
        dst = LineEnd()
    end

    eqs = [[connect(src.terminal, inj[2]) for inj in branches]...,
           [connect(dst.terminal, inj[3]) for inj in branches]...]

    ODESystem(eqs, t; systems=[ systems..., first.(branches)...], name)
end

# @mtkmodel LineEnds begin
#     @components begin
#         src = Terminal()
#         dst = Terminal()
#     end
#     @variables begin
#         src_u_r(t), [description="src end d-voltage", input=true]
#         src_u_i(t), [description="src end q-voltage", input=true]
#         src_i_r(t), [description="src end d-current", output=true]
#         src_i_i(t), [description="src end d-current", output=true]
#         dst_u_r(t), [description="dst end d-voltage", input=true]
#         dst_u_i(t), [description="dst end q-voltage", input=true]
#         dst_i_r(t), [description="dst end d-current", output=true]
#         dst_i_i(t), [description="dst end d-current", output=true]
#     end
#     @equations begin
#         src_u_r ~ src.u_r
#         src_u_i ~ src.u_i
#         src_i_r ~ src.i_r
#         src_i_i ~ src.i_i
#         dst_u_r ~ dst.u_r
#         dst_u_i ~ dst.u_i
#         dst_i_r ~ dst.i_r
#         dst_i_i ~ dst.i_i
#     end
# end
