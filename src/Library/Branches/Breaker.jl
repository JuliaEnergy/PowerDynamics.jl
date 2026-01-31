@mtkmodel Breaker begin
    @parameters begin
        closed=1, [description="Breaker closed (1) or open (0)"]
    end
    @components begin
        src = Terminal()
        dst = Terminal()
    end
    @variables begin
        i_r(t), [guess=0, description="Current real part through breaker"]
        i_i(t), [guess=0, description="Current imaginary part through breaker"]
    end
    @equations begin
        ## When closed: enforce u_dst = u_src (with current as implicit output)
        ## When open: enforce i = 0
        0 ~ ifelse(closed == 1, dst.u_r - src.u_r + NetworkDynamics.implicit_output(i_r), i_r)
        0 ~ ifelse(closed == 1, dst.u_i - src.u_i + NetworkDynamics.implicit_output(i_i), i_i)
        ## Current flow equations (standard for all edge components)
        dst.i_r ~ i_r
        dst.i_i ~ i_i
        src.i_r ~ -dst.i_r
        src.i_i ~ -dst.i_i
    end
end
