using LinearAlgebra: I


begin
    @doc """
    ```Julia
    RLLine(; from, to, R, L, ω0)
    ```
    dynamic line with series resistance R and series inductance L.

    When setting initial conditions for this line type be aware that every line carries four variables. The first two for the current from source to destination and the second from destination to source. Since i_src_dst = - i_dst_src, the inital condtions need to be chosen as negatives of each other as well, i.e. ic[1:2] == -ic[3:4].

    # Keyword Arguments
    - `from`: start node of the line
    - `to`: end node of the line
    - `R`: series resistance R
    - `L`: series inductance L
    - `ω0`: rated frequency in [rad/s] of the power grid, often ``2π50``Hz

    # Mathematical Representation
    Using `RLLine` for the line `from`--`to` applies Eqn. (2) from
    __Brouillon, J. S., Colombino, M., Groß, D., & Dörfler, F. (2018). The effect of transmission-line dynamics on a globally synchronizing controller for power inverters. In 2018 European Control Conference (ECC) (pp. 2242-2247). IEEE.__
    """
    struct RLLine <: AbstractLine
            from
            to
            R
            L
            ω0
            RLLine(; from, to, R, L, ω0) = new(from, to, R, L, ω0)
        end
    function construct_edge(par::RLLine)
        from = par.from
        to = par.to
        R = par.R
        L = par.L
        ω0 = par.ω0
        @assert R > 0 "resistance (R) should be positive"
        @assert L > 0 "inductance (L) should be positive"
        @assert ω0 > 0 "rated frequency (ω0) should be positive"
        Z = [R -ω0*L; ω0*L R]
        function rhs!(de, e, v_s, v_d, p, t)
            i_left_right = [e[1]; e[2]]
            # i_right_left = [e[3]; e[4]]
            v_left = [v_s[1]; v_s[2]]
            v_right = [v_d[1]; v_d[2]]
            di_left_right = (-(Z * i_left_right) .+ v_right .- v_left) ./ L
            de .= [di_left_right; -di_left_right]
        end
        return ODEEdge(f! = rhs!, dim=4, mass_matrix=I, sym=Symbol[:id, :iq, :id_r, :iq_r], coupling = :fiducial)
    end
    symbolsof(::RLLine) = begin
            [:id, :iq, :id_r, :iq_r]
    end
    dimension(::RLLine) = begin
            4
    end
end

export RLLine
