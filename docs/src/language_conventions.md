PowerDynBasePowerDynBase# Language & Conventions

Generally, variables are miniscule (e.g. `u`, `i`, `ω`) and parameters are capitalized (e.g. `H`, `D`, `P`, `Ω`). As it is common to use greek letters for modeling equations and Julia supports
Unicode, greek letters are used within the Code, e.g. `Ω` and `ω` in [`PowerDynBase.SwingEq`](@ref). If you don't want to use greek keyboard (which I am currently switching to) you can simply type the latex representating `\Omega` and Atom can complete it with `Ω` using Tab.

## List of symbols and corresponding names

| Symbol (Code) | Symbol (Math) | Name within `PowerDynamics.jl` | Common alternative names |
|:-------------:|:-------------:|:---------------------:|:------------------------:|
|   |   | node  | bus, vertex  |
|   |   | grid  | network, power grid, power network  |
|  | $y_{ab} = y_{ba}$  | admittance between nodes $a$ and $b$  |   |
| `LY` | $Y^L$  | admittance laplacian  | (nodal) admittance matrix  |
| `t`  | $t$  | time  |  |
| `im`  | $j$  | imaginary element  | $\sqrt{-1}$  |
| `u = v \cdot exp(im*φ)`  | $u = v \cdot e^{jφ}$  | complex voltage  |   |
| `v`   | $v$  | voltage magnitude  | absolute voltage  |
| `φ`  | $\phi$  | voltage angle  |   |
| `i_c = i \cdot exp(im*δ)`  | $i_c = i \cdot e^{j\delta}$  | nodal complex current  |   |
| `i`  | $i$  | magnitude of the current  |   |
| `δ`  | $\delta$  |  angle of the current |   |
|  `s = p + im*q` | $s = p + jq$  | complex power   |   |
| `p`   |$p$   | real power   | active power  |
|`q`   |  $q$ | imaginary power  | reactive power  |

## List of modeling conventions

- Counting of nodes starts at 1.
- Ranges of nodes are mathematical, i.e. they include the first and the last element. For example $\sum_k=3^6$ sums over $3$, $4$, $5$, and $6$.
- For now, no selfadmittance is allowed, i.e. $y_{aa} = 0$ for all nodes $a$.
- The admittance laplacian uses the following definition ([convention from wikipedia](https://en.wikipedia.org/wiki/Nodal_admittance_matrix#Construction))
```math
Y^L_{ab} = \begin{cases}
  \sum_{c} y_{ac} & \text{if } a=b, \\
  -y_{ab} & \text{otherwise.}
\end{cases}
```
- The nodal complex current is calculated as
```math
{i_c}_a = \sum_{b} LY_{ab} u_b .
```
- The complex power is calculated as (with ``^*`` as complex comjucation)
```math
s_a = u_a \cdot {i_c}_a^*.
```
