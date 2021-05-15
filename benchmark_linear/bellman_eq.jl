"""
bellman_eq(a1::Float64, a0::Float64,
           grid::Array{Float64,1}, vfcn::Array{Float64,1},
           β::Float64, γ::Float64, r::Float64, w::Float64)

Purpose:
Bellman equation.
Given next period's value function, return the current value.
Use linear interpolation.

Input:
a1 → next period's asset (OPTIMIZE),
a0 → current asset,
grid → asset grid,
vfcn → next period's value function,
β → discount factor,
γ → inverse of ies,
r → interest rate,
w → wage.

Output:
v0 → current period's value.
"""
function bellman_eq(a1::Float64, a0::Float64,
                    grid::Vector{Float64}, vfcn::Vector{Float64},
                    β::Float64, γ::Float64, r::Float64, w::Float64)
    # consumption
    cons = (1 + r)*a0 + w - a1
    # CRRA utility function
    cons > 0 ? util = crra(cons, γ) : util = -1000000.0
    # use linear interpolation to compute next period's value
    # need to install and load "Interpolations"
    li = LinearInterpolation(grid, vfcn; extrapolation_bc=Line())
    v1 = li(a1)
    # value function
    v0 = util + β*v1
    # use Minimization
    v0 = -1.0*v0
    return v0
end
