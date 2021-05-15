"""
bellman_eq_sp(a1::Float64, a0::Float64,
           spl::Spline1D,
           β::Float64, γ::Float64, r::Float64, w::Float64)

Purpose:
Bellman equation.
Given next period's value function, return the current value.
Use "spline interpolation".

Input:
a1 → next period's asset (OPTIMIZE),
a0 → current asset,
vfcn → next period's value function,
β → discount factor,
γ → inverse of ies,
r → interest rate,
w → wage.

Output:
v0 → current period's value.
"""
function bellman_eq_sp(a1::Float64, a0::Float64,
                       spl::Spline1D,
                       β::Float64, γ::Float64, r::Float64, w::Float64)
    # consumption
    cons = (1 + r)*a0 + w - a1
    # CRRA utility function
    cons > 0 ? util = crra(cons, γ) : util = -1000000.0
    # use cubic spline interpolation to compute next period's value
    # need to install and load "Dierckx"
    v1 = spl(a1)
    # value function
    v0 = util + β*v1
    # use Minimization
    v0 = -1.0*v0
    return v0
end
