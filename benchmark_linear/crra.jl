"""
crra(cons::Float64, γ::Float64)

Purpose:
CRRA type utility function.

Input:
cons → consumption,
γ → relative risk aversion (inverse of the intertemporal elasticity of substitution).

Output:
util → current utility level.
"""
function crra(cons::Float64, γ::Float64)
    if γ == 1.0
        util = log(cons)
    else
        util = cons^(1-γ) / (1-γ)
    end
    return util
end
