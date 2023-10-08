using DifferentialEquations
#] add DiffEqProblemLibrary
using DiffEqProblemLibrary.ODEProblemLibrary
# load problems

prob = ODEProblemLibrary.prob_ode_linear
sol = solve(prob)
