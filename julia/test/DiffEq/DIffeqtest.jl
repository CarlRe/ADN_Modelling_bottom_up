using DifferentialEquations
#] add DiffEqProblemLibrary
using DiffEqProblemLibrary.ODEProblemLibrary
# load problems
#ODEProblemLibrary.importodeproblems()
using Optimization
using DiffEqParamEstim
using OptimizationOptimJL
prob1 = ODEProblemLibrary.prob_ode_pollution
sol1 = solve(prob)

function g(du,u,p,t)
    du[1] = dx = p[1]*u[1] - u[1]*u[2]
    du[2] = dy = -3*u[2] + u[1]*u[2]
end
  
u0 = [1.0;1.0]
tspan = (0.0,10.0)
p = [1.5]
prob = ODEProblem(g,u0,tspan,p)

sol = solve(prob,Tsit5())
t = collect(range(0,stop=10,length=200))
using RecursiveArrayTools # for VectorOfArray
randomized = VectorOfArray([(sol(t[i]) + .01randn(2)) for i in 1:length(t)])
data = convert(Array,randomized)
newprob = remake(prob, p = [1.42])
newsol = solve(newprob,Tsit5())
plot(sol)
plot!(newsol)

cost_function = build_loss_objective(prob, Tsit5(), L2Loss(t,data),
                                     Optimization.AutoForwardDiff(),
                                     maxiters=10000,verbose=false)

vals = 0.0:0.1:10.0
plot(vals,[cost_function(i) for i in vals],yscale=:log10,
     xaxis = "Parameter", yaxis = "Cost", title = "1-Parameter Cost Function",
     lw = 3)

optprob = Optimization.OptimizationProblem(cost_function, [1.42])
optsol = solve(optprob, BFGS())

newprob = remake(prob, p = optsol.u)
newsol = solve(newprob,Tsit5())
plot(sol)
plot!(newsol)