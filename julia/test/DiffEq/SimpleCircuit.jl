using DifferentialEquations
using Interpolations
using Plots
using CalculusWithJulia
using Optimization, Optim, OptimizationOptimJL, OptimizationBBO
using Sundials
using OrdinaryDiffEq

function RL(du,u,p,t)
    R, L = p
    V(t) = 1.24 -(0.01)*sign(t-5.)-(0.)*sign(t-7.)
    du[1] = (1/L)*(V(t)-R*u[1])
    #out[1] = (1/L)*(V(t)-R*u[1]) - du[1]
    #out[2] = R*u[1]-u[2]
end
# p = [admittance, P,Q]
#p = [0.5,0.8,0.1,0.1,0.5,0.5,0.3,0.2]
p = [0.5,0.,0.1,0.1]
u₀ = [1.,0.]
du₀ = [0.,0.]
tspan = (0.0, 10.)
function RLC(out,du,u,p,t)
    R, L, C , V_c = p
    V(t) = 1. +(0.5)*sign(t-3.)-(0.1)*sign(t-7.)
    
    du[1]  = u[2]
    out[1] = (1/L)*(V(t)-u[1]*R-V_c)-du[1]
    out[2] = (-1/L)*(R*u[2]+u[1]/C)-du[2]

    #du[1] = (1/L)*(-R*du[1]-u[1]/C)
    
end


#Vabs = VectorOfArray([(0.7 + .01randn()) for i in 1:22])

        
differential_vars = [ false,true]
f = DAEFunction(RL;syms=[:i,:di,:ddi])
prob = DAEProblem(f, du₀, u₀, tspan, differential_vars = differential_vars,p)
#f = ODEFunction(RL;syms=[:i,:di])
#prob = ODEProblem(f,u₀,tspan,p)
#sol = solve(prob, Tsit5())
sol = solve(prob, IDA())
#=
iabs = Float64[]
for i=1 : size(sol.t)[1]
    push!(iabs,sqrt(sol[:i_d][i]^2+sol[:i_q][i]^2))
end
plot(sol.t,iabs)
=#