using DifferentialEquations
using Interpolations
using Plots
using CalculusWithJulia
using Optimization, Optim, OptimizationOptimJL, OptimizationBBO
using Sundials
#https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1429557

#https://a-lab.ee/sites/default/files/PS_Lecture_2.pdf

#= 
Power Load :


Infinite bus 

Zwei bus, eine Leitung, Admittanz ist erster Parameter,
Spannung Bus1 ist Eingang, Zustände sind Strom und Spannung Bus 2
=#
function network(out,du,u,p,t)
    #PQ Constraint 
    out[1] = 1 
end

function network1(du,u,p,t)
    du[1] = 2*pi*50*u[2]+2*(p[1]/p[2])*((u[1])/(u[1]^2+u[2]^2))-V(t)/p[2]
    du[2] = -2*pi*50*u[1] + 2*(p[1]/p[2])*((u[1])/(u[1]^2+u[2]^2))
end

function ZIP(out,du,u,p,t)
    out[1] = (u[1] -p[1]*(p[2] + p[3]*(V(t)/1.)+p[4]*(V(t)/1.)^2))-du[1]
    out[2] = (u[2] -p[5]*(p[6] + p[7]*(V(t)/1.)+p[8]*(V(t)/1.)^2))-du[2]
end

function PTD1(du,u,p,t)
    V(t) = 1. +(0.1)*sign(t-10.)-(0.)*sign(t-10.)
    
    du[1] = (1/p[1])*(p[2]*V(t)-u[1])
    du[2] = (1/p[3])*(p[4]*V(t)-u[2])
end

function RLC(du,u,p,t)
    R, L , C= p
    V(t) = 1. +(0.1)*sign(t-10.)-(0.)*sign(t-10.)
    du[1] = u[2]/C
    du[2] = -u[1]/L-(R/L)*u[2]+V(t)/L
end
# p = [admittance, P,Q]
#p = [0.5,0.8,0.1,0.1,0.5,0.5,0.3,0.2]
p = [0.9,0.9,0.5]
u₀ = [0.9, 0.]
du₀ = [0., 0.]
tspan = (0.0, 20.)



#Vabs = VectorOfArray([(0.7 + .01randn()) for i in 1:22])

        
#differential_vars = [false, false]
#f = DAEFunction(ZIP;syms=[:p,:q])
#prob = DAEProblem(f, du₀, u₀, tspan, differential_vars = differential_vars,p)
f = ODEFunction(RLC;syms=[:v_c,:i])
prob = ODEProblem(f,u₀,tspan,p)
sol = solve(prob, Tsit5())
#sol = solve(prob, IDA())
#=
iabs = Float64[]
for i=1 : size(sol.t)[1]
    push!(iabs,sqrt(sol[:i_d][i]^2+sol[:i_q][i]^2))
end
plot(sol.t,iabs)
=#