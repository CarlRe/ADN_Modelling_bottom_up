using DifferentialEquations
using Interpolations
using Plots
using CalculusWithJulia

using Sundials
using Optim

#using Optimization, OptimizationOptimJL, OptimizationBBO
#=
function dynamic_load(du,u,p,t)
    V(t) = 0.9+-(0.01)*sign(t-5.)
    du[1] = (1/p[1])*(-u[1]+P0*((V(t)/V0))^p[2]- P0*((V(t)/V0)^p[3]))
    du[2] = (1/p[4])*(-u[2]+Q0*((V(t)/V0))^p[5]- Q0*((V(t)/V0)^p[6]))
    out[3]= (u[1]+((V(t)/V0)^(p[3])))-p[7]
    out[4]= (u[2]+((V(t)/V0)^(p[6])))-p[8] 
end
=#
time = Data.t
V_input = Data.v_fit
V_mean = Data.v_mean
P_meas = Data.p_mean
Q_meas = Data.q_mean
step_size = 4/40960
p2 = [1.,1.,1.]
p1=[1.,0.57,1.,2.07]
p = [1.247,0.068,1.,2.07,0.17]
u0 = [0.1]
tspan = (0.,4.)
function dynamic_load_new(p,V)
    V0,P0,a,b,c = p
    function ODE(du,u,p,t)
    V(t) = 1.247-(0.012)*sign(t-0.125)+(0.012)*sign(t-0.25)
    du[1] = (1/a)*(-u[1]+ P0*((V(t)/V0))^b - P0*((V(t)/V0)^c))
    end
    f = ODEFunction(dynamic_state;syms=[:xp])
    prob = ODEProblem(f,u0,tspan,p)
    sol = solve(prob, Tsit5(),saveat=step_size)
    xp = sol[:xp][1:end-1]
    P_sim = xp + (V./V0).^c
    return P_sim
end
function simple_expo(p,V)
    P0,V0,a = p
    return P0*(V/V0).^a

end
function dynamic_load_WithOutC(p,V)
    V0,P0,a,b = p
    function ODE(du,u,p,t)
    U = 1.247-(0.012)*sign(t-0.125)+(0.012)*sign(t-0.25)
    du[1] = (1/a)*(-u[1]+ P0*((U/V0))^b - P0*((U/V0)^1))
    end
    f = ODEFunction(ODE;syms=[:xp])
    prob = ODEProblem(f,u0,tspan,p)
    sol = solve(prob, Tsit5(),saveat=step_size)
    xp = sol[:xp][1:end-1]
    P_sim = xp + (V/V0).^1.
    return P_sim
end

function objective(p)
    data_sim = simple_expo(p,V_mean)
    residuals = data_sim - P_meas 
    sumres = sum(residuals.^2)
    return sumres
end
lower = [0.0, -1.0, 0.0, 0.0]
upper = [1.0, 1.0, 1.0, 1.0]

result = optimize(objective,p2,LBFGS())
p_fit6 = Optim.minimizer(result)



function dynamic_state(du,u,p,t)
    V0,P0,a,b,c = p
    V(t) = 0.9+-(0.01)*sign(t-5.)
    V = V(t)
    
end

#=
f = ODEFunction(dynamic_state;syms=[:xp])
prob = ODEProblem(f,u0,tspan,p)
sol = solve(prob, Tsit5(),saveat=step_size)
xp = sol[:xp][1:end-1]

x_p = Float64[]
for i=1:40960
    push!(x_p,sol[:xp][i])
end
=#
function dynamic_load(p,V)
    #=
    P_sim = Float64[]
    for i=1:40960
        push!(P_sim,P0*(V[i]/V0)^p[3]+sol[:xp][i])
    end
    =#
    #V0,P0,a,b,c = p
    f = ODEFunction(dynamic_state;syms=[:xp])
    prob = ODEProblem(f,u0,tspan,p)
    sol = solve(prob, Tsit5(),saveat=step_size)
    xp = sol[:xp][1:end-1]
    P_sim = xp + (V/p[0]).^p[5]
    return P_sim
end
#=
function dynamic_load(p,V)
    prob = ODEProblem(f,u₀,tspan,p)
    sol = solve(prob, Tsit5(),saveat=step_size)
    V0,P0, a,b,c = p
    xp = sol[:xp][1:40960]
    P_sim =P0.*(V/V0).^c .+ xp
    return P_sim
end
=#


#=
#=
p = [0.57,1.06,2.07,0.17,1.17,2.07,1.,0.1]
u0 = [0.,0.,0.5,0.1]
du0 = [0.,0.,0.,0.]
=#
tspan = (0.,10.)
differential_vars = [true, false]
f = DAEFunction(dynamic_load;syms=[:xp,:P])
prob = DAEProblem(f, du0, u0, tspan, differential_vars = differential_vars,p)
sol = solve(prob, IDA())

#=
#Vabs = VectorOfArray([(0.7 + .01randn()) for i in 1:22])
V(t) = 0.9+-(0.01)*sign(t-5.)

#Vabs_fun = time_interpolation(Vabs,dt)
prob_init = ODEProblem((du,u,p,t) -> dynamic_load(du,u,p,t),u0,tspan,p)

#sol_init = solve(prob_init,Rosenbrock23())

function find_steady_state()
for i=1:size(sol_init.t)[1]-1
    if abs((sol_init[1,i+1]-sol_init[1,i])/(sol_init.t[i+1]-sol_init.t[i]))<= 1e-6
        display(sol_init[i])
        display(i)
        return u_steady = sol_init[i]
    end
end
end

u_steady = find_steady_state()
V(t) = (0.9+0.00001randn())-(0.1)*sign(t-2)
prob = ODEProblem((du,u,p,t) -> dynamic_load(du,u,p,t),u_steady,tspan,p)

sol = solve(prob,Rosenbrock23())
#plot(sol)

using RecursiveArrayTools # for VectorOfArray
randomized = VectorOfArray([(sol(sol.t[i]) + .01randn(4)) for i in 1:length(sol.t)])
data = convert(Array,randomized)
cost_function = build_loss_objective(prob, Tsit5(), L2Loss(sol.t,data),
                                     Optimization.AutoForwardDiff(),
                                     maxiters=10000,verbose=false)
#x_p, x_q
power_active = Float64[]
power_reactive = Float64[]
for i=1:size(sol.u)[1] 
    push!(power_active,sol[1,i]+P0*((V(sol.t[i])/V0))^p[3])
    push!(power_reactive,sol[2,i]+Q0*((V(sol.t[i])/V0))^p[3])
end

plot(sol.t[:],power_active[:],label= "P")
plot!(sol.t[:],power_reactive[:], label = "Q")
plot!(sol.t[:],sol[3,:],label= "P aus State")
plot!(sol.t[:],sol[4,:], label = "Q aus State")

=#
=#