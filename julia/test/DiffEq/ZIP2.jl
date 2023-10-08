using DifferentialEquations
using Interpolations
using Plots
using CalculusWithJulia
using Optimization, Optim, OptimizationOptimJL, OptimizationBBO
using Sundials
using .Data
using CSV
using DataFrames
using Dictionaries: Dict
using LineSearches
t = Data.t
V_input = Data.v_fit
V_mean = Data.v_mean
P_meas = Data.p_mean
Q_meas = Data.q_mean

pP = [0.5,0.5]
pQ = [0.5,0.5]
function ZIP_P(p,V)
    P0_1, P0_2 = p
    #P0_1 + P0_2 = 1
    V0_1 = 1.0
    V0_2 = 1.0
    Zp1 = 1.5
    Zp2 = 0.4 
    Ip1 = -2.31
    Ip2 = -0.41
    Pp1 = 1.81
    Pp2 = 1.01
    P =  P0_1*(Zp1*(V/V0_1).^2 .+Ip1*(V/V0_1) .+Pp1)+P0_2*(Zp2*(V/V0_2).^2 .+Ip2*(V/V0_2) .+Pp2)
    return P
end
function ZIP_Q(p,V)
    Q0_1,Q0_2 = p
    #Q0_1 + Q0_2 = 1 
    V0_1 = 1.0
    V0_2 = 1.0
    Zq1 = 1.5
    Zq2 = 0.4 
    Iq1 = -2.31
    Iq2 = -0.41
    Qq1 = 1.81
    Qq2 = 1.01
    Q =  Q0_1*(Zq1*(V/V0_1).^2 .+Iq1*(V/V0_1) .+Qq1)+Q0_2*(Zq2*(V/V0_2).^2 .+Iq2*(V/V0_2) .+Qq2)
    return Q
end

function objective_function_P(p)
    sim_data = ZIP_P(p,V_mean)
    residuals = sim_data - P_meas
    sum_res = sum(residuals.^2)
    return sum_res
end

function objective_function_Q(p)
    sim_data = ZIP_Q(p,V_mean)
    residuals = sim_data - Q_meas
    sum_res = sum(residuals.^2)
    return sum_res
end

lower = [0., 0.]
upper = [1., 1.]
initial_x = [0.5, 0.5]
# requires using LineSearches
inner_optimizer = GradientDescent(linesearch=LineSearches.BackTracking(order=3))
result_P = optimize(objective_function_P, lower, upper, initial_x, Fminbox(inner_optimizer))


#result_P = optimize(objective_function_P,pP,LBFGS())
p_fit_P = Optim.minimizer(result_P) 

result_Q = optimize(objective_function_Q,pQ,LBFGS())
p_fit_Q = Optim.minimizer(result_Q)


fitting = Dict(
    "p_fit_P" => p_fit_P,
    "p_fit_Q" => p_fit_Q)

file = DataFrame(fitting)

path = raw"C:\Users\carlr\Desktop\UniStuttgart\Bachelor\Simulationsergebnisse\ersteDaten.csv"
sim_P = ZIP_P(p_fit_P,V_mean)
visual = plot(time,P_meas,label="measurement");plot!(time,sim_P,label="sim_P")
display(visual)
#CSV.write(fitting,path)

#=
tspan = (0.,10.)
differential_vars = false
f = DAEFunction(ZIP;syms=[:P])
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