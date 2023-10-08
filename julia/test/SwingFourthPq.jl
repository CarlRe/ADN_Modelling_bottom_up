using .PowerDynamics
using OrderedCollections: OrderedDict
using Plots
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\components.jl")
  
buses = OrderedDict(
      "bus1" => SwingEqLVS(H=5., P=0., D=1., Ω=50., Γ=20., V=1.0),
      "bus2" => SwingEqLVS(H=5., P=1., D=1., Ω=50., Γ=20., V=1.0), 
      "bus3" => PQAlgebraic(P=-1., Q= 0.),
      );  

buses1 = OrderedDict(
      "bus1" => SlackAlgebraic(U=1),
      "bus2" => FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=1., H=5.148, E_f=1),
      "bus3" => PQAlgebraic(P=-1., Q= 0.),
      );  

buses2 = OrderedDict(
        "bus1" => SwingEqLVS(H=5., P=0., D=1., Ω=50., Γ=20., V=1.0),
        "bus2" => FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=1., H=5.148, E_f=1),
        "bus3" => PQAlgebraic(P=-1., Q= 0.),
        );  

Z₀ = 16 # per unit Impedanz   [Ω] 
length = 1.
#Leitungsparameter 
R = 0.25; #Ω/km
X = 0.12;  
branches=OrderedDict(
    "branch1" => StaticLine(from= "bus1", to = "bus2", Y = Z₀/(length*(R+X*im))),
    "branch2" => StaticLine(from= "bus1", to = "bus3", Y = Z₀/(length*(R+X*im)))
     );

pg1 = PowerGrid(buses, branches)
pg1 = PowerGrid(buses1, branches)
pg2 = PowerGrid(buses2, branches) 
#=   
    data, result = power_flow(pg)
        v = [result["solution"]["bus"][string(k)]["vm"] for k in 1:3]
        va = [result["solution"]["bus"][string(k)]["va"] for k in 1:3]
=#
        
operationpoint = find_operationpoint(pg,solve_powerflow=true,sol_method =:dynamic)

#state = State(pg1, operationpoint[:])  

timespan= (0.0,25.)
fault = PowerPerturbation(node= "bus3",fault_power=-0.,tspan_fault = (10.,15.),var=:P)
    
sol= simulate(fault, pg1, insert!(operationpoint[:],5 ,0.0), timespan)
sol1= simulate(fault, pg2, insert!(insert!(operationpoint[:],5 ,0.0),3 ,0.0 ), timespan)

v_sim=sol(sol.dqsol.t,:,:v)
v_sim1=sol1(sol1.dqsol.t,:,:v)
p_sim=sol(sol.dqsol.t,:,:p)
p_sim1=sol1(sol1.dqsol.t,:,:p)

plot_v= plot(sol.dqsol.t,v_sim',label=["Slack" "Swing" "PQ"])
plot_v1=plot(sol1.dqsol.t,v_sim1',label=["Swing" "Fourth" "PQ"],title="V")
plot_p = plot(sol.dqsol.t,p_sim',label=["Slack" "Swing" "PQ"])
plot_p1=plot(sol1.dqsol.t,p_sim1',label=["Swing" "Fourth" "PQ"],title="P")

plot( plot_v,plot_v1,plot_p,plot_p1;
        layout=(2,2),
        size = (1000, 500),
        lw=3,
        plot_title = "SwingFourthPq",
        xlabel="t[s]")