using .PowerDynamics
using OrderedCollections: OrderedDict
using Plots
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\components.jl")


buses = OrderedDict(
      "bus1" => SlackAlgebraic(U=1),
      "bus2" => SwingEqLVS(H=5., P=1., D=2., Ω=50., Γ=20., V=1.0), 
      );  

buses2 = OrderedDict(
        "bus1" => SlackAlgebraic(U=1),
        "bus2" => FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=1., H=5., E_f=1),
        );  

buses1 = OrderedDict(
      "bus1" => SwingEqLVS(H=100, P=-1., D=20., Ω=50., Γ=20., V=1.0), 
      "bus2" => FourthOrderEq(T_d_dash=5.4, D=2, X_d=0.8979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=1., X_q_dash=0.646, P=1., H=5., E_f=1.),
      );  

Z₀ = 16 # per unit Impedanz   [Ω] 
length = 10.
#Leitungsparameter 
R = 0.25; #Ω/km
X = 0.12; 

branches=OrderedDict(
    "branch1" => StaticLine(from= "bus1", to = "bus2", Y = Z₀/(length*(R+X*im))),
     );

pg = PowerGrid(buses, branches)
pg1 = PowerGrid(buses1, branches)
pg2 = PowerGrid(buses2, branches)
  
operationpoint = find_operationpoint(pg,solve_powerflow=true,sol_method =:dynamic)

timespan= (0.0,20.)
fault = PowerPerturbation(node= "bus2",fault_power=0.5,tspan_fault = (5.,10),var=:P)
fault2 = NodeParameterChange(node = "bus1", value = 0.9, tspan_fault = (10.,20.),var=:V)
#state = State(pg1, operationpoint[:])  
#solution = simulate(fault, state , timespan)

sol= simulate(fault, pg2, (insert!(operationpoint[:],5,0.0)), timespan)
sol1= simulate(fault2, pg1, insert!(insert!(operationpoint[:],5,0.0),3,0.0), timespan)

v_sim=sol(sol.dqsol.t,:,:v)
v_sim1=sol1(sol1.dqsol.t,:,:v)
p_sim=sol(sol.dqsol.t,:,:p)
p_sim1=sol1(sol1.dqsol.t,:,:p)
#ω_sim=sol(sol.dqsol.t,:,:ω)
ω_sim1=sol1(sol1.dqsol.t,:,:ω)  

plot_v= plot(sol.dqsol.t,v_sim',label=["Slack" "Fourth"])
plot_v1=plot(sol1.dqsol.t,v_sim1',label=["Swing" "Fourth"],title="V")
plot_p = plot(sol.dqsol.t,p_sim',label=["Slack" "Fourth"])
plot_p1=plot(sol1.dqsol.t,p_sim1',label=["Swing" "Fourth"],title="P")
plot_ω1=plot(sol1.dqsol.t,ω_sim1'/(2*pi),label=["Swing" "Fourth"],title="f in Hz")
plot( plot_v,plot_v1,plot_p,plot_p1,plot_ω1;
        layout=(3,2),
        size = (1000, 500),
        lw=1,
        plot_title = "PvFourthOrder",
        xlabel="t[s]")