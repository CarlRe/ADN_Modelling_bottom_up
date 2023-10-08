using .PowerDynamics
using OrderedCollections: OrderedDict
using Plots
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\components.jl")


buses = OrderedDict(
      "bus1" => SlackAlgebraic(U=1),
      "bus2" => SwingEqLVS(H=1., P=1., D=1., Ω=50., Γ=20., V=1.0), 
      );  

buses2 = OrderedDict(
        "bus1" => SwingEqLVS(H=5., P=-1., D=1., Ω=50., Γ=20., V=1.0), 
        "bus2" => FourthOrderEq(T_d_dash=2.4, D=1., X_d=0.8979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=1., H=1., E_f=1),
        );  

buses1 = OrderedDict(
      "bus1" => SlackAlgebraic(U=1),
      "bus2" => FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=1., H=5.148, E_f=1.1),
      );  

Z₀ = 16 # per unit Impedanz   [Ω] 
length = 10.
#Leitungsparameter 
R = 0.25; #Ω/km
X = 0.5; 

branches=OrderedDict(
    "branch1" => StaticLine(from= "bus1", to = "bus2", Y = Z₀/(length*(R+X*im))),
     );

pg = PowerGrid(buses, branches)
pg1 = PowerGrid(buses1, branches)
pg2 = PowerGrid(buses2, branches)
  
operationpoint = find_operationpoint(pg1,solve_powerflow=true,sol_method =:dynamic)

timespan= (0.0,5.)
fault = PowerPerturbation(node= "bus2",fault_power=0.5,tspan_fault = (12.,13.),var=:P)
fault1 = ChangeInitialConditions(node="bus1", var=:v, f=Inc(0.2))
fault2 = NodeParameterChange(node = "bus1", value = 1, tspan_fault = (5.,20.),var=:U)
fault3 = NodeParameterChange(node = "bus1", value = 0.9, tspan_fault = (1.,5.),var=:V)
fault4 = NodeParameterChange(node = "bus2", value = 0.9, tspan_fault = (1.,20.),var=:E_f)
#sol= simulate(fault2, pg2, insert!(operationpoint[:],6,0.0), timespan)
#sol1= simulate(fault3, pg2, insert!(operationpoint[:],size(operationpoint[:])[1],-0.1), timespan)
sol1= simulate(fault3, pg2,insert!(operationpoint[:],3,0.0), timespan)
#v_sim=sol(sol.dqsol.t,:,:v)
v_sim1=sol1(sol1.dqsol.t,:,:v)
#p_sim=sol(sol.dqsol.t,:,:p)
p_sim1=sol1(sol1.dqsol.t,:,:p)
ω_sim1=sol1(sol1.dqsol.t,["bus1" ,"bus2"],:ω)   

#plot_v= plot(sol.dqsol.t,v_sim',label=["Slack" "Fourth"])
plot_v1=plot(sol1.dqsol.t,v_sim1',label=["Swing" "Fourth"],title="V in p.u.")
#plot_p = plot(sol.dqsol.t,p_sim',label=["Slack" "Fourth"])
plot_p1=plot(sol1.dqsol.t,p_sim1',label=["Swing" "Fourth"],title="P in p.u.")
plot_ω1=plot(sol1.dqsol.t,(ω_sim1')/(2*pi),label=["Swing" "Fourth"],title="Δf in Hz")
plot( plot_v1,plot_p1,plot_ω1;
        layout=(3,1),
        size = (1000, 500),
        lw=1,
        plot_title = "SwingFourthOrder",
        xlabel="t[s]")

savefig("SwingFourth.png")