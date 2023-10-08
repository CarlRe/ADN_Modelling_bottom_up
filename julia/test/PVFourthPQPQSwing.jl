using .PowerDynamics
using OrderedCollections: OrderedDict
using Plots
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\components.jl")
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\test\LineParameters.jl")  

buses = OrderedDict(
      "bus1" => SlackAlgebraic(U=1),
      "bus2" => SwingEqLVS(H=5., P=0.5, D=1., Ω=50., Γ=20., V=1.0), 
      "bus3" => PQAlgebraic(P=-0.3, Q= 0.),
      "bus4" => PQAlgebraic(P=-0.1, Q= 0.),
      "bus5" => SwingEqLVS(H=5., P=-0.1, D=1., Ω=50., Γ=20., V=1.0), 
      );  

buses1 = OrderedDict(
            "bus1" => SlackAlgebraic(U=1),
            "bus2" => FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=0.55, H=5.148, E_f=1),
            "bus3" => PQAlgebraic(P=-0.3, Q= 0.),
            "bus4" => PQAlgebraic(P=-0.1, Q= 0.),
            "bus5" => SwingEqLVS(H=5., P=-0.1, D=1., Ω=50., Γ=20., V=1.0),
            );  

buses2 = OrderedDict(
      "bus1" => PVAlgebraic(P=0.,V=1.0),
      "bus2" => FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=0.55, H=5.148, E_f=1),
      "bus3" => PQAlgebraic(P=-0.3, Q= 0.),
      "bus4" => PQAlgebraic(P=-0.1, Q= 0.),
      "bus5" => SwingEqLVS(H=5., P=-0.1, D=1., Ω=50., Γ=20., V=1.0),
      );  

length = 1.
branches=OrderedDict(
    "branch1" => StaticLine(from= "bus1", to = "bus2", Y = Y₀/length),
    "branch2" => StaticLine(from= "bus2", to = "bus3", Y = Y₀/length),
    "branch3" => StaticLine(from= "bus2", to = "bus4", Y = Y₀/length),
    "branch4" => StaticLine(from= "bus4", to = "bus5", Y = Y₀/length)
     );

pg = PowerGrid(buses, branches)
pg1 = PowerGrid(buses1, branches)
pg2 = PowerGrid(buses2, branches)  
#=   
    data, result = power_flow(pg)
        v = [result["solution"]["bus"][string(k)]["vm"] for k in 1:3]
        va = [result["solution"]["bus"][string(k)]["va"] for k in 1:3]
=#
        
operationpoint = find_operationpoint(pg,solve_powerflow=true,sol_method =:dynamic)

#state = State(pg1, operationpoint[:])  

timespan= (0.0,10.)
fault = PowerPerturbation(node= "bus3",fault_power=-0.,tspan_fault = (2.,3.),var=:P)
    
sol= simulate(fault, pg1, insert!(operationpoint[:],5 ,0.0), timespan)
sol1= simulate(fault, pg1, insert!(operationpoint[:],5 ,0.0), timespan)

v_sim=sol(sol.dqsol.t,:,:v)
v_sim1=sol1(sol1.dqsol.t,:,:v)
p_sim=sol(sol.dqsol.t,:,:p)
p_sim1=sol1(sol1.dqsol.t,:,:p)
ω_sim=sol(sol.dqsol.t,["bus2"],:ω)  
ω_sim1=sol1(sol1.dqsol.t,["bus2"],:ω)  
plot_v= plot(sol.dqsol.t,v_sim',label=["Slack" "Fourth" "PQ"])
plot_v1=plot(sol1.dqsol.t,v_sim1',label=["PV" "Fourth" "PQ"],title="V")
plot_p = plot(sol.dqsol.t,p_sim',label=["Slack" "Fourth" "PQ"])
plot_p1=plot(sol1.dqsol.t,p_sim1',label=["PV" "Fourth" "PQ"],title="P")
plot_ω=plot(sol.dqsol.t,(ω_sim')/(2*pi),label=[ "Fourth"],title="f")
plot_ω1=plot(sol1.dqsol.t,(ω_sim1')/(2*pi),label=["Fourth"],title="f")
plot( plot_v,plot_v1,plot_p,plot_p1,plot_ω,plot_ω1;
        layout=(3,2),
        size = (1500, 750),
        lw=3,
        plot_title = "PvFourthPq",
        xlabel="t[s]")