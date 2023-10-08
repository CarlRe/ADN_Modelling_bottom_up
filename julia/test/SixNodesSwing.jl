using .PowerDynamics
using OrderedCollections: OrderedDict
using Plots
# for schleife 
buses = OrderedDict(
      "bus1" => SlackAlgebraic(U=1),
      "bus2" => SwingEqLVS(H=5., P=0.1, D=0.5, Ω=50., Γ=20., V=1.0),
      "bus3" => PQAlgebraic(P=-0.2,Q=-0.001),
      "bus4" => PQAlgebraic(P=-0.3,Q=-0.001),
      "bus5" => SwingEqLVS(H=2., P=0.2, D=0.5, Ω=50., Γ=15., V=1.0),
      "bus6" => SwingEqLVS(H=2., P=-0.2, D=0.5, Ω=50., Γ=15., V=1.0),
      );  

buses1 = OrderedDict(
      "bus1" => SwingEqLVS(H=5., P=0.6, D=1., Ω=50., Γ=20., V=1.0),
      "bus2" => FourthOrderEq(T_d_dash=4.4, D=0.5, X_d=0.8979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=0.1, H=2.148, E_f=1),
      "bus3" => PQAlgebraic(P=-0.2,Q=-0.001),
      "bus4" => PQAlgebraic(P=-0.3,Q=-0.001),
      "bus5" => SwingEqLVS(H=5., P=0.2, D=0.5, Ω=50., Γ=5., V=1.0),
      "bus6" => SwingEqLVS(H=2., P=-0.2, D=0.5, Ω=50., Γ=5., V=1.0),
      );  

Z₀ = 16 # per unit Impedanz   [Ω] 
length = 5.
#Leitungsparameter 
R = 0.25; #Ω/km
X = 0.125;     

branches=OrderedDict(
    "branch1" => StaticLine(from= "bus1", to = "bus2", Y = Z₀/(length*(R+X*im))),
    "branch2" => StaticLine(from= "bus2", to = "bus3", Y = Z₀/(length*(R+X*im))),
    "branch3" => StaticLine(from= "bus2", to = "bus4", Y = Z₀/(length*(R+X*im))),
    "branch4" => StaticLine(from= "bus4", to = "bus5", Y = Z₀/(length*(R+X*im))),
    "branch6" => StaticLine(from= "bus4", to = "bus6", Y = Z₀/(length*(R+X*im))),
   
     );

pg = PowerGrid(buses, branches)
pg1 = PowerGrid(buses1, branches)
  
operationpoint = find_operationpoint(pg,solve_powerflow=true,sol_method =:dynamic)

timespan= (0.0,30.)
fault = PowerPerturbation(node= "bus3",fault_power=-0.,tspan_fault = (15.,16.),var=:P)
fault1 = ChangeInitialConditions(node ="bus1", var=:ω,f=Inc(-0.1))

sol= simulate(fault, pg, operationpoint, timespan)
sol1= simulate(fault, pg1, insert!(insert!( operationpoint[:],3,0.0),6,0.0), timespan)

v_sim=sol(sol.dqsol.t,:,:v)
v_sim1=sol1(sol1.dqsol.t,:,:v)
i_sim=sol(sol.dqsol.t,:,:i)
i_sim1=sol1(sol1.dqsol.t,:,:i)
u_sim=sol(sol.dqsol.t,:,:u)
u_sim1=sol1(sol1.dqsol.t,:,:u)
p_sim=sol(sol.dqsol.t,:,:p)
p_sim1=sol1(sol1.dqsol.t,:,:p)
ω_sim =sol(sol.dqsol.t,["bus2", "bus5", "bus6"],:ω)   
ω_sim1 =sol1(sol1.dqsol.t,["bus1" ,"bus2" ,"bus5" ,"bus6"],:ω) 
plot_v= plot(sol.dqsol.t,v_sim',label=["Slack" "Fourth" "PQ" "PQ" "Swing" "SwingLoad"],title="V")
plot_v1=plot(sol1.dqsol.t,v_sim1',label=["SwingSlack" "Fourth" "PQ" "PQ" "Swing" "SwingLoad"])
plot_p = plot(sol.dqsol.t,p_sim',label=["Slack" "Fourth" "PQ" "PQ" "Swing" "SwingLoad"],title="P")
plot_p1=plot(sol1.dqsol.t,p_sim1',label=["SwingSlack" "Fourth" "PQ" "PQ" "Swing" "SwingLoad"])
plot_ω = plot(sol.dqsol.t,ω_sim'/(2*pi),label=[ "Fourth" "Swing" "SwingLoad"],title="ω")
plot_ω1=plot(sol1.dqsol.t,ω_sim1'/(2*pi),label=["SwingSlack" "Fourth" "Swing" "SwingLoad"])
plot( plot_v,plot_v1,plot_p,plot_p1,plot_ω,plot_ω1;
        layout=(3,2),
        size = (1000, 500),
        lw=0.5,
        plot_title = "SixBus",
        xlabel="t[s]")