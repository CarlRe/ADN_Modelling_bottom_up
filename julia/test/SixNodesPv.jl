using .PowerDynamics
using OrderedCollections: OrderedDict
using Plots
# for schleife 
buses = OrderedDict(
      "bus1" => SlackAlgebraic(U=1),
      "bus2" => SwingEqLVS(H=5., P=0.5, D=1., Ω=50., Γ=100., V=1.0),
      "bus3" => PQAlgebraic(P=-0.4,Q=-0.001),
      "bus4" => PQAlgebraic(P=-0.3,Q=-0.001),
      "bus5" => SwingEqLVS(H=2., P=0.2, D=0.5, Ω=50., Γ=25., V=1.0),
      "bus6" => SwingEqLVS(H=2., P=-0.2, D=0.5, Ω=50., Γ=25., V=1.0),
      "bus7" => SwingEqLVS(H=1., P=0.2, D=0.1, Ω=50., Γ=25., V=1.0),
      );  

buses1 = OrderedDict(
      "bus1" => PVAlgebraic(P=0.,V=1.),
      "bus2" => SwingEqLVS(H=5., P=0.51, D=1., Ω=50., Γ=100., V=1.0),
      "bus3" => PQAlgebraic(P=-0.4,Q=-0.001),
      "bus4" => PQAlgebraic(P=-0.3,Q=-0.001),
      "bus5" => SwingEqLVS(H=2., P=0.2, D=0.5, Ω=50., Γ=25., V=1.0),
      "bus6" => SwingEqLVS(H=1., P=-0.2, D=0.1, Ω=50., Γ=25., V=1.0),
      "bus7" => FourthOrderEq(T_d_dash=0.4, D=0.3, X_d=0.7979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=0.5, X_q_dash=0.646, P=0.2, H=1., E_f=1)
      );  

Z₀ = 16 # per unit Impedanz   [Ω] 
length1 = 1.5
length2 = 2.5
#Leitungsparameter 
R = 0.25; #Ω/km
X = 0.5;     

branches=OrderedDict(
    "branch1" => StaticLine(from= "bus1", to = "bus2", Y = Z₀/(length1*(R+X*im))),
    "branch2" => StaticLine(from= "bus2", to = "bus3", Y = Z₀/(length2*(R+X*im))),
    "branch3" => StaticLine(from= "bus2", to = "bus4", Y = Z₀/(length2*(R+X*im))),
    "branch4" => StaticLine(from= "bus4", to = "bus5", Y = Z₀/(length1*(R+X*im))),
    "branch6" => StaticLine(from= "bus4", to = "bus6", Y = Z₀/(length1*(R+X*im))),
    "branch7" => StaticLine(from= "bus2", to = "bus7", Y = Z₀/(length1*(R+X*im))),
     );

pg = PowerGrid(buses, branches)
pg1 = PowerGrid(buses1, branches)
  
operationpoint = find_operationpoint(pg,solve_powerflow=true,sol_method =:dynamic)

timespan= (0.0,10.)
fault = PowerPerturbation(node= "bus6",fault_power=-0.,tspan_fault = (5.,10.),var=:P)
fault1 = ChangeInitialConditions(node ="bus2", var=:v,f=Inc(-0.01))
fault2 = NodeParameterChange(node = "bus2", value = 0.99, tspan_fault = (5.,10.),var=:V)

sol= simulate(fault, pg, operationpoint, timespan)
sol1= simulate(fault2, pg1, insert!(operationpoint[:],size(operationpoint[:],1),0.0), timespan)

v_sim=sol(sol.dqsol.t,:,:v)
v_sim1=sol1(sol1.dqsol.t,:,:v)
i_sim=sol(sol.dqsol.t,:,:i)
i_sim1=sol1(sol1.dqsol.t,:,:i)
u_sim=sol(sol.dqsol.t,:,:u)
u_sim1=sol1(sol1.dqsol.t,:,:u)
p_sim=sol(sol.dqsol.t,:,:p)
p_sim1=sol1(sol1.dqsol.t,:,:p)
ω_sim =sol(sol.dqsol.t,["bus2", "bus5", "bus6", "bus7"],:ω)   
ω_sim1 =sol1(sol1.dqsol.t,["bus2" ,"bus5" ,"bus6", "bus7"],:ω)   

plot_v= plot(sol.dqsol.t,v_sim',label=["Slack" "Swing" "PQ" "PQ" "Swing" "SwingLoad" "Fourth"],title="V in p.u.")
plot_v1=plot(sol1.dqsol.t,v_sim1',label=["Pv" "Swing" "PQ" "PQ" "Swing" "SwingLoad" "Fourth"],title="V in p.u.")
plot_p = plot(sol.dqsol.t,p_sim',label=["Slack" "Swing" "PQ" "PQ" "Swing" "SwingLoad" "Fourth"],title="P in p.u.")
plot_p1=plot(sol1.dqsol.t,p_sim1',label=["Pv" "Swing" "PQ" "PQ" "Swing" "SwingLoad" "Fourth"],title="P in p.u.")

plot_ω = plot(sol.dqsol.t,((ω_sim'/(2*pi)).+50),label=["Swing" "Swing" "SwingLoad" "Fourth"],title="f in Hz")
plot_ω1=plot(sol1.dqsol.t,((ω_sim1'/(2*pi)).+50),label=["Swing" "Swing" "SwingLoad" "Fourth"])
plot( plot_v,plot_v1,plot_p,plot_p1,plot_ω,plot_ω1;
        layout=(3,2),
        size = (1000, 500),
        lw=0.5,
        plot_title = "SixBus",
        xlabel="t[s]")

