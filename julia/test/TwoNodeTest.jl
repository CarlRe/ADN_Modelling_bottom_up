using .PowerDynamics
using OrderedCollections: OrderedDict
using Plots
# for schleife 
loads = OrderedDict(
      "bus3" => PQAlgebraic(P=-0.2,Q=-0.),
      #"bus6" => VoltageDependentLoad(P=-0.2,Q=-0.001,U=1.,A=0.3,B=0.3),
     # "bus7" => SwingEqLVS(H=5., P=-0.2, D=0.5, Ω=50., Γ=5., V=1.0),
      );  
gens  = OrderedDict(
    "bus2" => SwingEqLVS(H=5., P=0.2, D=0.5, Ω=50., Γ=5., V=1.0),
    "bus3" => PQAlgebraic(P=0.2,Q=0.001),
    "bus4" => CSIMinimal(I_r=0.2),
    "bus5" => FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=0.2, H=5.148, E_f=1)
    );    
Z₀ = 16 # per unit Impedanz   [Ω] 
length = 1.
#Leitungsparameter 
R = 0.25; #Ω/km
X = 0.12;     
    


for value in loads 
  
    buses_slack = OrderedDict(
      "bus1" => SlackAlgebraic(U=1),
      "bus2" => SwingEqLVS(H=5., P=0.2, D=1., Ω=50., Γ=20., V=1.0),
      "bus3" => value[2],
    ); 
    buses_pv = OrderedDict(
        "bus1" => PVAlgebraic(P=0.,V=1.0),
        "bus2" => SwingEqLVS(H=5., P=0.2, D=1., Ω=50., Γ=20., V=1.0),
        "bus3" => value[2],  
    );
    
    branches=OrderedDict(
    "branch1" => StaticLine(from= "bus1", to = "bus2", Y = Z₀/(length*(R+X*im))),
    "branch2" => StaticLine(from= "bus2", to = "bus3", Y = Z₀/(length*(R+X*im)))
    );

    pg_slack = PowerGrid(buses_slack, branches)
    pg_pv = PowerGrid(buses_pv, branches)
    operationpoint = find_operationpoint(pg_slack,solve_powerflow=true,sol_method =:dynamic)
    timespan= (0.0,10.)
    fault = PowerPerturbation(node= "bus3",fault_power=-0.1,tspan_fault = (5.,6.),var=:P)
    sol_slack= simulate(fault, pg_slack, operationpoint, timespan)  
    sol_pv= simulate(fault, pg_pv, operationpoint , timespan)
    v_slack=(sol_slack(sol_slack.dqsol.t,:,:v))
    v_pv=sol_pv(sol_pv.dqsol.t,:,:v)
    p_slack=sol_slack(sol_slack.dqsol.t,:,:p)
    p_pv=sol_pv(sol_pv.dqsol.t,:,:p)
    plot_v= plot(sol_slack.dqsol.t,v_slack');plot!(sol_pv.dqsol.t,sol_pv')
    plot_p = plot(sol_slack.dqsol.t,p_slack');plot!(sol_pv.dqsol.t,sol_pv')
    plot( plot_v,plot_p;
        layout=(1,2),
        size = (1000, 500),
        lw=3,
       # plot_title = typeof(value),
        xlabel="t[s]")




end



