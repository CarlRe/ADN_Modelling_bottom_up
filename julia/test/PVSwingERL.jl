using .PowerDynamics
using OrderedCollections: OrderedDict
using Plots
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\components.jl")


buses = OrderedDict(
      "bus1" => SlackAlgebraic(U=1),
      "bus2" => SwingEqLVS(H=5., P=0.2, D=1., Ω=50., Γ=20., V=1.0), 
      "bus3" => PQAlgebraic(P=-0.2, Q= 0.),
      );  

buses1 = OrderedDict(
      "bus1" => PVAlgebraic(P=0.,V=1.0),
      "bus2" => SwingEqLVS(H=5., P=0.2, D=1., Ω=50., Γ=20., V=1.0),
      "bus3" => ExponentialRecoveryLoad(P0=-0.2, Q0=-0., Nps=0.1, Npt=0.1, Nqs=0.1, Nqt=0.1, Tp=0.1, Tq=0.1, V0=1.)
      );  

Z₀ = 16 # per unit Impedanz   [Ω] 
length = 1.
#Leitungsparameter 
R = 0.2; #Ω/km
X = 0.35; 

branches=OrderedDict(
    "branch1" => StaticLine(from= "bus1", to = "bus2", Y = Z₀/(length*(R+X*im))),
    "branch2" => StaticLine(from= "bus1", to = "bus3", Y = Z₀/(length*(R+X*im)))
     );

pg = PowerGrid(buses, branches)
pg1 = PowerGrid(buses1, branches)
  
operationpoint = find_operationpoint(pg,solve_powerflow=true,sol_method =:dynamic)

timespan= (0.0,10.)
fault = PowerPerturbation(node= "bus3",fault_power=-0.1,tspan_fault = (2.,3.),var=:P)
fault1 = PowerPerturbation(node= "bus3",fault_power=-0.1,tspan_fault = (2.,3.),var=:P0) #ERL has P0 instead of P
sol= simulate(fault, pg, operationpoint, timespan)
sol1= simulate(fault1, pg1, insert!(insert!( operationpoint[:],8,0.0),9,0.0), timespan)

v_sim=sol(sol.dqsol.t,:,:v)
v_sim1=sol1(sol1.dqsol.t,:,:v)
p_sim=sol(sol.dqsol.t,:,:p)
p_sim1=sol1(sol1.dqsol.t,:,:p)
  

plot_v= plot(sol.dqsol.t,v_sim');plot!(sol1.dqsol.t,v_sim1',title="V")
plot_p = plot(sol.dqsol.t,p_sim');plot!(sol1.dqsol.t,p_sim1',title="P")

plot( plot_v,plot_p;
        layout=(1,2),
        size = (500, 250),
        lw=3,
        plot_title = "PvSwingERL",
        xlabel="t[s]")