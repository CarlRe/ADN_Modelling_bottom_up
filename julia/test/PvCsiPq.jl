using .PowerDynamics
using OrderedCollections: OrderedDict
using Plots
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\components.jl")


buses = OrderedDict(
      "bus1" => SlackAlgebraic(U=1),
      "bus2" => PQAlgebraic(P = -0.1,Q=0.),
      "bus3" => PQAlgebraic(P = 0.1,Q=0.),
      );  

buses1 = OrderedDict(
      "bus1" => PVAlgebraic(P=0.,V=1.0),
      "bus2" => PQAlgebraic(P = -0.1,Q=0.) ,
      "bus3" => CSIMinimal(I_r = 0.1)
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
  
operationpoint = find_operationpoint(pg,solve_powerflow=true,sol_method =:nlsolve)

timespan= (0.0,10.)
fault = PowerPerturbation(node= "bus2",fault_power=0. ,tspan_fault = (2.,3.),var=:P)
#state = State(pg1, operationpoint[:])  
#solution = simulate(fault, state , timespan), deleteat!(operationpoint[:],7) 

sol= simulate(fault, pg, operationpoint, timespan)
sol1= simulate(fault, pg1, operationpoint, timespan)

v_sim=sol(sol.dqsol.t,:,:v)
v_sim1=sol1(sol1.dqsol.t,:,:v)
p_sim=sol(sol.dqsol.t,:,:p)
p_sim1=sol1(sol1.dqsol.t,:,:p)
  

plot_v= plot(sol.dqsol.t,v_sim');plot!(sol1.dqsol.t,v_sim1',title="V")
plot_p = plot(sol.dqsol.t,p_sim');plot!(sol1.dqsol.t,p_sim1',title="P")

plot( plot_v,plot_p;
        layout=(1,2),
        size = (1000, 500),
        lw=3,
        plot_title = "PvCsiPQ",
        xlabel="t[s]")