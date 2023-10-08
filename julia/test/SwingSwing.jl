using .PowerDynamics
using OrderedCollections: OrderedDict
using Plots
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\components.jl")
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\LineParameters.jl")

buses = OrderedDict(
      "bus1" => SlackAlgebraic(U=1),
      "bus2" => SwingEqLVS(H=5., P=-1., D=1., Ω=50., Γ=200., V=1.0), 
      );  

buses1 = OrderedDict(
      "bus1" => SwingEqLVS(H=5., P=1., D=1., Ω=50., Γ=20., V=1.0), 
      "bus2" => SwingEqLVS(H=5., P=-1., D=1., Ω=50., Γ=100., V=1.0), 
      );  

length = 1.

branches=OrderedDict(
    "branch1" => StaticLine(from= "bus1", to = "bus2", Y = Z₀/(length*(R+X*im))),
     );

pg = PowerGrid(buses, branches)
pg1 = PowerGrid(buses1, branches)
  
operationpoint = find_operationpoint(pg,solve_powerflow=true,sol_method =:dynamic)

timespan= (0.0,10.)
fault = PowerPerturbation(node= "bus2",fault_power=-0.,tspan_fault = (2.,5.),var=:P)
#fault = ChangeInitialConditions(node="bus1", var=:v, f=Inc(-0.2))

sol= simulate(fault, pg, operationpoint, timespan)
sol1= simulate(fault, pg1, insert!(operationpoint[:],3,0.0), timespan)

v_sim=sol(sol.dqsol.t,:,:v)
v_sim1=sol1(sol1.dqsol.t,:,:v)
p_sim=sol(sol.dqsol.t,:,:p)
p_sim1=sol1(sol1.dqsol.t,:,:p)

ω_sim1=sol1(sol1.dqsol.t,:,:ω)  

plot_v= plot(sol.dqsol.t,v_sim',label=["Slack" "Swing"],title="V")
plot_v1=plot(sol1.dqsol.t,v_sim1',label=["Swing" "Swing"],title="V")
plot_p = plot(sol.dqsol.t,p_sim',label=["Slack" "Swing"],title="P")
plot_p1=plot(sol1.dqsol.t,p_sim1',label=["Swing" "Swing"],title="P")

plot_ω1=plot(sol1.dqsol.t,(ω_sim1')/(2*pi),label=["Swing" "Swing"],title="f")

plot( plot_v,plot_v1,plot_p,plot_p1,plot_ω1;
        layout=(3,2),
        size = (1000, 500),
        lw=3,
        plot_title = "SwingSwing",
        xlabel="t[s]")