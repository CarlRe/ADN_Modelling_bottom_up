using .PowerDynamics
using OrderedCollections: OrderedDict
using Plots
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\components.jl")


buses = OrderedDict(
      "bus1" => SlackAlgebraic(U=1.),
      "bus2" => VoltageDependentLoad(P=-0.1, Q=-0.1,U=1.,A=1.,B=0.3)
      );  

Z₀ = 16 # per unit Impedanz   [Ω] 
length = 1.
#Leitungsparameter 
R = 0.2; #Ω/km
X = 0.35; 

branches=OrderedDict(
    "branch1" => StaticLine(from= "bus1", to = "bus2", Y = Z₀/(length*(R+X*im))),
     );

pg = PowerGrid(buses, branches)
states = [1.0,0.,0.995,0.]
state = State(pg,states)
operationpoint = find_operationpoint(pg,sol_method =:dynamic)

timespan= (0.0,5.)

fault = NodeParameterChange(node = "bus1", value = 0.99, tspan_fault = (1.5,3.5),var=:U)  
sol= simulate(fault, pg, operationpoint, timespan)


v_sim=sol(sol.dqsol.t,"bus2",:v)

p_sim=sol(sol.dqsol.t,["bus2"],:p)

  

plot_v= plot(sol.dqsol.t,v_sim)
plot_p = plot(sol.dqsol.t,p_sim')

plot( plot_v,plot_p;
        layout=(1,2),
        size = (1000, 500),
        lw=3,
        plot_title = "ERL",
        xlabel="t[s]")