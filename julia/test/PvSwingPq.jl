using .PowerDynamics
using OrderedCollections: OrderedDict
using Plots
#include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\PowerDynamics.jl\src\PowerDynamics_New.jl")
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\components.jl")

function create_plot(sol)
  #generator_indices = findall(bus -> typeof(bus) == SwingEqLVS,powergrid.nodes)
  #labels = reshape(generator_indices,(1,length(generator_indices)))

  pl_v = plot(sol, "bus2", :v, legend = (0.8, 0.7), ylabel="V [p.u.]")
  pl_p = plot(sol, "bus2", :p,ylims=(-10.,10.) ,legend = (0.8, 0.7), ylabel="p [p.u.]" )
  pl_q = plot(sol, "bus2", :q, ylims=(-10.,10.), legend = (0.8, 0.7), ylabel="q [p.u.]")
  pl_ω = plot(sol, "bus2", :ω,ylims=(-5.,5.), legend = (0.8, 0.7), ylabel="ω [rad/s]")
  
  pl = plot( pl_v, pl_ω, pl_p, pl_q;
          layout=(2,2),
          size = (1000, 500),
          lw=3,
          xlabel="t[s]")
end


#=
buses = OrderedDict(
    "bus1" => SlackAlgebraic(U=1),
    "bus2" => SwingEqLVS(H=5., P=0.1, D=1., Ω=1., Γ=20., V=1.0),
    "bus3" => FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=1., X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=0.7, H=5.148, E_f=1),
    "bus4" => PQAlgebraic(P=-0.8, Q= -0.1),
    );

   =# buses = OrderedDict(
      "bus1" => SlackAlgebraic(U=1),
      #"bus1" => PVAlgebraic(P=0.,V=1.0),
      "bus2" => SwingEqLVS(H=5., P=0.2, D=1., Ω=50., Γ=20., V=1.0),
     
     #"bus3" => FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=50., X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=0.1, H=5.148, E_f=1),
      "bus3" => PQAlgebraic(P=-0.2, Q= 0.),
      );  

      buses1 = OrderedDict(
      #"bus1" => SlackAlgebraic(U=1),
      "bus1" => PVAlgebraic(P=0.,V=1.0),
      "bus2" => SwingEqLVS(H=5., P=0.2, D=1., Ω=50., Γ=20., V=1.0),
     
     #"bus3" => FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=50., X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=0.1, H=5.148, E_f=1),
      "bus3" => PQAlgebraic(P=-0.2, Q= 0.),
      );  

  #=    
buses = OrderedDict(
      "bus1" => SwingEqLVS(H=5., P=0.81, D=1., Ω=1., Γ=10., V=1.0),
      "bus2" =>PQAlgebraic(P=-0.8, Q= -0.0),    
    );
   =#  

    Z₀ = 16 # per unit Impedanz   [Ω] 
length = 1.
#Leitungsparameter 
R = 0.2; #Ω/km
X = 0.35; 
branches=OrderedDict(

    "branch1" => StaticLine(from= "bus1", to = "bus2", Y = Z₀/(length*(R+X*im))),
    "branch2" => StaticLine(from= "bus2", to = "bus3", Y = Z₀/(length*(R+X*im)))
   #"branch3" => StaticLine(from= "bus2", to = "bus4", Y = Z₀/(length*(R+X*im))),
     #"branch4" => StaticLine(from= "bus1", to = "bus3", Y = Z₀/(length*(0.3+0.3*im)))
     );
    pg = PowerGrid(buses, branches)
    pg1 = PowerGrid(buses1, branches)
  
   
    data, result = power_flow(pg)
        v = [result["solution"]["bus"][string(k)]["vm"] for k in 1:3]
        va = [result["solution"]["bus"][string(k)]["va"] for k in 1:3]
  operationpoint = find_operationpoint(pg,solve_powerflow=true,sol_method =:dynamic)

    state = State(pg1, operationpoint[:])  
    timespan= (0.0,10.)
    fault = PowerPerturbation(node= "bus3",fault_power=-0.0,tspan_fault = (.5,5.),var=:P)
    solution = simulate(fault, state , timespan)

    sol= simulate(fault, pg, operationpoint, timespan)
    sol1= simulate(fault, pg1, operationpoint, timespan)

    pl_v = plot(sol, "bus2", :v, legend = (0.8, 0.7), ylabel="V [p.u.]")
    pl_p = plot(sol, "bus2", :p,ylims=(-10.,10.) ,legend = (0.8, 0.7), ylabel="p [p.u.]" )
    pl_q = plot(sol, "bus2", :q, ylims=(-10.,10.), legend = (0.8, 0.7), ylabel="q [p.u.]")

