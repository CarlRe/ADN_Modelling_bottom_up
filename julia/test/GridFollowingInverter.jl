using PowerDynamics
using OrderedCollections: OrderedDict
using Plots
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\components.jl")
#=
buses = OrderedDict(
    "bus1" => SlackAlgebraic(U=1),
    "bus2" => SwingEqLVS(H=5., P=0.1, D=1., Ω=1., Γ=20., V=1.0),
    "bus3" => FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=1., X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=0.1, H=5.148, E_f=1),
    "bus4" => PQAlgebraic(P=-0.2, Q= -0.01),
    #"bus5" => GridFollowingTecnalia(tau_u=1.9998,omega_ini=50,K_pomega=0.001,K_iomega=0.02,K_omega=40000,K_v=0.8,omega_r=50.,V_r=1,P=0.2,Q_r=0.01)
    );

#GridFollowingTecnalia(tau_u=1.9998,omega_ini=50,K_pomega=0.001,K_iomega=0.02,K_omega=40000,K_v=0.8,omega_r=50.,V_r=1,P=0.1,Q_r=0.2)
    buses = OrderedDict(
      "bus1" => SwingEqLVS(H=5., P=0.1, D=1., Ω=50., Γ=20., V=1.0),
      "bus2" => FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=50., X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=0.1, H=5.148, E_f=1),
      "bus3" => PQAlgebraic(P=-0.3, Q= -0.01),
      "bus4" => GridFollowingTecnalia(tau_u=1.9998,omega_ini=50,K_pomega=0.001,K_iomega=0.02,K_omega=1/40000,K_v=0.8,omega_r=50.,V_r=1,P=0.1,Q_r=0.01)
      );  

  =#  
buses = OrderedDict(
      "bus1" => SlackAlgebraic(U=1),
      "bus2" => GridFollowingTecnalia(tau_u=1.9998,omega_ini=50,K_pomega=0.001,K_iomega=0.02,K_omega=1/40000,K_v=0.8,omega_r=50.,V_r=1,P=0.1,Q_r=0.01),
      "bus3" => Motor(P = -0.1, ω0 = 50, Xs = 5, Xt = 0.5,T = 10, H = 1 )  
    );
    

    Z₀ = 16 # per unit Impedanz   [Ω] 
length = 1
#Leitungsparameter 
R = 0.2; #Ω/km
X = 0.35; 
branches=OrderedDict(

    "branch1" => StaticLine(from= "bus1", to = "bus2", Y = Z₀/(length*(R+X*im))),
    "branch2" => StaticLine(from= "bus2", to = "bus3", Y = Z₀/(length*(R+X*im))),
   #"branch3" => StaticLine(from= "bus2", to = "bus4", Y = Z₀/(length*(R+X*im))),
     #"branch4" => StaticLine(from= "bus2", to = "bus5", Y = Z₀/(length*(0.3+0.3*im)))
     );
    pg2 = PowerGrid(buses, branches)
  # operationpoint = find_operationpoint(pg1,sol_method =:nlsolve)
  #operationpoint9 = find_operationpoint(pg1,solve_powerflow=true,sol_method =:nlsolve)
  _, result9 = power_flow(pg2)
 #=
    _, result8 = power_flow(pg1)
        v = [result["solution"]["bus"][string(k)]["vm"] for k in 1:3]
        va = [result["solution"]["bus"][string(k)]["va"] for k in 1:3]
=#
#=
states = [0.9999948357533613, -0.0032137953722651304, 1.5129565341213817e-11, 1.008221872929932, 0.012320413854144787, 0.012219334833354257, 0.0, 0.9873320872360581, -0.019423183599188334,0.98,0.01,0.0,0.0,0.0,0.0]
    state = State(pg1, states)  
    timespan= (0.0,200.)
    fault = PowerPerturbation(node= "bus3",fault_power=-0.0,tspan_fault = (10.5,11.),var=:P)
    solution = simulate(fault, state , timespan)

=#