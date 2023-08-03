using PowerDynamics
using OrderedCollections: OrderedDict
using Plots
import PowerDynamics: dimension, symbolsof, construct_vertex

include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\components.jl")
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\plot.jl")

buses = OrderedDict(
    # "bus1" => SlackAlgebraic(U=1),
    "bus1" => SwingEqLVS(H=5., P=0.1, D=0.1, Ω=50., Γ=20., V=1.), #Ersatz für Slack
   
    # Erzeugung  
    "bus2" =>FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=0.1, H=5.148, E_f=1), #Synchron Einspeisung
   #"bus3" => GridFollowingTecnalia(tau_u=1.9998,omega_ini=50,K_pomega=0.001,K_iomega=0.02,K_omega=40000,K_v=0.8,omega_r=50.,V_r=1,P=0.1,Q_r=0.2),
    #Lasten
    #"bus3" => Motor(P = -0.1, ω0 = 50, Xs = 5, Xt = 0.5,T = 10, H = 1 ), #Motor Last 
    #"bus3" => ZIP(P0=-0.1,Q0=-0.01, A=0.1, B=0.1,C=0.8,D=2,E=-0.2,F=-0.8),
     "bus3" => PQAlgebraic(P=-0.2,Q=0.),
     
    );
#Leitungsparameter 
R = 0.2; #Ω/km
X = 0.35; 
# p.u. Impedanz
Z₀=16 #Ω
#Leitungslänge als variable 
length1 = 1;
length2 = 1;
length3 = 1;
length4 = 1;
length5 = 1;
length6 = 1;
length7 = 1;
#Leitungsbeläge
#Y_pu = 1/Z_pu = 1/(Z/Z₀) = Z₀ / (length*R + 1*im*length*X)

branches=OrderedDict(
   
  # "branch1" => StaticLine(from= "bus1", to = "bus2", Y = Z₀ / (length1*R + 1*im*length1*X)),
   "branch1" => StaticLine(from= "bus1", to = "bus2", Y =  Z₀ / (length2*R + 1*im*length2*X)),
   "branch2" => StaticLine(from= "bus1", to = "bus3", Y =  Z₀/ (length3*R + 1*im*length3*X)),
   # "branch4" => StaticLine(from= "bus3", to = "bus4", Y =  Z₀ / (length4*R +1*im*length4*X))
    );
#
powergrid = PowerGrid(buses, branches)
#write_powergrid(powergrid, joinpath(@__DIR__,"grid.json"), Json)

#=
_, result = power_flow(powergrid)
v = [result["solution"]["bus"][string(k)]["vm"] for k in 1:3]
va = [result["solution"]["bus"][string(k)]["va"] for k in 1:3]
=#
_, result = power_flow(powergrid)
#ic = ( v .* exp.(1im .* va))
operationpoint = find_operationpoint(powergrid, solve_powerflow = true, sol_method = :dynamic)


#=
timespan= (0.0,200.)
fault1 = ChangeInitialConditions(node = "bus0", var=:v, f = Inc(0.1))
fault2 = PowerPerturbation(node= "bus1",fault_power=0.15,tspan_fault = (120.,125.),var=:P)
fault3 = PowerPerturbation(node= "bus4",fault_power=-0.39,tspan_fault = (120.,150.),var=:P0)
solution = simulate(fault3, powergrid, operationpoint, timespan)
=#


#state_init=fill(0.8,(systemsize(powergrid)))
#state = State(powergrid,state_init)
#states =  rand(systemsize(powergrid))
#state = State(powergrid,states )
#solution = simulate(fault3, state, timespan) 


