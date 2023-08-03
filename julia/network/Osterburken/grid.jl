using PowerDynamics
using OrderedCollections: OrderedDict
using Plots
import PowerDynamics: dimension, symbolsof, construct_vertex

include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\components.jl")
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\plot.jl")

buses = OrderedDict(
    "bus0" => SwingEqLVS(H=1., P=1.0, D=2., Ω=50, Γ=20, V=1), #Ersatz für Slack
    #"bus0" => SlackAlgebraic(U=1),
    # Erzeugung 
    "bus1" => PQAlgebraic(P=0.1,Q=0.1),#Wind Einspeisung 
    "bus2" => FourthOrderEq(T_d_dash=1.1, D=2, X_d=1.05, X_q=0.98,  Ω=50, X_d_dash=0.185, T_q_dash=0.4,  X_q_dash=0.36, P=0.1, H=6.54, E_f= 1), #Synchron Einspeisung
    "bus3" => GridFollowingTecnalia(tau_u=1.9998,omega_ini=50,K_pomega=0.001,K_iomega=0.02,K_omega=40000,K_v=0.8,omega_r=50,V_r=1,P=0.2,Q_r=0.2),
    
    #Lasten
    "bus4" => Motor(P = -0.2, ω0 = 50, Xs = 5, Xt = 0.5, T = 10, H = 1 ), #Motor Last 
    "bus5" => ZIP(P0=-0.2,Q0=-0.1, A=0.1, B=0.1,C=0.8,D=2,E=-0.2,F=-0.8),
    "bus6" => PQAlgebraic(P=-0.1, Q = -0.2) 
    );
#Leitungsparameter 
R = 0.3;
X = 0.3;
#Leitungslänge als variable 
length1 = 1;
length2 = 1;
length3 = 1;
length4 = 1;
length5 = 1;
length6 = 1;
length7 = 1;
#Leitungsbeläge
#Y = 1 / (length*R + 1*im*length*X)

branches=OrderedDict(
   
    "branch1" => StaticLine(from= "bus0", to = "bus1", Y = 1 / (length1*R + 1*im*length1*X)),
    "branch2" => StaticLine(from= "bus0", to = "bus2", Y =  1 / (length2*R + 1*im*length2*X)),
   # "branch3" => StaticLine(from= "bus0", to = "bus3", Y =  1 / (length3*R + 1*im*length3*X)),
  #  "branch4" => StaticLine(from= "bus0", to = "bus4", Y =  1 / (length4*R +1*im*length4*X)),
    "branch5" => StaticLine(from= "bus0", to = "bus5", Y =  1 / (length5*R+ 1*im*length5*X)),
    "branch6" => StaticLine(from = "bus0", to = "bus6",Y =  1 / (length6*R + 1*im*length6*X))
    #"branch7" => StaticLine(from= "bus0", to ="bus7", Y =  1 / (length7*R+1*im*length7*X))
    );
#

powergrid = PowerGrid(buses, branches)
write_powergrid(powergrid, joinpath(@__DIR__,"grid.json"), Json)
#operationpoint = find_operationpoint(powergrid,sol_method = :nlsolve)
timespan= (0.0,200.)
fault1 = ChangeInitialConditions(node = "bus2", var=:v, f = Inc(0.1))
fault2 = PowerPerturbation(node= "bus1",fault_power=0.15,tspan_fault = (120.,125.),var=:P)
fault3 = PowerPerturbation(node= "bus5",fault_power=-0.19,tspan_fault = (2.,5.),var=:P0)
#solution = simulate(fault3, powergrid, operationpoint, timespan)
#state_init = [1,0,0.1,1,0,1,0,0.1,1,0,1,0,1,0,1,0]
#state = State(powergrid,state_init)
states =  rand(systemsize(powergrid))
state = State(powergrid,states )
solution = simulate(fault2, state, timespan) 
#plot = create_plot(solution)
#v
#plot(solution.dqsol)


