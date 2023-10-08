using .PowerDynamics
using OrderedCollections: OrderedDict
using Plots
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\components.jl")
  
buses = OrderedDict(
      "bus1" => SlackAlgebraic(U=1.),
      #"bus2" => Motor(P = 0.4 ,ω0= 50, Xs= 1., Xt=0.5, T= 5., H= 5. ) , 
      #"bus2" => MotorFirst(P=0.1, V₀=1. , Tₚ=2. , αₛ=0.3, αᵣ=0.1) 
     # "bus2" => MotorThird(P=-0.1, ω0 = 50.)
      "bus2" => MotorSlipModel1(P=-0.1, X_s=0.1 , X_r=0.2 ,R_r=0.1, E=1.)
      #"bus2" => MotorSingleCageNew(τ_m0=0.1, ω0 = 50.)
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
operationpoint = find_operationpoint(pg,sol_method =:dynamic)
states = [1.,0.,1.,0.,0.,0.1]
state = State(pg, states) 
timespan= (0.0,10.)
fault = NodeParameterChange(node = "bus1", value = 0.99, tspan_fault = (8.,10.),var=:U)  
sol= simulate(fault, pg, state, timespan)


v_sim=sol(sol.dqsol.t,:,:v)

p_sim=sol(sol.dqsol.t,:,:p)


plot_v= plot(sol.dqsol.t,v_sim',label= ["Slack" "Motor"])
plot_p = plot(sol.dqsol.t,p_sim',label= ["Slack" "Motor"])

plot( plot_v,plot_p;
        layout=(1,2),
        size = (1000, 500),
        lw=3,
        plot_title = "SlackMotor",
        xlabel="t[s]")