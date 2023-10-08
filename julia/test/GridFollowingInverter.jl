using .PowerDynamics
using OrderedCollections: OrderedDict
using Plots
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\components.jl")
  
buses = OrderedDict(
      "bus1" => SlackAlgebraic(U=1),
      #"bus2" => GridFollowingTecnalia(tau_u=0.1,omega_ini=50.,K_pomega=0.01,K_iomega=1,K_omega=0.5,K_v=0.001,omega_r=50.,V_r=0.1,P=0.5,Q_r=0.01),
      #"bus2" => ConverterPLL(P = 0.1,Q= 0.1,Ω = 50.,Kₚ= 0.01,Kᵢ= 0.1,K_ω= 0.01,Kᵥ=0.01,Vᵣ=1.)
      "bus2" => CSIPLL2(P = 0.1,Q= 0.1,Ω = 50.,Kₚ= 0.01,Kᵢ= 0.1,K_ω= 0.01,Vᵣ=1.)
      );  

Z₀ = 16 # per unit Impedanz   [Ω] 
length = 1.
#Leitungsparameter 
R = 0.1; #Ω/km
X = 0.2; 
branches=OrderedDict(
    "branch1" => StaticLine(from= "bus1", to = "bus2", Y = Z₀/(0.1*(R+X*im))),
     );

pg = PowerGrid(buses, branches)
states = [1.,0.,1.,0.,100.,0.]
state = State(pg, states)     
operationpoint = find_operationpoint(pg,sol_method =:dynamic)



timespan= (0.0,10.)
fault = PowerPerturbation(node= "bus2",fault_power=0.,tspan_fault = (2.,3.),var=:P)
#sol = simulate(fault, state , timespan)
sol= simulate(fault, pg, operationpoint, timespan)


v_sim=sol(sol.dqsol.t,:,:v)
v_real=sol(sol.dqsol.t,:,:u_r)
v_imag=sol(sol.dqsol.t,:,:u_i)
p_sim=sol(sol.dqsol.t,:,:p)
#p_sim1=sol1(sol1.dqsol.t,:,:p)

plot_v= plot(sol.dqsol.t,v_sim',xlims=(0.,10.),label= ["Slack" "GridFollow"])
plot_v_real= plot(sol.dqsol.t,v_real',xlims=(0.,10.),label= ["Slack" "GridFollow"])
plot_v_imag= plot(sol.dqsol.t,v_imag',xlims=(0.,10.),label= ["Slack" "GridFollow"])
plot_p = plot(sol.dqsol.t,p_sim',xlims=(0.,10.),label= ["Slack" "GridFollow"])

plot( plot_v,plot_v_real,plot_v_imag,plot_p;
        layout=(4,1),
        size = (1000, 500),
        lw=3,
        plot_title = "SlackGridFollowing",
        xlabel="t[s]")