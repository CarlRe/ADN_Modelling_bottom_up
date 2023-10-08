using .PowerDynamics
using OrderedCollections: OrderedDict
using Plots
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\components.jl")
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\LineParameters.jl")  

buses = OrderedDict(
      "bus1" => SlackAlgebraic(U=1.),
      "bus2" => SwingEqLVS(H=0.01, P=0., D=0.1, Ω=50., Γ=20., V=1.0), 
      "bus3" => PQAlgebraic(P=0.2, Q= 0.1),
      "bus4" => PQAlgebraic(P=-0.1, Q= -0.1),
      "bus5" => SwingEqLVS(H=2., P=-0.2, D=0.5, Ω=50., Γ=10., V=1.0), 
      "bus6" => VoltageDependentLoad(P=-0.2,Q= -0.1, U=1.0,A=0.3,B=0.3),
      "bus7" => FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=0.1, H=5., E_f=1),
      

      );  

buses1 = OrderedDict(
      "bus1" => SlackAlgebraic(U=1.),
      "bus2" => PVAlgebraic(P=0.,V=1.),
      "bus3" => PQAlgebraic(P=0.2, Q= 0.1),
      "bus4" => PQAlgebraic(P=-0.1, Q= -0.1),
      "bus5" => SwingEqLVS(H=5., P=-0.1, D=0.5, Ω=50., Γ=20., V=1.0), 
      "bus6" => VoltageDependentLoad(P=-0.2,Q= -0.1, U=1.0,A=0.9,B=0.3),
      "bus7" => FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=0.1, X_q_dash=0.646, P=0.3, H=5., E_f=1),
     
      );  


length1 = 50.
length2 = 5.
length3 = 8.
branches=OrderedDict(
    "branch1" => StaticLine(from= "bus1", to = "bus2", Y = Y₀/length1),
    "branch2" => StaticLine(from= "bus2", to = "bus3", Y = Y₀/length2),
    "branch3" => StaticLine(from= "bus2", to = "bus4", Y = Y₀/length2),
    "branch4" => StaticLine(from= "bus4", to = "bus5", Y = Y₀/length3),
    "branch5" => StaticLine(from= "bus4", to = "bus6", Y = Y₀/length3),
    "branch6" => StaticLine(from= "bus2", to = "bus7", Y = Y₀/length2),
     );

pg = PowerGrid(buses, branches)
pg1 = PowerGrid(buses1, branches)
        
operationpoint = find_operationpoint(pg,solve_powerflow=true,sol_method =:dynamic)

timespan= (0.0,10.)
fault = PowerPerturbation(node= "bus6",fault_power=0.,tspan_fault = (2.,5.),var=:P)
fault1 = NodeParameterChange(node = "bus1", value = 0.9, tspan_fault = (1.,10.),var=:U)   
sol= simulate(fault1, pg, operationpoint, timespan)


v_sim=sol(sol.dqsol.t,:,:v)

p_sim=sol(sol.dqsol.t,:,:p)

q_sim=sol(sol.dqsol.t,:,:q)

ω_sim=sol(sol.dqsol.t,["bus2", "bus5", "bus7"],:ω)  
 
plot_v= plot(sol.dqsol.t,v_sim',label=["Slack" "Swing" "PQGen" "PQLoad" "SwingLoad" "VoltageDependentLoad" "Fourth"],title="V in p.u.")

plot_q = plot(sol.dqsol.t,q_sim',label=["Slack" "Swing" "PQGen" "PQLoad" "SwingLoad" "VoltageDependentLoad" "Fourth"],title="Q")

plot_p = plot(sol.dqsol.t,p_sim',label=["Slack" "Swing" "PQGen" "PQLoad" "SwingLoad" "VoltageDependentLoad" "Fourth"],title="P")

plot_ω=plot(sol.dqsol.t,(ω_sim')/(2*pi),label=[ "Swing" "SwingLoad" "Fourth"],title="Δf in Hz")

plot( plot_v,plot_ω,plot_p,plot_q;
        layout=(2,2),
        size = (1500, 1100),
        lw=1,
        plot_title = "7bus",
        xlabel="t[s]")

#savefig("7busWithSlack.png")