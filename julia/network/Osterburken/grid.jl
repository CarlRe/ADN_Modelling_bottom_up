using .PowerDynamics
using OrderedCollections: OrderedDict
using Plots
using LinearAlgebra
using JLD2
#import .PowerDynamics: dimension, symbolsof, construct_vertex

include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\components.jl")
#include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\plot.jl")
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\LineParameters.jl")
#include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\DataPrePro\EMT_DQ.jl")
V₀ = 1.24
V₁ = 1.
P_Load = -0.8667
Q_Load = 0.025
P_Gen = 0.08
buses = OrderedDict(
    "bus1" => SlackAlgebraic(U=V₀),
   # "bus2" => SwingEqLVS(H=10., P=0., D=2., Ω=50., Γ=10., V=V₁), 
    "bus3" => ExponentialRecoveryLoad(P0=P_Load, Q0=Q_Load, Nps=2., Npt=1., Nqs=3., Nqt=1., Tp=0.1, Tq=1., V0=V₀),
    #"bus4" => PQAlgebraic(P=0.1*P_Load, Q= 0.5*Q_Load),
    #"bus5" => SwingEqLVS(H=2., P=0.1*P_Load, D=2.5, Ω=49.5, Γ=10., V=V₁), 
    #"bus6" => VoltageDependentLoad(P=0.2*P_Load,Q= 0.5*Q_Load, U=V₁,A=0.3,B=0.3),
    #"bus7" => FourthOrderEq(T_d_dash=7.4, D=2, X_d=0.8979, X_q=0.646, Ω=50, X_d_dash=0.2995, T_q_dash=0.4, X_q_dash=0.646, P=P_Gen, H=5., E_f=V₁),
    );
#Leitungsparameter 
R = 0.2;
X = 0.5;
#Leitungslänge als variable 

#Leitungsbeläge
#Y = 1 / (length*R + 1*im*length*X)
length0 = 10.
length1 = 0.1
length2 = 3.
length3 = 5.
branches=OrderedDict(
  #  "branch1" => StaticLine(from= "bus1", to = "bus2", Y = 16/(length0*(R+im*X))),
    "branch2" => StaticLine(from= "bus1", to = "bus3", Y = Y₀/length1),
    #"branch3" => StaticLine(from= "bus1", to = "bus4", Y = Y₀/length1),
    #"branch4" => StaticLine(from= "bus4", to = "bus5", Y = Y₀/length2),
    #"branch5" => StaticLine(from= "bus4", to = "bus6", Y = Y₀/length3),
    #"branch6" => StaticLine(from= "bus1", to = "bus7", Y = Y₀/length2),
     );
#

powergrid = PowerGrid(buses, branches)
write_powergrid(powergrid, joinpath(@__DIR__,"grid.json"), Json)
#operationpoint = find_operationpoint(powergrid,solve_powerflow = true, sol_method = :dynamic)
operationpoint = find_operationpoint(powergrid, sol_method = :dynamic)
timespan= (0.0,4.)
fault = NodeParameterChange(node = "bus1", value = 1.23, tspan_fault = (0.1,4.),var=:U)  
solution = simulate(fault, powergrid, operationpoint, timespan)


np_powergrid = fault(powergrid)
regular = rhs(powergrid)
error = rhs(np_powergrid)

if regular.mass_matrix != error.mass_matrix || length(regular.syms) != length(error.syms)
    error("Change of MassMatrix or system size in abstract pertubation not supported!")
end

# wrap f and introduce parameter: if p=true no error, if p=false errorstate
_f = (dx, x, p, t) -> p ? regular(dx,x,nothing,t) : error(dx,x,nothing,t)

f = ODEFunction(_f, mass_matrix = regular.mass_matrix, syms = regular.syms)

problem = ODEProblem{true}(f, operationpoint, timespan, true)
save("sim.jld2", "problem", problem )
#display("MassMatrix is singular?: "*string(issingular(f.mass_matrix)))
display("Rang(M): "*string(rank(f.mass_matrix)))
display("Size(M): "*string(size(f.mass_matrix)))


v_sim=solution(solution.dqsol.t,["bus1", "bus3"],:v)
p_sim=solution(solution.dqsol.t,[ "bus3"],:p)
q_sim=solution(solution.dqsol.t,[ "bus3"],:q)
sim = [v_sim,p_sim,q_sim]
save("sim.jld2", "sim", sim )
#ω_sim=solution(solution.dqsol.t,["bus2"],:ω)

plot_v= plot(solution.dqsol.t,v_sim',label=["Slack" " Swing" "ERL"],title="V in p.u.")
plot_q = plot(solution.dqsol.t,q_sim',label=[ "ERL"],title="Q")
plot_p = plot(solution.dqsol.t,p_sim',label=[ "ERL"],title="P")
#plot_ω=plot(solution.dqsol.t,(ω_sim')/(2*pi),label="Swing",title="Δf in Hz")

plot( plot_v,plot_p,plot_q;
        layout=(3,1),
        size = (1500, 1100),
        lw=1,
        plot_title = "grid",
        xlabel="t[s]")

#savefig("ERL.png")