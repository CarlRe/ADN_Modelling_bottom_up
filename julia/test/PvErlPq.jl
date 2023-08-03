#using PowerDynamics: SlackAlgebraic, ExponentialRecoveryLoad, FourthOrderEq, StaticLine, PowerGrid, write_powergrid, Json
using OrderedCollections: OrderedDict
#using PowerDynamics: read_powergrid, Inc, find_operationpoint, ChangeInitialConditions, LineFailure, PowerPerturbation, simulate
using Plots
using PowerDynamics

buses = OrderedDict(
    "bus0" => SlackAlgebraic(U=1),
    "bus1" => FourthOrderEq(T_d_dash=6.1, D=2, X_d=1.05, X_q=0.98, 
    Ω=50, X_d_dash=0.185, T_q_dash=0.4,
     X_q_dash=0.36, P=1., H=6.54, E_f= 1),
    "bus2" => ExponentialRecoveryLoad(P0=-1, Q0=-0.1, Nps=0.6, Npt = 3,Nqs=2,Nqt=18,Tp=0.01,Tq=0.1,V0=1.));
#Leitungsparameter 
R = 0.3;
X = 0.3;
#Leitungsbeläge
#Y = 1 / (length*R + 1*im*length*X)
branches=OrderedDict(
    "branch1" => StaticLine(from= "bus0", to = "bus1", Y = 1 / (R + 1*im*X)),
    "branch2" => StaticLine(from= "bus0", to = "bus2", Y = 1 / (R + 1*im*X)));  
#

powergrid = PowerGrid(buses, branches)
operationpoint = find_operationpoint(powergrid,sol_method = :nlsolve)
timespan= (0.0,5.)
fault1 = ChangeInitialConditions(node = "bus1", var=:ω, f = Inc(0.01))
solution1 = simulate(fault1, powergrid, operationpoint, timespan)
#plot(solution1, "bus2", :iabs, legend = (1.2, 0.0), ylabel="V [p.u.]",label = "generator_v")
plot(solution1, "bus2", :v, legend = (1.2, 0.8), ylabel="V [p.u.]")