using DifferentialEquations, RecursiveArrayTools, Plots, DiffEqParamEstim
using Optimization, ForwardDiff, OptimizationOptimJL, OptimizationBBO
using JLD2
using RecursiveArrayTools # for VectorOfArray
#include(raw"julia\DataPrePro\EMT_DQ.jl")
#import .Data: v_mean, p_mean, q_mean, v_max, v_min, steps


#=
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\network\Osterburken\grid.jl")

problem = load("sim.jld2")["problem"]
sim_data = load("sim.jld2")["sim"]
t = collect(range(0,stop=4,length=10240))

cost_function = build_loss_objective(problem, Tsit5(), L2Loss(t,sim_data),
                                     Optimization.AutoForwardDiff(),
                                     maxiters=10000,verbose=false)

optprob = Optimization.OptimizationProblem(cost_function, [1.42])
optsol = solve(optprob, BFGS())
=#
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
operationpoint = find_operationpoint(powergrid, sol_method = :dynamic)
timespan= (0.0,4.)
fault = NodeParameterChange(node = "bus1", value = 1., tspan_fault = (0.1,4.),var=:U)  
solution = simulate(fault, powergrid, operationpoint, timespan)
v_sim=solution(solution.dqsol.t,["bus1", "bus3"],:v)
p_sim=solution(solution.dqsol.t,[ "bus3"],:p)
q_sim=solution(solution.dqsol.t,[ "bus3"],:q)

np_powergrid = fault(powergrid)
regular = rhs(powergrid)
error = rhs(np_powergrid)

# wrap f and introduce parameter: if p=true no error, if p=false errorstate
#_f = (dx, x, p, t) -> p ? regular(dx,x,nothing,t) : error(dx,x,nothing,t)
_f = (dx, x, p, t) -> p ? regular(dx,x,nothing,t) : error(dx,x,nothing,t)
f = ODEFunction(_f, mass_matrix = regular.mass_matrix, syms = regular.syms)

problem = ODEProblem{true}(f, operationpoint[:], timespan, true)

function errorState(integrator)
    integrator.p = false
    # reset the adaptive timestepping
    if integrator.opts.adaptive
        auto_dt_reset!(integrator)
        set_proposed_dt!(integrator, integrator.dt)
    end
end

function regularState(integrator)
    integrator.p = true
    # reset the adaptive timestepping
    if integrator.opts.adaptive
        auto_dt_reset!(integrator)
        set_proposed_dt!(integrator, integrator.dt)
    end
end

t1 = fault.tspan_fault[1]
t2 = fault.tspan_fault[2]

cb1 = PresetTimeCallback([t1], errorState)
cb2 = PresetTimeCallback([t2], regularState)

#Synthetic Data
sol_sim = solve(problem,Rodas4(); callback = CallbackSet(cb1, cb2))
sol_syn = solve(problem,Rodas4(); callback = CallbackSet(cb1, cb2))
t_syn = collect(range(0,stop=4.,length=40960))
randomized = VectorOfArray([(sol_syn(t_syn[i]) + .01randn(6)) for i in 1:length(t_syn)])
data_syn = convert(Array,randomized)
normalized = VectorOfArray([(sol_syn(t_syn[i]) ) for i in 1:length(t_syn)])
data_sim = convert(Array,normalized)


cost_function = build_loss_objective(problem, Rodas4(), L2Loss(t,data_syn),
                                Optimization.AutoForwardDiff())

optprob = Optimization.OptimizationProblem(cost_function,[1.,1.,1.,1.,1.,1.])
optsol = solve(optprob,BFGS())

println()
println(optsol.u)