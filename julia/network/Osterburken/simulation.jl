#Test für erste Versuche mit UW Osterburken Trafo 1

using PowerDynamics: read_powergrid, Json, Inc, find_operationpoint, ChangeInitialConditions, LineFailure, PowerPerturbation, simulate
using PowerDynamics: PVInverterWithFrequencyControl
using Plots
#using DifferentialEquations
#include(joinpath(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\InductionMotor.jl"))
powergrid = read_powergrid(joinpath(@__DIR__,"grid.json"), Json)
operationpoint = find_operationpoint(powergrid,sol_method = :nlsolve)
timespan= (0.0,10.0)

# simulating a Power perturbation at bus2
#fault1 = PowerPerturbation(node = "bus1",fault_power = -0.12, tspan_fault = (2.,6.),var=:P0)
fault1 = ChangeInitialConditions(node = "bus2", var=:ω, f = Inc(-0.1))
solution1 = simulate(fault1, powergrid, operationpoint, timespan)


#plot(solution1, "bus2", :iabs, legend = (1.5, 0.7), ylabel="ω [rad/s]",label = "generator_omega")
plot(solution1, "bus4", :v, legend = (1.1, 0.9), ylabel="V [p.u.]",label = "generator_v")


plot(solution1.dqsol)