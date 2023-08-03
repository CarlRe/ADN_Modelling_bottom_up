using PowerDynamics: SlackAlgebraic,WindTurbineGenType4,PVInverterWithFrequencyControl, CSIMinimal, VoltageDependentLoad, ExponentialRecoveryLoad, FourthOrderEq, PQAlgebraic, PiModelLine, StaticLine, Transformer, PowerGrid, write_powergrid, Json
using OrderedCollections: OrderedDict
import PowerDynamics: PVInverterWithFrequencyControl, WindTurbineGenType4
#include(joinpath(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models", raw"LinearERL.jl"))
using PowerDynamics: read_powergrid, Json, Inc, find_operationpoint, ChangeInitialConditions, LineFailure, PowerPerturbation, simulate
using PowerDynamics: PVInverterWithFrequencyControl
using Plots
using DifferentialEquations

buses = OrderedDict(
    "bus0" => SlackAlgebraic(U=1),
    "bus1" => CSIMinimal(I_r=0.1));

#Leitungsparameter 
R = 0.3;
X = 0.3;
#Leitungslänge als variable 
length1 = 20;

#Leitungsbeläge
#Y = 1 / (length*R + 1*im*length*X)

branches=OrderedDict(
   
    "branch1" => StaticLine(from= "bus0", to = "bus1", Y = 1 / (length1*R + 1*im*length1*X)));
    
#

powergrid = PowerGrid(buses, branches)
#write_powergrid(powergrid, joinpath(@__DIR__,"gridCSI.json"), Json)

#powergrid = read_powergrid(joinpath(@__DIR__,"gridCSI.json"), Json)
operationpoint = find_operationpoint(powergrid,sol_method = :dynamic)
timespan= (0.0,10.)
fault1 = ChangeInitialConditions(node = "bus1", var=:v, f = Inc(0.01))
solution1 = simulate(fault1, powergrid, operationpoint, timespan)
plot(solution1, "bus1", :iabs, legend = (1.2, 0.0), ylabel="V [p.u.]",label = "generator_v")
#plot(solution1,"bus1",:v)