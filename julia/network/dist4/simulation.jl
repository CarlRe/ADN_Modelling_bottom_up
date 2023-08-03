using PowerDynamics
using OrderedCollections: OrderedDict
using Plots
import PowerDynamics: dimension, symbolsof, construct_vertex

include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\components.jl")
include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\plot.jl")


powergrid =read_powergrid(joinpath(@__DIR__,"grid.json"), Json)
#operationpoint = find_operationpoint(powergrid,sol_method = :nlsolve)
timespan= (0.0,200.)
fault1 = ChangeInitialConditions(node = "bus0", var=:v, f = Inc(0.1))
fault2 = PowerPerturbation(node= "bus1",fault_power=0.15,tspan_fault = (120.,125.),var=:P)
fault3 = PowerPerturbation(node= "bus4",fault_power=-0.39,tspan_fault = (120.,150.),var=:P0)
#solution = simulate(fault3, powergrid, operationpoint, timespan)
#state_init = [1,0,0.1,1,0,1,0,0.1,1,0,1,0,1,0,1,0]
#state = State(powergrid,state_init)
states =  rand(systemsize(powergrid))
state = State(powergrid,states )
solution = simulate(fault3, state, timespan) 
#plot = create_plot(solution)
#v
#plot(solution.dqsol)


