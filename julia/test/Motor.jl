using PowerDynamics
using OrderedCollections:OrderedDict
import PowerDynamics: dimension, symbolsof, construct_vertex, MassMatrix
using Plots
using Symbolics

include(raw"C:\Users\carlr\Documents\GitHub\ADN_Modelling_bottom_up\julia\models\components.jl")
 
buses = OrderedDict(
   #"bus1" => SwingEqLVS(H=100., P=0.01, D=10, Ω=50, Γ=20, V=1),
    "bus1" => SlackAlgebraic(U=1),
   # "bus2" => Motor(P = -1., ω0 = 50., Xs = 10, Xt = 0.5, T = 1, H = 1.2 )
   # "bus2" => MotorThird(P = -1., ω0 = 50)
   #"bus2" => MotorSingleCage(P=-0.1 ,ω0=1., H =1.5 )
   #"bus2" => ZIP(P0=-0.01,Q0=-0.1, A=0.1, B=0.1,C=0.8,D=2,E=-0.2,F=-0.8)
   "bus2" =>  PQAlgebraic(P=-0.011,Q=-0.01)
);

branches = OrderedDict( "branch1" =>StaticLine(from= "bus1", to = "bus2", Y = 1+ 1*im*1.));
#"branch2" =>StaticLine(from= "bus1", to = "bus3", Y = 0.5 + 1*im*0.5));

powergrid = PowerGrid(buses, branches)



#states = fill(0.5,systemsize(powergrid))
#states = rand(systemsize(powergrid))
#states = insert!(rand(systemsize(powergrid)-1),8,1.) 
#Stabil mit slack
#states= [0.024828497946497996,0.6202424671594247, 0.8242896080709682,0.45348987929787643,0.07534943289346785,0.3362188216770363,0.3983382833345466,0.4275644274982413,0.12110193238543121]
#states = [1.,0.,0.5,0.,0.5,0.1,0.1,50.,0.1]    
#state = State(powergrid, states)

timespan= (0.0,20.)
fault = PowerPerturbation(node= "bus2",fault_power=-0.02,tspan_fault = (10.,12.),var=:P0)
#fault = ChangeInitialConditions(node = "bus1", var=:v, f = Inc(0.1))
#
#solution = simulate(fault, state , timespan)
operationpoint1 = find_operationpoint(powergrid,solve_powerflow=true,sol_method = :nlsolve)
#solution = simulate(fault, powergrid,operationpoint, timespan)





function noError()
    while true
        try
          states = rand(systemsize(powergrid))
          state = State(powergrid, states)
          solution = simulate(fault, state , timespan)
          #mySim()  
          sleep(1.)
          break
        catch e
            println("nosucces")
        end
    end
    println("succes")
    display(states)
end

#=
function mySim()
    if rand() < 0.5
        throw(ErrorException("Random error occurred!"))
    end
    return "Simulation Result"
end
=#
