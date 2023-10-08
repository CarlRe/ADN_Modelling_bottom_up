using PowerDynamics
using JuMP
using OrderedCollections: OrderedDict
using EasyFit
using .Data
# Definition of simple Network 
#Slack, ERL with Admittanz 

t = Data.t # time 
Vabs = Data.v_abs  #xData
P_meas = Data.p_mean #yData

#smoothed p_meas 

p_mean = movavg(P_meas,1000)
p_mean = p_mean.x.*(-1)

v_mean = movavg(Vabs,1000)
v_mean = v_mean.x

V₀ = 1.24
V₁ = 1.2
P_Load = -0.01
Q_Load = -0.02

# parameter vektor 
# 6 parameter für ERL 
p = [2.;1.;3.;1.;0.5;0.5]


buses = OrderedDict(
    "bus1" => SlackAlgebraic(U=V₀),
    "bus2" =>  ExponentialRecoveryLoad(P0=P_Load, Q0=Q_Load, Nps=p[1], Npt=p[2], Nqs=p[3], Nqt=p[4], Tp=p[5], Tq=p[6], V0=V₀),
)

branches = OrderedDict(
    "branch1" => StaticLine(from = "bus1",to = "bus2",Y=0.5+0.5*im)
)

powergrid = PowerGrid(buses,branches)
operationpoint = find_operationpoint(powergrid, sol_method = :dynamic)
timespan= (0.0,4.)
fault = NodeParameterChange(node = "bus1", value = V₁, tspan_fault = (0.1,4.),var=:U)  
solution = simulate(fault, powergrid, operationpoint, timespan)

#activePower = solution[]
t = solution.dqsol.t
p_sim = solution(t,:,:p)[2,:]  #Leistungsverlauf bus2
q_sim = solution(t,:,:q)[2,:]
powerplot = plot(t,p_sim,label = "Wirkleistung in p.u.");plot!(t,q_sim,label = "Blindleistung in p.u.")

#Vektor interpolieren für gleiche Länge
function dotting(smallVector, longVector)
    dottedVector = smallVector
    length_old = length(smallVector)
    length_new = length(longVector)
    missingElements = length_new-length_old 
    dotsPerSpace = Int(round(missingElements/(length_old-1)))
    println("dotsPerSpace: "*string(dotsPerSpace))
    for i=1:((length_old-1))
        position = i+dotsPerSpace*(i-1)
         mean = (smallVector[position]+smallVector[position+1])/2
         println("Schritt: "*string(i))
         println("Mean: "*string(mean))
         println("Länge: "*string(length(dottedVector)))
         for j=1:dotsPerSpace
            insert!(dottedVector,position+j,mean)
         end
    end
    return dottedVector 
end 
p_sim_new = dotting(p_sim,p_mean)
p_sim_smooth = movavg(p_sim_new,1516)
p_sim_smooth = p_sim_smooth.x
powerplot_new =  plot(range(0.,4.,40960),p_sim_new,label = "Wirkleistung in p.u.[dotting]");plot!(range(0.,4.,40960),p_sim_smooth,label = "Wirkleistung in p.u.[smoothed]")






#Bis hierhin reine Simulation

function optimization(sim,meas)
    model = Model(with_optimizer(Ipopt.Optimizer));

    @variable(model, Nps)
    @variable(model, Npt)
    @variable(model, Nqs) 
    @variable(model, Nqt)
    @variable(model, Tp)
    @variable(model, Tq) 

    set_start_value(Nps, 2.)
    set_start_value(Npt, 1.)
    set_start_value(Nqs, 3.)
    set_start_value(Nqt, 1.)
    set_start_value(Tp, 0.5)
    set_start_value(Tq, 0.5)
    
    @objective(model, Min, sum((sim - meas).^2 ))
    
    optimize!(model)
    
    param_fit = [value(Nps),value(Npt),value(Nqs), value(Nqt),value(Tp),value(Tq)]

    return param_fit
end

function sim(params)
    p = params    
    buses = OrderedDict(
        "bus1" => SlackAlgebraic(U=V₀),
        "bus2" =>  ExponentialRecoveryLoad(P0=P_Load, Q0=Q_Load, Nps=p[1], Npt=p[2], Nqs=p[3], Nqt=p[4], Tp=p[5], Tq=p[6], V0=V₀),
    )

    branches = OrderedDict(
        "branch1" => StaticLine(from = "bus1",to = "bus2",Y=0.5+0.5*im)
    )

    powergrid = PowerGrid(buses,branches)
    operationpoint = find_operationpoint(powergrid, sol_method = :dynamic)
    timespan= (0.0,4.)
    fault = NodeParameterChange(node = "bus1", value = V₁, tspan_fault = (0.1,4.),var=:U)  
    solution = simulate(fault, powergrid, operationpoint, timespan)

    p_sim = solution(t,:,:p)[2,:]

    p_sim_new = dotting(p_sim,p_mean)
    p_sim_smooth = movavg(p_sim_new,1516)
    p_sim_smooth = p_sim_smooth.x
    return p_sim_smooth
end

function runOpt(initial,meas)
    sol = sim(initial)
    param = optimization(sol,meas)
    sol = sim(param)
    
    return param
end

parameter = runOpt(p,p_mean)